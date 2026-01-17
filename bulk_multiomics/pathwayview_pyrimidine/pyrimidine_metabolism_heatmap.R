#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(stringr)
  library(KEGGREST)
  library(httr)
  library(jsonlite)
  library(ComplexHeatmap)
  library(circlize)
  library(xml2)
  library(png)
  library(grid)
})

options(stringsAsFactors = FALSE, scipen = 999)

# ======================
# 0) Paths and settings
# ======================
base_dir <- "/data1/dyh/GDF11/difftable"
out_dir  <- file.path(base_dir, "pathveiw")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Genotypes / groups
GENOTYPES <- c("sh", "C", "OE") # column order

# Significance thresholds (stars)
P_ADJ_CUTOFF <- 0.05  # transcript/protein: adjusted P (padj / p.adj)
P_MET_CUTOFF <- 0.05  # metabolite: p-value
VIP_CUTOFF   <- 1.0   # metabolite VIP cutoff (set to NA to disable)

# Heatmap style (red-green as requested; negative=green, positive=red)
COL_DOWN <- "#1A9850" # green  (negative log2FC)
COL_MID  <- "white"
COL_UP   <- "#D73027" # red    (positive log2FC)

# If NULL, auto by 95% quantile of |log2FC| (min 2.0)
LOG2FC_CAP_ABS <- NULL

# ======================
# 1) KEGG pathway members
# ======================
# KEGG REST can be flaky; prefer local cached KGML if present
pathway_id <- "path:cfa00240" # Canis familiaris pyrimidine metabolism
kgml_path <- file.path(out_dir, "cfa00240.kgml")
png_path  <- file.path(out_dir, "cfa00240.png")
txt_path  <- file.path(out_dir, "cfa00240.txt")

# Keep KEGG original layout: use KEGG high-resolution (2x) pathway background,
# and draw heatmaps on top (vector overlays). This avoids the "redraw" layout drift.
USE_KEGG_BASE_2X <- TRUE
KEGG_BASE_PNG_2X <- file.path(out_dir, "map00240@2x.png")
KEGG_BASE_URL_2X <- "https://www.kegg.jp/kegg/pathway/map/map00240@2x.png"

parse_kegg_flat <- function(txt_file) {
  lines <- readLines(txt_file, warn = FALSE)

  gene_start <- which(startsWith(lines, "GENE"))[1]
  comp_start <- which(startsWith(lines, "COMPOUND"))[1]
  if (is.na(gene_start) || is.na(comp_start) || comp_start <= gene_start) {
    return(list(gene_map = data.table(), compounds = character()))
  }

  gene_lines <- lines[gene_start:(comp_start - 1)]
  gene_lines[1] <- sub("^GENE\\s+", "", gene_lines[1])
  gene_lines <- gene_lines[str_detect(gene_lines, "^\\s*\\d+\\s+")]
  gene_dt <- data.table(raw = gene_lines)
  gene_dt[, entrez := str_match(raw, "^\\s*(\\d+)")[, 2]]
  gene_dt[, gene := str_match(raw, "^\\s*\\d+\\s+([^;\\s]+)\\s*;")[, 2]]
  gene_dt <- gene_dt[!is.na(entrez) & !is.na(gene) & gene != ""]
  gene_dt <- unique(gene_dt[, .(entrez = as.character(entrez), gene = as.character(gene))])

  comp_lines <- lines[comp_start:length(lines)]
  comp_lines[1] <- sub("^COMPOUND\\s+", "", comp_lines[1])
  comp_lines <- comp_lines[str_detect(comp_lines, "^\\s*C\\d{5}\\s+")]
  comp_dt <- data.table(raw = comp_lines)
  comp_dt[, compound_id := str_match(raw, "^\\s*(C\\d{5})")[, 2]]
  comp_dt[, compound_name := str_trim(str_replace(raw, "^\\s*C\\d{5}\\s+", ""))]
  comp_dt <- comp_dt[!is.na(compound_id) & compound_id != ""]
  compounds <- setNames(comp_dt$compound_name, comp_dt$compound_id)

  list(gene_map = gene_dt, compounds = compounds)
}

gene_map <- data.table()
gene_order <- character()
kegg_compounds <- character()
kegg_compound_ids <- character()

if (file.exists(txt_path)) {
  flat <- tryCatch(parse_kegg_flat(txt_path), error = function(e) NULL)
  if (!is.null(flat) && nrow(flat$gene_map) > 0 && length(flat$compounds) > 0) {
    gene_map <- flat$gene_map
    gene_order <- gene_map$gene
    kegg_compounds <- flat$compounds
    kegg_compound_ids <- names(kegg_compounds)
  }
}

pw <- NULL
if (nrow(gene_map) == 0 || length(kegg_compound_ids) == 0) {
  pw <- tryCatch(keggGet(pathway_id)[[1]], error = function(e) NULL)
}

if (is.null(pw) && file.exists(kgml_path) && (nrow(gene_map) == 0 || length(kegg_compound_ids) == 0)) {
  # Fallback: parse KGML for node names/labels (genes/compounds) without KEGG API call
  message("KEGG API timeout; falling back to local KGML: ", kgml_path)
  doc <- xml2::read_xml(kgml_path)
  entries <- xml2::xml_find_all(doc, ".//entry")
  get_attr <- function(node, name) xml2::xml_attr(node, name)

  # Gene symbols from graphics label
  gene_nodes_tmp <- rbindlist(lapply(entries, function(e) {
    if (get_attr(e, "type") != "gene") return(NULL)
    g <- xml2::xml_find_first(e, ".//graphics")
    data.table(label = get_attr(g, "name"))
  }), fill = TRUE)
  gene_nodes_tmp <- gene_nodes_tmp[!is.na(label)]
  gene_syms <- unique(str_trim(sub("\\.\\.\\.$", "", gene_nodes_tmp$label)))
  gene_syms <- gene_syms[gene_syms != ""]

  # Gene Entrez IDs from entry "name" (cfa:xxxx ...)
  gene_ids_tmp <- rbindlist(lapply(entries, function(e) {
    if (get_attr(e, "type") != "gene") return(NULL)
    data.table(name = get_attr(e, "name"))
  }), fill = TRUE)
  gene_ids_tmp <- gene_ids_tmp[!is.na(name)]
  entrez_ids <- unique(str_extract(gene_ids_tmp$name, "\\d+"))
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]

  # Compound IDs from entry "name" (cpd:Cxxxxx)
  comp_ids <- rbindlist(lapply(entries, function(e) {
    if (get_attr(e, "type") != "compound") return(NULL)
    data.table(cid = str_extract(get_attr(e, "name"), "C\\d{5}"))
  }), fill = TRUE)
  comp_ids <- unique(comp_ids[!is.na(cid)]$cid)

  gene_map <- data.table(entrez = entrez_ids, gene = gene_syms[seq_len(min(length(gene_syms), length(entrez_ids)))])
  gene_map <- unique(gene_map[!is.na(entrez) & entrez != "" & !is.na(gene) & gene != ""])
  gene_order <- gene_map$gene

  kegg_compound_ids <- comp_ids
  kegg_compounds <- setNames(kegg_compound_ids, kegg_compound_ids)
} else if (!is.null(pw) && (nrow(gene_map) == 0 || length(kegg_compound_ids) == 0)) {
  # GENE is alternating: Entrez ID, then "SYMBOL; description ..."
  gene_ids  <- pw$GENE[seq(1, length(pw$GENE), by = 2)]
  gene_desc <- pw$GENE[seq(2, length(pw$GENE), by = 2)]
  kegg_genes <- str_trim(sub(";.*$", "", gene_desc))
  kegg_genes <- kegg_genes[kegg_genes != ""]

  gene_map <- data.table(
    entrez = as.character(gene_ids),
    gene   = as.character(kegg_genes)
  )
  gene_map <- unique(gene_map[gene != "" & !is.na(gene)])
  gene_order <- gene_map$gene

  # COMPOUND is a named vector: names = CXXXXX, values = compound name/abbrev
  kegg_compounds <- pw$COMPOUND
  kegg_compound_ids <- names(kegg_compounds)
}

if (nrow(gene_map) == 0 || length(kegg_compound_ids) == 0) {
  stop("Failed to obtain KEGG pathway members (gene/compound). Please ensure at least one of these exists: ",
       txt_path, " or ", kgml_path, " (or KEGG REST is reachable).")
}

# Save KEGG members for traceability
members_dt <- rbind(
  data.table(type = "gene", id = gene_map$entrez, name = gene_map$gene),
  data.table(type = "compound", id = kegg_compound_ids, name = unname(kegg_compounds))
)
fwrite(members_dt, file.path(out_dir, "pyrimidine_metabolism_kegg_members.tsv"), sep = "\t")

# Lookup from compound label -> ID (e.g. "UMP" -> "C00105")
compound_label_to_id <- setNames(names(kegg_compounds), unname(kegg_compounds))

# ======================
# 2) ID mapping (MyGene.info)
# ======================
mygene_entrez_to_ensembl <- function(entrez_ids) {
  entrez_ids <- unique(as.character(entrez_ids))
  entrez_ids <- entrez_ids[!is.na(entrez_ids) & entrez_ids != ""]
  if (length(entrez_ids) == 0) return(data.table(entrez = character(), gene_id = character()))

  url <- "https://mygene.info/v3/gene"
  chunks <- split(entrez_ids, ceiling(seq_along(entrez_ids) / 200))

  out <- rbindlist(lapply(chunks, function(ids) {
    r <- tryCatch(
      httr::POST(url, body = list(ids = ids, fields = "symbol,ensembl.gene"), encode = "json", httr::timeout(20)),
      error = function(e) NULL
    )
    if (is.null(r)) return(data.table(entrez = character(), gene_id = character()))

    ok <- tryCatch({httr::stop_for_status(r); TRUE}, error = function(e) FALSE)
    if (!ok) return(data.table(entrez = character(), gene_id = character()))

    js <- httr::content(r, as = "text", encoding = "UTF-8")
    dt <- as.data.table(jsonlite::fromJSON(js))
    setnames(dt, c("query", "ensembl"), c("entrez", "gene_id"), skip_absent = TRUE)
    dt <- dt[!is.na(entrez) & !is.na(gene_id) & gene_id != ""]
    dt[, .(entrez = as.character(entrez), gene_id = as.character(gene_id))]
  }), use.names = TRUE, fill = TRUE)

  unique(out)
}

ens_map <- mygene_entrez_to_ensembl(gene_map$entrez)
ens_map <- merge(ens_map, gene_map, by = "entrez", all.x = TRUE)
ens_map <- ens_map[!is.na(gene) & gene != ""]
ens_ids <- unique(ens_map$gene_id)
fwrite(ens_map, file.path(out_dir, "pyrimidine_ensembl_mapping.tsv"), sep = "\t")

mygene_entrez_to_uniprot <- function(entrez_ids) {
  entrez_ids <- unique(as.character(entrez_ids))
  entrez_ids <- entrez_ids[!is.na(entrez_ids) & entrez_ids != ""]
  if (length(entrez_ids) == 0) return(data.table(entrez = character(), uniprot = character()))

  url <- "https://mygene.info/v3/gene"
  chunks <- split(entrez_ids, ceiling(seq_along(entrez_ids) / 200))

  out <- rbindlist(lapply(chunks, function(ids) {
    r <- tryCatch(
      httr::POST(url, body = list(ids = ids, fields = "uniprot"), encode = "json", httr::timeout(20)),
      error = function(e) NULL
    )
    if (is.null(r)) return(data.table(entrez = character(), uniprot = character()))

    ok <- tryCatch({httr::stop_for_status(r); TRUE}, error = function(e) FALSE)
    if (!ok) return(data.table(entrez = character(), uniprot = character()))

    js <- httr::content(r, as = "text", encoding = "UTF-8")
    lst <- jsonlite::fromJSON(js, simplifyVector = FALSE)

    rbindlist(lapply(lst, function(x) {
      if (is.null(x$query) || isTRUE(x$notfound) || is.null(x$uniprot)) return(NULL)

      unis <- c(x$uniprot$`Swiss-Prot`, x$uniprot$TrEMBL)
      unis <- unis[!is.na(unis) & unis != ""]
      if (length(unis) == 0) return(NULL)

      data.table(entrez = as.character(x$query), uniprot = as.character(unis))
    }), fill = TRUE)
  }), use.names = TRUE, fill = TRUE)

  unique(out[!is.na(entrez) & entrez != "" & !is.na(uniprot) & uniprot != ""])
}

# ======================
# 3) Load transcript (DEG) tables
# ======================
read_deg <- function(file, genotype_label) {
  dt <- as.data.table(readxl::read_excel(file.path(base_dir, file)))
  id_col <- intersect(c("Unnamed: 0", "...1"), names(dt))
  if (length(id_col) == 0) id_col <- names(dt)[1]
  setnames(dt, id_col[1], "gene_id")
  stopifnot(all(c("log2FoldChange", "pvalue", "padj") %in% names(dt)))

  dt <- dt[gene_id %in% ens_ids]
  if (nrow(dt) == 0) return(data.table(gene = character(), genotype = character(), log2FC = numeric(), p_adj = numeric()))

  dt <- merge(dt, ens_map[, .(gene_id, gene)], by = "gene_id", all.x = TRUE)
  dt <- dt[!is.na(gene) & gene != ""]

  dt[, padj := suppressWarnings(as.numeric(padj))]
  dt[, log2FoldChange := suppressWarnings(as.numeric(log2FoldChange))]
  dt[, p_adj2 := fifelse(is.na(padj), 1, padj)]
  dt <- dt[order(p_adj2, -abs(log2FoldChange))]
  dt <- dt[, .SD[1], by = gene]

  dt[, `:=`(
    genotype = genotype_label,
    log2FC   = as.numeric(log2FoldChange),
    p_adj    = as.numeric(padj)
  )]
  dt[, .(gene, genotype, log2FC, p_adj)]
}

rna_all <- rbindlist(list(
  read_deg("DEG_S.xlsx", "sh"),
  read_deg("DEG_C.xlsx", "C"),
  read_deg("DEG_O.xlsx", "OE")
), use.names = TRUE, fill = TRUE)

# ======================
# 4) Load protein (DEP) table
# ======================
dep <- as.data.table(readxl::read_excel(file.path(base_dir, "DEP.xlsx")))
stopifnot(all(c(
  "ID",
  "C_H_vs_C_NC_p.adj", "OE_H_vs_OE_NC_p.adj", "SH_H_vs_SH_NC_p.adj",
  "C_H_vs_C_NC_ratio", "OE_H_vs_OE_NC_ratio", "SH_H_vs_SH_NC_ratio"
) %in% names(dep)))

uniprot_map <- mygene_entrez_to_uniprot(gene_map$entrez)
uniprot_map <- uniprot_map[uniprot %in% dep$ID]
uniprot_map <- merge(uniprot_map, gene_map, by = "entrez", all.x = TRUE)
uniprot_map <- uniprot_map[!is.na(gene) & gene != ""]
fwrite(uniprot_map, file.path(out_dir, "pyrimidine_uniprot_mapping.tsv"), sep = "\t")

dep <- merge(dep, uniprot_map[, .(uniprot, gene)], by.x = "ID", by.y = "uniprot", all.x = FALSE)
dep <- dep[gene %in% gene_order]

# Select one representative protein row per gene (lowest min p.adj, then largest mean |log2FC|)
if (nrow(dep) > 0) {
  dep[, min_padj := pmin(`C_H_vs_C_NC_p.adj`, `OE_H_vs_OE_NC_p.adj`, `SH_H_vs_SH_NC_p.adj`, na.rm = TRUE)]
  dep[is.infinite(min_padj), min_padj := NA_real_]
  dep[, mean_abs_fc := rowMeans(abs(as.matrix(.SD)), na.rm = TRUE),
      .SDcols = c("C_H_vs_C_NC_ratio", "OE_H_vs_OE_NC_ratio", "SH_H_vs_SH_NC_ratio")]

  dep[, min_padj2 := fifelse(is.na(min_padj), 1, min_padj)]
  dep <- dep[order(min_padj2, -mean_abs_fc)]
  dep <- dep[, .SD[1], by = gene]
}

prot_fc <- dep[, .(
  gene,
  sh = as.numeric(`SH_H_vs_SH_NC_ratio`),
  C  = as.numeric(`C_H_vs_C_NC_ratio`),
  OE = as.numeric(`OE_H_vs_OE_NC_ratio`)
)]
prot_padj <- dep[, .(
  gene,
  sh = as.numeric(`SH_H_vs_SH_NC_p.adj`),
  C  = as.numeric(`C_H_vs_C_NC_p.adj`),
  OE = as.numeric(`OE_H_vs_OE_NC_p.adj`)
)]

# ======================
# 5) Load metabolite (DEM) tables
# ======================
parse_compound_id <- function(kegg_id, name) {
  # Prefer explicit KEGG compound ID in the table
  cid <- str_extract(as.character(kegg_id), "C\\d{5}")
  if (!is.na(cid)) return(cid)

  # Try abbreviation in parentheses, e.g. "Uridine monophosphate (UMP)"
  abbr <- str_match(as.character(name), "\\\\(([A-Za-z0-9\\\\-]+)\\\\)")[, 2]
  if (!is.na(abbr) && abbr %in% names(compound_label_to_id)) return(compound_label_to_id[[abbr]])

  # Try exact match against KEGG compound labels (often short like UMP/UDP/CTP)
  nm <- as.character(name)
  if (!is.na(nm) && nm %in% names(compound_label_to_id)) return(compound_label_to_id[[nm]])

  NA_character_
}

read_dem <- function(file, genotype_label) {
  dt <- as.data.table(readxl::read_excel(file.path(base_dir, file)))
  stopifnot(all(c("log2FC", "pval", "VIP", "Name", "Kegg_ID") %in% names(dt)))

  dt[, compound_id := mapply(parse_compound_id, Kegg_ID, Name)]
  dt <- dt[!is.na(compound_id) & compound_id %in% kegg_compound_ids]

  if (nrow(dt) == 0) {
    return(data.table(compound_id = character(), display = character(), genotype = character(),
                      log2FC = numeric(), pval = numeric(), VIP = numeric()))
  }

  # Prefer KEGG label for display when available
  dt[, display := fifelse(compound_id %in% names(kegg_compounds), unname(kegg_compounds[compound_id]), as.character(Name))]

  # Deduplicate per compound (lowest p, then highest |FC|, then VIP)
  dt <- dt[order(pval, -abs(log2FC), -VIP)]
  dt <- dt[, .SD[1], by = compound_id]

  dt[, `:=`(
    genotype = genotype_label,
    log2FC   = as.numeric(log2FC),
    pval     = as.numeric(pval),
    VIP      = as.numeric(VIP)
  )]

  dt[, .(compound_id, display, genotype, log2FC, pval, VIP)]
}

met_all <- rbindlist(list(
  read_dem("shH_shN_DEM.xlsx", "sh"),
  read_dem("CH_CN_DEM.xlsx", "C"),
  read_dem("OEH_OEN_DEM.xlsx", "OE")
), use.names = TRUE, fill = TRUE)

# NOTE: keep all pathway-mapped metabolites found in DEM tables.
# Significance is shown by stars (p-value + optional VIP threshold).

# ======================
# 6) Build matrices (genes: RNA+Protein; metabolites)
# ======================
make_wide <- function(dt, value_col, row_id_col) {
  if (nrow(dt) == 0) return(data.table())
  wide <- dcast(dt, as.formula(paste(row_id_col, "~ genotype")), value.var = value_col)
  for (g in GENOTYPES) if (!(g %in% names(wide))) wide[, (g) := NA_real_]
  setcolorder(wide, c(row_id_col, GENOTYPES))
  wide
}

# RNA matrices (gene symbol row key)
rna_fc_wide   <- make_wide(rna_all, "log2FC", "gene")
rna_padj_wide <- make_wide(rna_all, "p_adj", "gene")

# Protein matrices (already wide)
for (g in GENOTYPES) if (!(g %in% names(prot_fc))) prot_fc[, (g) := NA_real_]
for (g in GENOTYPES) if (!(g %in% names(prot_padj))) prot_padj[, (g) := NA_real_]
prot_fc <- prot_fc[, c("gene", GENOTYPES), with = FALSE]
prot_padj <- prot_padj[, c("gene", GENOTYPES), with = FALSE]

# Gene order: KEGG order, keep those that exist in RNA or protein
genes_present <- union(rna_fc_wide$gene, prot_fc$gene)
gene_order <- gene_order[gene_order %in% genes_present]

if (length(gene_order) == 0) {
  stop("No pyrimidine-metabolism genes found in RNA/protein tables after mapping.")
}

to_mat <- function(wide_dt, row_key, row_order) {
  if (is.null(row_order)) row_order <- character()
  if (nrow(wide_dt) == 0 || length(row_order) == 0) {
    m <- matrix(NA_real_, nrow = length(row_order), ncol = length(GENOTYPES))
    rownames(m) <- row_order
    colnames(m) <- GENOTYPES
    return(m)
  }
  key_dt <- data.table()
  key_dt[, key := row_order]
  tmp <- merge(key_dt, wide_dt, by.x = "key", by.y = row_key, all.x = TRUE, sort = FALSE)
  m <- as.matrix(tmp[, ..GENOTYPES])
  rownames(m) <- row_order
  colnames(m) <- GENOTYPES
  m
}

rna_fc_mat   <- to_mat(rna_fc_wide, "gene", gene_order)
rna_padj_mat <- to_mat(rna_padj_wide, "gene", gene_order)
prot_fc_mat  <- to_mat(prot_fc, "gene", gene_order)
prot_p_mat   <- to_mat(prot_padj, "gene", gene_order)

# Interleave RNA/Protein rows: 2 x 3 per gene
mat_gene <- do.call(rbind, lapply(gene_order, function(g) rbind(rna_fc_mat[g, ], prot_fc_mat[g, ])))
p_gene   <- do.call(rbind, lapply(gene_order, function(g) rbind(rna_padj_mat[g, ], prot_p_mat[g, ])))

omics_vec <- rep(c("Transcript", "Protein"), times = length(gene_order))
row_labels_gene <- rep("", nrow(mat_gene))
row_labels_gene[seq(1, nrow(mat_gene), by = 2)] <- gene_order

rownames(mat_gene) <- paste0(rep(gene_order, each = 2), "|", omics_vec)
rownames(p_gene) <- rownames(mat_gene)
colnames(mat_gene) <- GENOTYPES
colnames(p_gene) <- GENOTYPES

# Metabolite matrices
met_fc_wide  <- make_wide(met_all, "log2FC", "compound_id")
met_p_wide   <- make_wide(met_all, "pval", "compound_id")
met_vip_wide <- make_wide(met_all, "VIP", "compound_id")

comp_present <- met_fc_wide$compound_id
if (is.null(comp_present)) comp_present <- character()
comp_order <- kegg_compound_ids[kegg_compound_ids %in% comp_present]

mat_met <- to_mat(met_fc_wide, "compound_id", comp_order)
p_met   <- to_mat(met_p_wide, "compound_id", comp_order)
vip_met <- to_mat(met_vip_wide, "compound_id", comp_order)

display_map <- unique(met_all[, .(compound_id, display)])
display_map <- display_map[order(match(compound_id, comp_order))]
met_labels <- display_map$display[match(comp_order, display_map$compound_id)]
met_labels[is.na(met_labels) | met_labels == ""] <- comp_order[is.na(met_labels) | met_labels == ""]
rownames(mat_met) <- met_labels
rownames(p_met) <- met_labels
rownames(vip_met) <- met_labels
colnames(mat_met) <- GENOTYPES
colnames(p_met) <- GENOTYPES
colnames(vip_met) <- GENOTYPES

# ======================
# 7) Stars, scaling, and plotting
# ======================
stars_from_p <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}

stars_gene <- stars_from_p(p_gene)

# Metabolite significance: p-value + optional VIP
stars_met <- stars_from_p(p_met)
if (!is.na(VIP_CUTOFF)) {
  stars_met[vip_met < VIP_CUTOFF] <- ""
}

round_up <- function(x, step = 0.5) ceiling(x / step) * step
all_vals <- c(as.numeric(mat_gene), as.numeric(mat_met))
all_vals <- all_vals[is.finite(all_vals)]
if (length(all_vals) == 0) stop("No finite log2FC values found for plotting.")

cap_abs <- if (is.null(LOG2FC_CAP_ABS)) {
  max(2.0, round_up(as.numeric(stats::quantile(abs(all_vals), probs = 0.95, na.rm = TRUE)), 0.5))
} else {
  as.numeric(LOG2FC_CAP_ABS)
}

clip_mat <- function(m) pmax(pmin(m, cap_abs), -cap_abs)
mat_gene_clip <- clip_mat(mat_gene)
mat_met_clip  <- clip_mat(mat_met)

col_fun <- circlize::colorRamp2(c(-cap_abs, 0, cap_abs), c(COL_DOWN, COL_MID, COL_UP))

ht_row_anno <- rowAnnotation(
  Omics = omics_vec,
  col = list(Omics = c("Transcript" = "grey25", "Protein" = "grey65")),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  width = unit(4, "mm")
)

split_gene <- rep(seq_along(gene_order), each = 2)

ht_gene <- Heatmap(
  mat_gene_clip,
  name = "log2FC",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = split_gene,
  row_title = NULL,
  row_gap = unit(0.9, "mm"),
  cluster_row_slices = FALSE,
  left_annotation = ht_row_anno,
  show_row_names = TRUE,
  row_labels = row_labels_gene,
  row_names_gp = gpar(fontsize = 9),
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_title = "Pyrimidine metabolism (H2O2 vs Normal): log2FC",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(
    title = "log2FC",
    title_position = "topcenter",
    legend_height = unit(35, "mm")
  ),
  cell_fun = function(j, i, x, y, w, h, fill) {
    lab <- stars_gene[i, j]
    if (lab != "") grid.text(lab, x, y, gp = gpar(fontsize = 10, fontface = "bold"))
  }
)

ht_met <- Heatmap(
  mat_met_clip,
  name = "log2FC_met",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  show_column_names = FALSE,
  row_title = "Metabolites",
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  row_title_rot = 0,
  show_heatmap_legend = FALSE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    lab <- stars_met[i, j]
    if (lab != "") grid.text(lab, x, y, gp = gpar(fontsize = 10, fontface = "bold"))
  }
)

ht_list <- ht_gene %v% ht_met

# Figure size (dynamic)
n_rows_total <- nrow(mat_gene) + nrow(mat_met)
fig_h <- max(6.5, min(14, 0.12 * n_rows_total + 2.5))
fig_w <- 7.2

out_pdf <- file.path(out_dir, "pyrimidine_metabolism_multiomics_heatmap.pdf")
out_png <- file.path(out_dir, "pyrimidine_metabolism_multiomics_heatmap.png")

pdf(out_pdf, width = fig_w, height = fig_h, useDingbats = FALSE)
draw(
  ht_list,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(2, 2, 2, 2), "mm")
)
grid.text(
  sprintf("Stars: RNA/Protein padj <= %.02f; Metabolite p <= %.02f%s. log2FC capped at ±%.1f.",
          P_ADJ_CUTOFF, P_MET_CUTOFF,
          if (!is.na(VIP_CUTOFF)) sprintf(" & VIP >= %.1f", VIP_CUTOFF) else "",
          cap_abs),
  x = unit(2, "mm"),
  y = unit(2, "mm"),
  just = c("left", "bottom"),
  gp = gpar(fontsize = 9, col = "grey30")
)
dev.off()

png(out_png, width = fig_w, height = fig_h, units = "in", res = 300)
draw(
  ht_list,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(2, 2, 2, 2), "mm")
)
grid.text(
  sprintf("Stars: RNA/Protein padj <= %.02f; Metabolite p <= %.02f%s. log2FC capped at ±%.1f.",
          P_ADJ_CUTOFF, P_MET_CUTOFF,
          if (!is.na(VIP_CUTOFF)) sprintf(" & VIP >= %.1f", VIP_CUTOFF) else "",
          cap_abs),
  x = unit(2, "mm"),
  y = unit(2, "mm"),
  just = c("left", "bottom"),
  gp = gpar(fontsize = 9, col = "grey30")
)
dev.off()

message("Saved: ", out_pdf)
message("Saved: ", out_png)

# ============================================================
# 8) Overlay mini-heatmaps on KEGG pathway diagram (KGML + PNG)
# ============================================================
kgml_path <- file.path(out_dir, "cfa00240.kgml")
png_path  <- file.path(out_dir, "cfa00240.png")

if (!file.exists(kgml_path) || !file.exists(png_path)) {
  message("KGML/PNG not found; download from KEGG REST and re-run to enable pathway overlay.")
  quit(status = 0)
}

read_kgml_nodes <- function(kgml) {
  doc <- xml2::read_xml(kgml)
  entries <- xml2::xml_find_all(doc, ".//entry")

  get_attr <- function(node, name) xml2::xml_attr(node, name)

  rbindlist(lapply(entries, function(e) {
    g <- xml2::xml_find_first(e, ".//graphics")
    if (is.na(get_attr(g, "x"))) return(NULL)
    data.table(
      id = as.integer(get_attr(e, "id")),
      type = get_attr(e, "type"),
      name = get_attr(e, "name"),
      label = get_attr(g, "name"),
      x = as.numeric(get_attr(g, "x")),
      y = as.numeric(get_attr(g, "y")),
      w = as.numeric(get_attr(g, "width")),
      h = as.numeric(get_attr(g, "height")),
      gtype = get_attr(g, "type"),
      fgcolor = get_attr(g, "fgcolor"),
      bgcolor = get_attr(g, "bgcolor")
    )
  }), fill = TRUE)
}

nodes <- read_kgml_nodes(kgml_path)

gene_nodes <- nodes[type == "gene"]
compound_nodes <- nodes[type == "compound"]

clean_label <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\n", "")
  str_trim(x)
}

# KGML gene nodes may represent multiple genes (label with "...")
# Use entry 'name' (cfa:XXXX) to map to gene symbols, then aggregate per node.
gene_nodes[, label_clean := clean_label(label)]
gene_nodes[, entrez_ids := str_extract_all(name, "\\d+")]
gene_nodes[, gene_syms := lapply(entrez_ids, function(ids) unique(gene_map[entrez %in% ids]$gene))]
gene_nodes[, gene_syms := lapply(seq_len(.N), function(i) {
  syms <- gene_syms[[i]]
  if (length(syms) > 0) return(syms)
  lab <- sub("\\.\\.\\.$", "", label_clean[i])
  if (!is.na(lab) && lab %in% gene_order) return(lab)
  character()
})]
gene_nodes <- gene_nodes[sapply(gene_syms, length) > 0]

# Compound node labels: use compound_id
compound_nodes[, compound_id := str_extract(name, "C\\d{5}")]
compound_nodes <- compound_nodes[!is.na(compound_id)]

# Map compound_id -> heatmap row label used in mat_met
met_id_to_label <- setNames(rownames(mat_met), comp_order)
compound_nodes[, met_label := met_id_to_label[compound_id]]
compound_nodes <- compound_nodes[!is.na(met_label)]

img <- png::readPNG(png_path)
img_ref <- img
ref_w <- dim(img_ref)[2]
ref_h <- dim(img_ref)[1]

# Use 2x base image for sharper publication output (keeps original KEGG layout)
base_png <- png_path
if (isTRUE(USE_KEGG_BASE_2X)) {
  if (!file.exists(KEGG_BASE_PNG_2X)) {
    message("Downloading KEGG 2x base image: ", KEGG_BASE_URL_2X)
    tryCatch(
      utils::download.file(KEGG_BASE_URL_2X, KEGG_BASE_PNG_2X, mode = "wb", quiet = TRUE),
      error = function(e) message("Warning: failed to download 2x base image; fallback to 1x: ", conditionMessage(e))
    )
  }
  if (file.exists(KEGG_BASE_PNG_2X)) {
    base_png <- KEGG_BASE_PNG_2X
  }
}

img <- png::readPNG(base_png)
img_w <- dim(img)[2]
img_h <- dim(img)[1]

scale_x <- img_w / ref_w
scale_y <- img_h / ref_h

# Helper: draw a matrix heatmap inside a rectangle in image coordinate system
draw_mini_heatmap <- function(mat, rect_x, rect_y, rect_w, rect_h, star_mat = NULL,
                              nrow_cells, ncol_cells, col_fun) {
  # mat is numeric matrix (rows x cols) but we draw in the order given
  # rect_* in KGML pixels with origin top-left
  if (nrow_cells <= 0 || ncol_cells <= 0) return(invisible(NULL))

  # Convert rect center-based to top-left
  x0 <- rect_x - rect_w / 2
  y0 <- rect_y - rect_h / 2

  cell_w <- rect_w / ncol_cells
  cell_h <- rect_h / nrow_cells

  for (i in seq_len(nrow_cells)) {
    for (j in seq_len(ncol_cells)) {
      v <- mat[i, j]
      if (!is.finite(v)) next
      fill <- col_fun(v)

      # Cell center in kgml pixels
      cx <- x0 + (j - 0.5) * cell_w
      cy <- y0 + (i - 0.5) * cell_h

      # Convert to grid NPC: x from left, y from bottom
      grid.rect(
        x = unit(cx / img_w, "npc"),
        y = unit(1 - cy / img_h, "npc"),
        width = unit(cell_w / img_w, "npc"),
        height = unit(cell_h / img_h, "npc"),
        gp = gpar(col = "grey30", lwd = 0.2, fill = fill)
      )

      if (!is.null(star_mat)) {
        s <- star_mat[i, j]
        if (!is.na(s) && s != "") {
          grid.text(
            s,
            x = unit(cx / img_w, "npc"),
            y = unit(1 - cy / img_h, "npc"),
            gp = gpar(fontsize = 6.5, fontface = "bold", col = "black")
          )
        }
      }
    }
  }
}

# Prepare per-gene 2x3 matrix (Transcript row then Protein row)
agg_fc <- function(m) {
  if (nrow(m) == 0) return(rep(NA_real_, length(GENOTYPES)))
  v <- colMeans(m, na.rm = TRUE)
  v[is.nan(v)] <- NA_real_
  as.numeric(v)
}
agg_pmin <- function(m) {
  if (nrow(m) == 0) return(rep(NA_real_, length(GENOTYPES)))
  v <- apply(m, 2, function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    min(x)
  })
  as.numeric(v)
}

gene_fc_lookup <- function(gene_syms) {
  gene_syms <- unique(as.character(gene_syms))
  gene_syms <- gene_syms[!is.na(gene_syms) & gene_syms != ""]

  rna <- rna_fc_mat[intersect(gene_syms, rownames(rna_fc_mat)), , drop = FALSE]
  pro <- prot_fc_mat[intersect(gene_syms, rownames(prot_fc_mat)), , drop = FALSE]
  rbind(agg_fc(rna), agg_fc(pro))
}
gene_star_lookup <- function(gene_syms) {
  gene_syms <- unique(as.character(gene_syms))
  gene_syms <- gene_syms[!is.na(gene_syms) & gene_syms != ""]

  rna_p <- rna_padj_mat[intersect(gene_syms, rownames(rna_padj_mat)), , drop = FALSE]
  pro_p <- prot_p_mat[intersect(gene_syms, rownames(prot_p_mat)), , drop = FALSE]
  stars_from_p(rbind(agg_pmin(rna_p), agg_pmin(pro_p)))
}

# Prepare metabolite 1x3 matrix
met_fc_lookup <- function(met_label) {
  mat_met[met_label, , drop = FALSE]
}
met_star_lookup <- function(met_label) {
  stars_from_p(p_met[met_label, , drop = FALSE])
}

# Clip values for consistent color
cfun <- circlize::colorRamp2(c(-cap_abs, 0, cap_abs), c(COL_DOWN, COL_MID, COL_UP))

gene_nodes_plot <- gene_nodes[sapply(gene_syms, function(syms) {
  any(is.finite(as.numeric(gene_fc_lookup(syms))))
})]

overlay_tag <- if (identical(normalizePath(base_png, winslash = "/"), normalizePath(png_path, winslash = "/"))) "1x" else "kegg2x"
overlay_pdf <- file.path(out_dir, sprintf("pyrimidine_metabolism_pathway_overlay_%s.pdf", overlay_tag))
overlay_png <- file.path(out_dir, sprintf("pyrimidine_metabolism_pathway_overlay_%s.png", overlay_tag))

draw_overlay <- function(device_fun) {
  device_fun()
  grid.newpage()
  grid.raster(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)

  # Draw gene mini-heatmaps
  for (k in seq_len(nrow(gene_nodes_plot))) {
    g <- gene_nodes_plot[k]
    mat <- gene_fc_lookup(g$gene_syms[[1]])
    mat <- clip_mat(mat)
    stars <- gene_star_lookup(g$gene_syms[[1]])
    draw_mini_heatmap(
      mat,
      g$x * scale_x, g$y * scale_y,
      g$w * scale_x, g$h * scale_y,
      star_mat = stars, nrow_cells = 2, ncol_cells = 3, col_fun = cfun
    )
  }

  # Draw metabolite mini-heatmaps (1x3)
  for (k in seq_len(nrow(compound_nodes))) {
    m <- compound_nodes[k]
    mat <- met_fc_lookup(m$met_label)
    mat <- clip_mat(mat)
    stars <- met_star_lookup(m$met_label)
    rect_w <- max(m$w, 46) * scale_x
    rect_h <- max(m$h, 14) * scale_y
    draw_mini_heatmap(
      mat,
      m$x * scale_x, m$y * scale_y,
      rect_w, rect_h,
      star_mat = stars, nrow_cells = 1, ncol_cells = 3, col_fun = cfun
    )
  }

  # Legend and annotation panel (compact, bottom-right)
  grid.rect(
    x = unit(0.98, "npc"), y = unit(0.02, "npc"),
    width = unit(0.28, "npc"), height = unit(0.12, "npc"),
    just = c("right", "bottom"),
    gp = gpar(fill = adjustcolor("white", alpha.f = 0.85), col = "grey80", lwd = 0.3)
  )

  # Mini heatmap key: rows = omics, cols = genotype
  key_x0 <- 0.74; key_y0 <- 0.03
  key_w <- 0.10; key_h <- 0.07
  key_cell_w <- key_w / 3
  key_cell_h <- key_h / 2

  # Draw empty grid and labels
  for (i in 1:2) {
    for (j in 1:3) {
      grid.rect(
        x = unit(key_x0 + (j - 0.5) * key_cell_w, "npc"),
        y = unit(key_y0 + (2 - i + 0.5) * key_cell_h, "npc"),
        width = unit(key_cell_w, "npc"),
        height = unit(key_cell_h, "npc"),
        gp = gpar(fill = "white", col = "grey50", lwd = 0.3)
      )
    }
  }
  grid.text("RNA", x = unit(key_x0 - 0.01, "npc"), y = unit(key_y0 + key_cell_h * 1.5, "npc"),
            just = "right", gp = gpar(fontsize = 7))
  grid.text("Prot", x = unit(key_x0 - 0.01, "npc"), y = unit(key_y0 + key_cell_h * 0.5, "npc"),
            just = "right", gp = gpar(fontsize = 7))
  for (j in 1:3) {
    grid.text(GENOTYPES[j], x = unit(key_x0 + (j - 0.5) * key_cell_w, "npc"),
              y = unit(key_y0 + key_h + 0.01, "npc"), gp = gpar(fontsize = 7, fontface = "bold"))
  }

  # Colorbar
  lgd <- Legend(
    title = "log2FC",
    col_fun = cfun,
    at = c(-cap_abs, 0, cap_abs),
    labels = c(paste0("-", cap_abs), "0", paste0("+", cap_abs)),
    direction = "horizontal",
    legend_width = unit(38, "mm"),
    title_position = "topcenter"
  )
  pushViewport(viewport(x = unit(0.92, "npc"), y = unit(0.09, "npc"), just = c("right", "bottom")))
  draw(lgd)
  popViewport()

  grid.text("* p<0.05  ** p<0.01  *** p<0.001", x = unit(0.74, "npc"), y = unit(0.02, "npc"),
            just = c("left", "bottom"), gp = gpar(fontsize = 7, col = "grey20"))
  grid.text("H2O2 vs Normal", x = unit(0.74, "npc"), y = unit(0.105, "npc"),
            just = c("left", "bottom"), gp = gpar(fontsize = 8, fontface = "bold"))

  dev.off()
}

overlay_w_in <- 7.2
overlay_h_in <- overlay_w_in * (img_h / img_w)
draw_overlay(function() grDevices::cairo_pdf(overlay_pdf, width = overlay_w_in, height = overlay_h_in))
draw_overlay(function() png(overlay_png, width = img_w, height = img_h, units = "px"))

message("Saved: ", overlay_pdf)
message("Saved: ", overlay_png)

# ============================================================
# 9) Vector redraw of pathway (publication-grade; no raster base)
# ============================================================
read_kgml_reactions <- function(kgml) {
  doc <- xml2::read_xml(kgml)
  rxns <- xml2::xml_find_all(doc, ".//reaction")

  rbindlist(lapply(rxns, function(r) {
    rid <- suppressWarnings(as.integer(xml2::xml_attr(r, "id")))
    rtype <- xml2::xml_attr(r, "type")

    subs <- xml2::xml_find_all(r, "./substrate")
    prods <- xml2::xml_find_all(r, "./product")

    sub_ids <- suppressWarnings(as.integer(xml2::xml_attr(subs, "id")))
    prod_ids <- suppressWarnings(as.integer(xml2::xml_attr(prods, "id")))

    data.table(
      enzyme_id = rid,
      reversible = identical(rtype, "reversible"),
      substrate_ids = list(sub_ids[!is.na(sub_ids)]),
      product_ids = list(prod_ids[!is.na(prod_ids)])
    )
  }), fill = TRUE)
}

draw_mini_heatmap_native <- function(mat, rect_x, rect_y, rect_w, rect_h, star_mat = NULL,
                                     nrow_cells, ncol_cells, col_fun,
                                     cell_border_col = "grey25", cell_border_lwd = 0.3,
                                     star_fontsize = 6.2) {
  if (nrow_cells <= 0 || ncol_cells <= 0) return(invisible(NULL))

  x0 <- rect_x - rect_w / 2
  y0 <- rect_y - rect_h / 2
  cell_w <- rect_w / ncol_cells
  cell_h <- rect_h / nrow_cells

  for (i in seq_len(nrow_cells)) {
    for (j in seq_len(ncol_cells)) {
      v <- mat[i, j]
      if (!is.finite(v)) next
      fill <- col_fun(v)

      cx <- x0 + (j - 0.5) * cell_w
      cy <- y0 + (i - 0.5) * cell_h

      grid.rect(
        x = unit(cx, "native"),
        y = unit(cy, "native"),
        width = unit(cell_w, "native"),
        height = unit(cell_h, "native"),
        gp = gpar(col = cell_border_col, lwd = cell_border_lwd, fill = fill)
      )

      if (!is.null(star_mat)) {
        s <- star_mat[i, j]
        if (!is.na(s) && s != "") {
          grid.text(
            s,
            x = unit(cx, "native"),
            y = unit(cy, "native"),
            gp = gpar(fontsize = star_fontsize, fontface = "bold", col = "black")
          )
        }
      }
    }
  }
}

grid_ellipse_native <- function(x, y, w, h, gp) {
  # Approximate ellipse with polygon (native units)
  theta <- seq(0, 2 * pi, length.out = 80)
  rx <- w / 2
  ry <- h / 2
  xs <- x + rx * cos(theta)
  ys <- y + ry * sin(theta)
  grid.polygon(x = unit(xs, "native"), y = unit(ys, "native"), gp = gp)
}

draw_kegg_shape <- function(node, fill = "white", col = "grey70", lwd = 0.35, alpha = 1.0) {
  gtype <- tolower(as.character(node$gtype))
  fill_col <- adjustcolor(fill, alpha.f = alpha)
  gp <- gpar(fill = fill_col, col = col, lwd = lwd)

  if (is.na(gtype) || gtype == "" || gtype == "rectangle") {
    grid.rect(
      x = unit(node$x, "native"),
      y = unit(node$y, "native"),
      width = unit(node$w, "native"),
      height = unit(node$h, "native"),
      gp = gp
    )
  } else if (gtype == "roundrectangle") {
    rr <- min(node$w, node$h) * 0.18
    grid.roundrect(
      x = unit(node$x, "native"),
      y = unit(node$y, "native"),
      width = unit(node$w, "native"),
      height = unit(node$h, "native"),
      r = unit(rr, "native"),
      gp = gp
    )
  } else if (gtype == "circle") {
    grid_ellipse_native(node$x, node$y, node$w, node$h, gp = gp)
  } else {
    # Fallback
    grid.rect(
      x = unit(node$x, "native"),
      y = unit(node$y, "native"),
      width = unit(node$w, "native"),
      height = unit(node$h, "native"),
      gp = gp
    )
  }
}

rect_bbox <- function(x, y, w, h) {
  list(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2)
}
rect_intersect <- function(a, b) {
  !(a$xmax <= b$xmin || a$xmin >= b$xmax || a$ymax <= b$ymin || a$ymin >= b$ymax)
}

bbox_outside <- function(bb, xmin, xmax, ymin, ymax) {
  # Amount (in points^2 proxy) outside plotting bounds; used as a penalty.
  dx <- max(0, xmin - bb$xmin) + max(0, bb$xmax - xmax)
  dy <- max(0, ymin - bb$ymin) + max(0, bb$ymax - ymax)
  dx + dy
}

place_labels_greedy <- function(dt, label_col, node_bbox_list, xmin, xmax, ymin, ymax,
                                gp, gap_seq = c(3, 5, 7, 9, 11, 14), order_pref = c("left", "right", "top", "bottom"),
                                pad_native = 1.2,
                                placed_bboxes_init = NULL,
                                attach_bboxes = TRUE,
                                allow_overlap_own_node = FALSE) {
  dt <- copy(dt)
  dt[, `:=`(label_x = NA_real_, label_y = NA_real_, label_pos = NA_character_)]

  placed_bboxes <- if (!is.null(placed_bboxes_init)) placed_bboxes_init else list()

  for (i in seq_len(nrow(dt))) {
    lab <- as.character(dt[[label_col]][i])
    if (is.na(lab) || lab == "") next

    # text size (native units)
    tg <- textGrob(lab, gp = gp)
    lw <- convertWidth(grobWidth(tg), "native", valueOnly = TRUE)
    lh <- convertHeight(grobHeight(tg), "native", valueOnly = TRUE)

    # own node bbox
    own_bb <- rect_bbox(dt$x[i], dt$y[i], dt$w[i], dt$h[i])

    best <- list(score = Inf, x = NA_real_, y = NA_real_, pos = NA_character_)

    left_x  <- dt$x[i] - dt$w[i] / 2
    right_x <- dt$x[i] + dt$w[i] / 2
    top_y   <- dt$y[i] - dt$h[i] / 2
    bot_y   <- dt$y[i] + dt$h[i] / 2

    for (gap in gap_seq) {
      for (pos in order_pref) {
        cx <- dt$x[i]; cy <- dt$y[i]
        if (pos == "left")  cx <- left_x  - gap - lw / 2
        if (pos == "right") cx <- right_x + gap + lw / 2
        if (pos == "top")   cy <- top_y   - gap - lh / 2
        if (pos == "bottom") cy <- bot_y  + gap + lh / 2

        bb <- rect_bbox(cx, cy, lw + 2 * pad_native, lh + 2 * pad_native)

        # overlaps with nodes (excluding own node)
        ov_nodes <- 0L
        for (nb in node_bbox_list) {
          if (!is.null(nb) && rect_intersect(bb, nb)) ov_nodes <- ov_nodes + 1L
        }
        if (isTRUE(allow_overlap_own_node) && rect_intersect(bb, own_bb)) {
          ov_nodes <- max(0L, ov_nodes - 1L)
        }

        # overlaps with already placed labels
        ov_labels <- 0L
        if (length(placed_bboxes) > 0) {
          for (pb in placed_bboxes) {
            if (!is.null(pb) && rect_intersect(bb, pb)) ov_labels <- ov_labels + 1L
          }
        }

        out_pen <- bbox_outside(bb, xmin, xmax, ymin, ymax)

        # heuristic scoring: strongly avoid overlaps; allow small outward movement if needed
        score <- 450 * out_pen + 600 * ov_labels + 150 * ov_nodes + 0.01 * gap

        if (score < best$score) {
          best <- list(score = score, x = cx, y = cy, pos = pos, lw = lw, lh = lh, bb = bb)
        }

        if (best$score == 0) break
      }
      if (best$score == 0) break
    }

    dt$label_x[i] <- best$x
    dt$label_y[i] <- best$y
    dt$label_pos[i] <- best$pos
    placed_bboxes <- c(placed_bboxes, list(best$bb))
  }

  if (isTRUE(attach_bboxes)) {
    attr(dt, "label_bboxes") <- placed_bboxes
  }
  dt
}

place_heatmap_boxes_greedy <- function(dt, anchor_w, anchor_h, rect_w, rect_h,
                                      node_bbox_list,
                                      xmin_ref, xmax_ref, ymin_ref, ymax_ref,
                                      gap_seq = c(0, 3, 5, 7, 9, 12, 15),
                                      order_pref = c("right", "left", "bottom", "top", "top-right", "bottom-right", "top-left", "bottom-left", "center"),
                                      pad_native = 1.2) {
  dt <- copy(dt)
  dt[, `:=`(
    hm_x = as.numeric(x),
    hm_y = as.numeric(y),
    hm_w = as.numeric(rect_w),
    hm_h = as.numeric(rect_h),
    hm_pos = "center"
  )]

  placed_bboxes <- list()

  for (i in seq_len(nrow(dt))) {
    ax <- as.numeric(dt$x[i]); ay <- as.numeric(dt$y[i])
    aw <- as.numeric(anchor_w[i]); ah <- as.numeric(anchor_h[i])
    rw <- as.numeric(rect_w[i]); rh <- as.numeric(rect_h[i])
    if (!is.finite(ax) || !is.finite(ay) || !is.finite(rw) || !is.finite(rh)) next

    r_anchor <- max(aw, ah, na.rm = TRUE) / 2
    if (!is.finite(r_anchor) || r_anchor <= 0) r_anchor <- 6

    own_bb <- rect_bbox(ax, ay, aw, ah)

    best <- list(score = Inf, x = ax, y = ay, pos = "center", bb = rect_bbox(ax, ay, rw, rh))

    for (gap in gap_seq) {
      for (pos in order_pref) {
        cx <- ax; cy <- ay
        dx <- r_anchor + gap + rw / 2
        dy <- r_anchor + gap + rh / 2

        if (pos == "right")        cx <- ax + dx
        if (pos == "left")         cx <- ax - dx
        if (pos == "bottom")       cy <- ay + dy
        if (pos == "top")          cy <- ay - dy
        if (pos == "top-right")    { cx <- ax + dx; cy <- ay - dy }
        if (pos == "bottom-right") { cx <- ax + dx; cy <- ay + dy }
        if (pos == "top-left")     { cx <- ax - dx; cy <- ay - dy }
        if (pos == "bottom-left")  { cx <- ax - dx; cy <- ay + dy }
        if (pos == "center")       { cx <- ax; cy <- ay }

        bb <- rect_bbox(cx, cy, rw + 2 * pad_native, rh + 2 * pad_native)

        ov_nodes <- 0L
        for (nb in node_bbox_list) {
          if (!is.null(nb) && rect_intersect(bb, nb)) ov_nodes <- ov_nodes + 1L
        }
        if (rect_intersect(bb, own_bb)) ov_nodes <- max(0L, ov_nodes - 1L)

        ov_hm <- 0L
        if (length(placed_bboxes) > 0) {
          for (pb in placed_bboxes) {
            if (!is.null(pb) && rect_intersect(bb, pb)) ov_hm <- ov_hm + 1L
          }
        }

        out_pen <- bbox_outside(bb, xmin_ref, xmax_ref, ymin_ref, ymax_ref)
        dist2 <- (cx - ax)^2 + (cy - ay)^2

        # Strongly discourage overlaps; allow distance only as a weak tiebreaker.
        score <- 1000 * out_pen + 500 * ov_hm + 120 * ov_nodes + 0.0005 * dist2 + 0.01 * gap
        if (score < best$score) {
          best <- list(score = score, x = cx, y = cy, pos = pos, bb = bb)
        }
        if (best$score == 0) break
      }
      if (best$score == 0) break
    }

    dt$hm_x[i] <- best$x
    dt$hm_y[i] <- best$y
    dt$hm_pos[i] <- best$pos
    placed_bboxes <- c(placed_bboxes, list(best$bb))
  }

  attr(dt, "heatmap_bboxes") <- placed_bboxes
  dt
}

overlay_vec_focus_pdf <- file.path(out_dir, "pyrimidine_metabolism_pathway_overlay_vector_focus.pdf")
overlay_vec_focus_png <- file.path(out_dir, "pyrimidine_metabolism_pathway_overlay_vector_focus.png")
overlay_vec_full_pdf  <- file.path(out_dir, "pyrimidine_metabolism_pathway_overlay_vector_full.pdf")
overlay_vec_full_png  <- file.path(out_dir, "pyrimidine_metabolism_pathway_overlay_vector_full.png")

rxn_dt <- read_kgml_reactions(kgml_path)

draw_overlay_vector <- function(device_fun,
                                mode = c("focus", "full"),
                                fig_w_in = 7.2,
                                pad = 18,
                                met_rect_w = 46,
                                met_rect_h = 12,
                                gene_keep_ids = NULL,
                                met_keep_ids = NULL,
                                show_title = TRUE,
                                show_legend = TRUE,
                                title_text = "Pyrimidine metabolism") {
  mode <- match.arg(mode)

  gene_plot <- copy(gene_nodes_plot)
  gene_plot[, label_draw := label_clean]
  gene_plot[duplicated(label_draw), label_draw := ""]

  met_plot <- copy(compound_nodes)

  if (!is.null(gene_keep_ids)) gene_plot <- gene_plot[id %in% gene_keep_ids]
  if (!is.null(met_keep_ids)) met_plot <- met_plot[id %in% met_keep_ids]

  comp_in_rxn <- unique(c(unlist(rxn_dt$substrate_ids), unlist(rxn_dt$product_ids)))
  comp_in_rxn <- comp_in_rxn[!is.na(comp_in_rxn)]

  if (mode == "focus") {
    # Focus view: crop tightly around nodes that actually carry data (genes + mapped metabolites),
    # then add minimal context compounds within that window (to avoid huge whitespace from remote branches).
    gene_ids_plot <- gene_plot$id
    met_ids_plot <- met_plot$id

    gene_core <- nodes[type == "gene" & id %in% gene_ids_plot]
    comp_core <- nodes[type == "compound" & id %in% met_ids_plot]

    nodes_core <- rbindlist(list(gene_core, comp_core), use.names = TRUE, fill = TRUE)
    if (nrow(nodes_core) == 0) stop("No drawable KGML nodes with mapped data found.")

    xmin0 <- min(nodes_core$x - nodes_core$w / 2, na.rm = TRUE) - pad
    xmax0 <- max(nodes_core$x + nodes_core$w / 2, na.rm = TRUE) + pad
    ymin0 <- min(nodes_core$y - nodes_core$h / 2, na.rm = TRUE) - pad
    ymax0 <- max(nodes_core$y + nodes_core$h / 2, na.rm = TRUE) + pad

    # Add 1-step reaction context compounds, but only keep those that fall inside the focus window
    ctx_comp_ids <- unique(c(
      met_ids_plot,
      unlist(rxn_dt[enzyme_id %in% gene_ids_plot, substrate_ids]),
      unlist(rxn_dt[enzyme_id %in% gene_ids_plot, product_ids])
    ))
    ctx_comp_ids <- ctx_comp_ids[!is.na(ctx_comp_ids)]

    gene_ctx <- gene_core
    compound_ctx <- nodes[type == "compound" & id %in% ctx_comp_ids]
    compound_ctx <- compound_ctx[x >= xmin0 & x <= xmax0 & y >= ymin0 & y <= ymax0]

    nodes_bbox <- rbindlist(list(nodes_core, compound_ctx), use.names = TRUE, fill = TRUE)
  } else {
    # Full view: still avoid drawing/boxing orphan compounds that are not part of any KGML reaction
    gene_ctx <- nodes[type == "gene"]
    compound_ctx <- nodes[type == "compound" & id %in% unique(c(comp_in_rxn, met_plot$id))]
    nodes_bbox <- rbindlist(list(gene_ctx, compound_ctx), use.names = TRUE, fill = TRUE)
  }

  if (nrow(nodes_bbox) == 0) stop("No drawable KGML nodes found.")

  # Base bounds (for greedy placement reference)
  xmin_base <- min(nodes_bbox$x - nodes_bbox$w / 2, na.rm = TRUE) - pad
  xmax_base <- max(nodes_bbox$x + nodes_bbox$w / 2, na.rm = TRUE) + pad
  ymin_base <- min(nodes_bbox$y - nodes_bbox$h / 2, na.rm = TRUE) - pad
  ymax_base <- max(nodes_bbox$y + nodes_bbox$h / 2, na.rm = TRUE) + pad

  # Node bboxes for collision checking (genes + compound anchors)
  node_bbox_list <- list()
  if (nrow(gene_ctx) > 0) {
    node_bbox_list <- c(node_bbox_list, lapply(seq_len(nrow(gene_ctx)), function(i) {
      rect_bbox(gene_ctx$x[i], gene_ctx$y[i], gene_ctx$w[i], gene_ctx$h[i])
    }))
  }
  if (nrow(compound_ctx) > 0) {
    node_bbox_list <- c(node_bbox_list, lapply(seq_len(nrow(compound_ctx)), function(i) {
      rect_bbox(compound_ctx$x[i], compound_ctx$y[i], compound_ctx$w[i], compound_ctx$h[i])
    }))
  }

  # Place metabolite heatmap boxes near their KEGG compound anchors to avoid overlaps.
  if (nrow(met_plot) > 0) {
    met_plot[, `:=`(
      hm_w = pmax(w, met_rect_w),
      hm_h = pmax(h, met_rect_h)
    )]
    met_plot <- place_heatmap_boxes_greedy(
      met_plot,
      anchor_w = met_plot$w, anchor_h = met_plot$h,
      rect_w = met_plot$hm_w, rect_h = met_plot$hm_h,
      node_bbox_list = node_bbox_list,
      xmin_ref = xmin_base, xmax_ref = xmax_base, ymin_ref = ymin_base, ymax_ref = ymax_base,
      gap_seq = c(0, 3, 5, 7, 9, 12, 15),
      order_pref = c("right", "left", "bottom", "top", "top-right", "bottom-right", "top-left", "bottom-left", "center"),
      pad_native = 1.2
    )
  }

  hm_bbox_list <- if (nrow(met_plot) > 0) attr(met_plot, "heatmap_bboxes") else list()

  # Final bounds include heatmap boxes to minimize whitespace without clipping
  all_bbox <- c(node_bbox_list, hm_bbox_list)
  xmin <- min(sapply(all_bbox, function(bb) bb$xmin), na.rm = TRUE) - pad
  xmax <- max(sapply(all_bbox, function(bb) bb$xmax), na.rm = TRUE) + pad
  ymin <- min(sapply(all_bbox, function(bb) bb$ymin), na.rm = TRUE) - pad
  ymax <- max(sapply(all_bbox, function(bb) bb$ymax), na.rm = TRUE) + pad

  x_range <- xmax - xmin
  y_range <- ymax - ymin
  fig_h_in <- fig_w_in * (y_range / x_range)

  device_fun(fig_w_in, fig_h_in)
  grid.newpage()

  # Invert y-axis to match KGML (y increases downward)
  pushViewport(viewport(xscale = c(xmin, xmax), yscale = c(ymax, ymin)))

  # Background
  grid.rect(x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            width = unit(1, "npc"), height = unit(1, "npc"),
            gp = gpar(fill = "white", col = NA))

  # Draw reaction edges (behind nodes)
  if (nrow(rxn_dt) > 0) {
    edge_gp <- gpar(col = "grey80", lwd = 0.35)
    arrow_irrev <- arrow(length = unit(1.2, "mm"), type = "closed")
    arrow_rev <- arrow(length = unit(1.2, "mm"), type = "closed", ends = "both")

    rxn_use <- copy(rxn_dt)
    if (mode == "focus") {
      rxn_use <- rxn_use[enzyme_id %in% gene_ctx$id]
    }

    for (i in seq_len(nrow(rxn_use))) {
      r <- rxn_use[i]
      enz <- gene_ctx[id == r$enzyme_id]
      if (nrow(enz) == 0) next

      subs <- compound_ctx[id %in% unlist(r$substrate_ids)]
      prods <- compound_ctx[id %in% unlist(r$product_ids)]
      if (nrow(subs) == 0 && nrow(prods) == 0) next

      # Substrate -> enzyme (no arrow to reduce clutter)
      if (nrow(subs) > 0) {
        for (k in seq_len(nrow(subs))) {
          grid.segments(
            x0 = unit(subs$x[k], "native"),
            y0 = unit(subs$y[k], "native"),
            x1 = unit(enz$x[1], "native"),
            y1 = unit(enz$y[1], "native"),
            gp = edge_gp
          )
        }
      }

      # Enzyme -> product (arrow indicates direction)
      if (nrow(prods) > 0) {
        for (k in seq_len(nrow(prods))) {
          grid.segments(
            x0 = unit(enz$x[1], "native"),
            y0 = unit(enz$y[1], "native"),
            x1 = unit(prods$x[k], "native"),
            y1 = unit(prods$y[k], "native"),
            gp = edge_gp,
            arrow = if (isTRUE(r$reversible)) arrow_rev else arrow_irrev
          )
        }
      }
    }
  }

  # Draw compound anchors as small circles (context)
  if (nrow(compound_ctx) > 0) {
    for (k in seq_len(nrow(compound_ctx))) {
      node <- compound_ctx[k]
      fill0 <- if (!is.na(node$bgcolor) && node$bgcolor != "") node$bgcolor else "white"
      col0  <- if (!is.na(node$fgcolor) && node$fgcolor != "") node$fgcolor else "grey70"
      draw_kegg_shape(node, fill = fill0, col = col0, lwd = 0.4, alpha = 1.0)
    }
  }

  # Draw gene boxes (outline) for context (KEGG shapes); keep background white to avoid confusion with red/green scale
  if (nrow(gene_ctx) > 0) {
    for (k in seq_len(nrow(gene_ctx))) {
      node <- gene_ctx[k]
      draw_kegg_shape(node, fill = "white", col = "grey85", lwd = 0.35, alpha = 1.0)
    }
  }

  # Draw gene mini-heatmaps (2x3) for genes with data
  if (nrow(gene_plot) > 0) {
    gene_label_gp <- gpar(fontsize = 6.8, fontface = "bold", col = "grey10", fontfamily = "sans")
    gene_outline_gp_col <- function(node) {
      fg <- as.character(node$fgcolor)[1]
      if (!is.na(fg) && fg != "") return(fg)
      "grey25"
    }

    for (k in seq_len(nrow(gene_plot))) {
      g <- gene_plot[k]
      mat <- clip_mat(gene_fc_lookup(g$gene_syms[[1]]))
      stars <- gene_star_lookup(g$gene_syms[[1]])
      draw_mini_heatmap_native(mat, g$x, g$y, g$w, g$h, star_mat = stars,
                               nrow_cells = 2, ncol_cells = 3, col_fun = cfun,
                               cell_border_lwd = 0.3, star_fontsize = 6.0)

      # Outline on top (match KEGG shape)
      node0 <- gene_ctx[id == g$id][1]
      if (nrow(node0) == 0) node0 <- g
      draw_kegg_shape(node0, fill = NA, col = gene_outline_gp_col(node0), lwd = 0.7, alpha = 1.0)
    }
  }

  # Draw metabolite mini-heatmaps (1x3), larger than anchor circles
  if (nrow(met_plot) > 0) {
    met_outline_gp <- gpar(col = "grey25", lwd = 0.55, fill = NA)

    for (k in seq_len(nrow(met_plot))) {
      m <- met_plot[k]
      rect_w <- max(m$hm_w, met_rect_w)
      rect_h <- max(m$hm_h, met_rect_h)

      mat <- clip_mat(met_fc_lookup(m$met_label))
      stars <- met_star_lookup(m$met_label)
      draw_mini_heatmap_native(mat, m$hm_x, m$hm_y, rect_w, rect_h, star_mat = stars,
                               nrow_cells = 1, ncol_cells = 3, col_fun = cfun,
                               cell_border_lwd = 0.3, star_fontsize = 6.0)

      # Outline for the metabolite heatmap box
      grid.rect(
        x = unit(m$hm_x, "native"),
        y = unit(m$hm_y, "native"),
        width = unit(rect_w, "native"),
        height = unit(rect_h, "native"),
        gp = met_outline_gp
      )

      # Connect anchor circle to heatmap box (if offset)
      if (is.finite(m$hm_x) && is.finite(m$hm_y) && ((abs(m$hm_x - m$x) + abs(m$hm_y - m$y)) > 1e-6)) {
        grid.segments(
          x0 = unit(m$x, "native"),
          y0 = unit(m$y, "native"),
          x1 = unit(m$hm_x, "native"),
          y1 = unit(m$hm_y, "native"),
          gp = gpar(col = "grey75", lwd = 0.35)
        )
      }
    }
  }

  # Label placement (greedy; avoid overlaps with nodes + heatmap boxes)
  all_node_bboxes <- c(node_bbox_list, hm_bbox_list)

  # Gene labels
  gene_labels <- gene_plot
  if (nrow(gene_labels) > 0) {
    gene_labels <- place_labels_greedy(
      gene_labels,
      label_col = "label_draw",
      node_bbox_list = all_node_bboxes,
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      gp = gpar(fontsize = 6.4, fontface = "bold", fontfamily = "sans"),
      gap_seq = c(2, 4, 6, 8, 10, 12, 15, 18, 22, 28, 34, 40),
      order_pref = c("left", "right", "top", "bottom"),
      pad_native = 1.2,
      placed_bboxes_init = NULL,
      attach_bboxes = TRUE,
      allow_overlap_own_node = FALSE
    )
  }
  gene_label_bboxes <- if (nrow(gene_labels) > 0) attr(gene_labels, "label_bboxes") else list()

  if (nrow(gene_labels) > 0) {
    for (k in seq_len(nrow(gene_labels))) {
      g <- gene_labels[k]
      if (is.na(g$label_draw) || g$label_draw == "" || !is.finite(g$label_x) || !is.finite(g$label_y)) next
      grid.text(g$label_draw, x = unit(g$label_x, "native"), y = unit(g$label_y, "native"),
                just = "center", gp = gpar(fontsize = 6.4, fontface = "bold", col = "grey10", fontfamily = "sans"))
    }
  }

  # Metabolite labels
  met_labels <- met_plot
  if (nrow(met_labels) > 0) {
    met_labels[, met_lab_draw := {
      met_id <- str_extract(name, "C\\d{5}")
      ifelse(!is.na(met_id) & met_id %in% names(kegg_compounds), unname(kegg_compounds[met_id]), met_label)
    }]
    # Place labels relative to the *heatmap box* (not the small anchor circle) to avoid covering the heatmap.
    met_labels[, `:=`(x = hm_x, y = hm_y, w = hm_w, h = hm_h)]
    met_labels <- place_labels_greedy(
      met_labels,
      label_col = "met_lab_draw",
      node_bbox_list = all_node_bboxes,
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      gp = gpar(fontsize = 6.0, fontfamily = "sans"),
      gap_seq = c(2, 4, 6, 8, 10, 12, 15, 18, 22, 28, 34),
      order_pref = c("bottom", "top", "right", "left"),
      pad_native = 1.2,
      placed_bboxes_init = gene_label_bboxes,
      attach_bboxes = TRUE,
      allow_overlap_own_node = FALSE
    )
  }

  if (nrow(met_labels) > 0) {
    for (k in seq_len(nrow(met_labels))) {
      m <- met_labels[k]
      if (is.na(m$met_lab_draw) || m$met_lab_draw == "" || !is.finite(m$label_x) || !is.finite(m$label_y)) next
      grid.text(m$met_lab_draw, x = unit(m$label_x, "native"), y = unit(m$label_y, "native"),
                just = "center", gp = gpar(fontsize = 6.0, col = "grey10", fontfamily = "sans"))
    }
  }

  if (isTRUE(show_title)) {
    grid.text(
      title_text,
      x = unit(0.02, "npc"),
      y = unit(0.985, "npc"),
      just = c("left", "top"),
      gp = gpar(fontsize = 10, fontface = "bold", col = "grey10", fontfamily = "sans")
    )
  }

  if (isTRUE(show_legend)) {
    # Minimal legend (top-right). For publication, a separate legend-only file is also generated.
    lgd <- Legend(
      title = "log2FC",
      col_fun = cfun,
      at = c(-cap_abs, 0, cap_abs),
      labels = c(paste0("-", cap_abs), "0", paste0("+", cap_abs)),
      direction = "horizontal",
      legend_width = unit(28, "mm"),
      title_position = "topcenter"
    )
    pushViewport(viewport(x = unit(0.985, "npc"), y = unit(0.985, "npc"), just = c("right", "top")))
    ComplexHeatmap::draw(lgd, just = c("right", "top"))
    popViewport()
  }

  popViewport()
  dev.off()
}

draw_overlay_vector(function(w, h) pdf(overlay_vec_focus_pdf, width = w, height = h, useDingbats = FALSE),
                    mode = "focus",
                    show_title = FALSE,
                    show_legend = FALSE)
draw_overlay_vector(function(w, h) png(overlay_vec_focus_png, width = w, height = h, units = "in", res = 600),
                    mode = "focus",
                    show_title = FALSE,
                    show_legend = FALSE)
draw_overlay_vector(function(w, h) pdf(overlay_vec_full_pdf, width = w, height = h, useDingbats = FALSE),
                    mode = "full",
                    show_title = FALSE,
                    show_legend = FALSE)
draw_overlay_vector(function(w, h) png(overlay_vec_full_png, width = w, height = h, units = "in", res = 600),
                    mode = "full",
                    show_title = FALSE,
                    show_legend = FALSE)

message("Saved: ", overlay_vec_focus_pdf)
message("Saved: ", overlay_vec_focus_png)
message("Saved: ", overlay_vec_full_pdf)
message("Saved: ", overlay_vec_full_png)

# Optional: split "focus" into two spatial clusters to reduce internal whitespace (main module vs distant branch).
split_focus_clusters <- function(gene_dt, met_dt) {
  pts <- rbindlist(list(
    gene_dt[, .(id, type = "gene", x, y)],
    met_dt[, .(id, type = "compound", x, y)]
  ), use.names = TRUE, fill = TRUE)
  pts <- pts[is.finite(x) & is.finite(y)]

  out <- list(
    main_gene_ids = unique(gene_dt$id),
    branch_gene_ids = character(),
    main_met_ids = unique(met_dt$id),
    branch_met_ids = character()
  )

  if (nrow(pts) < 6) return(out)

  set.seed(1)
  km <- tryCatch(stats::kmeans(scale(pts[, .(x, y)]), centers = 2, nstart = 10), error = function(e) NULL)
  if (is.null(km)) return(out)

  pts[, cluster := km$cluster]
  medx <- pts[, .(medx = median(x, na.rm = TRUE)), by = cluster][order(medx)]
  main_cluster <- medx$cluster[1]

  out$main_gene_ids <- unique(pts[type == "gene" & cluster == main_cluster]$id)
  out$branch_gene_ids <- unique(pts[type == "gene" & cluster != main_cluster]$id)
  out$main_met_ids <- unique(pts[type == "compound" & cluster == main_cluster]$id)
  out$branch_met_ids <- unique(pts[type == "compound" & cluster != main_cluster]$id)

  # If the split degenerates, fall back to single focus.
  if (length(out$main_gene_ids) + length(out$main_met_ids) == 0) return(out)
  if (length(out$branch_gene_ids) + length(out$branch_met_ids) == 0) return(out)

  out
}

clusters <- split_focus_clusters(gene_nodes_plot, compound_nodes)

overlay_vec_main_pdf <- file.path(out_dir, "pyrimidine_metabolism_pathway_overlay_vector_main.pdf")
overlay_vec_main_png <- file.path(out_dir, "pyrimidine_metabolism_pathway_overlay_vector_main.png")
draw_overlay_vector(function(w, h) pdf(overlay_vec_main_pdf, width = w, height = h, useDingbats = FALSE),
                    mode = "focus",
                    fig_w_in = 10.5,
                    gene_keep_ids = clusters$main_gene_ids,
                    met_keep_ids = clusters$main_met_ids,
                    show_title = FALSE,
                    show_legend = FALSE)
draw_overlay_vector(function(w, h) png(overlay_vec_main_png, width = w, height = h, units = "in", res = 600),
                    mode = "focus",
                    fig_w_in = 10.5,
                    gene_keep_ids = clusters$main_gene_ids,
                    met_keep_ids = clusters$main_met_ids,
                    show_title = FALSE,
                    show_legend = FALSE)
message("Saved: ", overlay_vec_main_pdf)
message("Saved: ", overlay_vec_main_png)

if (length(clusters$branch_gene_ids) + length(clusters$branch_met_ids) > 0) {
  overlay_vec_branch_pdf <- file.path(out_dir, "pyrimidine_metabolism_pathway_overlay_vector_branch.pdf")
  overlay_vec_branch_png <- file.path(out_dir, "pyrimidine_metabolism_pathway_overlay_vector_branch.png")

  draw_overlay_vector(function(w, h) pdf(overlay_vec_branch_pdf, width = w, height = h, useDingbats = FALSE),
                      mode = "focus",
                      fig_w_in = 10.5,
                      gene_keep_ids = clusters$branch_gene_ids,
                      met_keep_ids = clusters$branch_met_ids,
                      show_title = FALSE,
                      show_legend = FALSE)
  draw_overlay_vector(function(w, h) png(overlay_vec_branch_png, width = w, height = h, units = "in", res = 600),
                      mode = "focus",
                      fig_w_in = 10.5,
                      gene_keep_ids = clusters$branch_gene_ids,
                      met_keep_ids = clusters$branch_met_ids,
                      show_title = FALSE,
                      show_legend = FALSE)
  message("Saved: ", overlay_vec_branch_pdf)
  message("Saved: ", overlay_vec_branch_png)
}

# Legend-only panel (for figure assembly)
draw_overlay_legend_only <- function(device_fun, fig_w_in = 3.6, fig_h_in = 1.45) {
  device_fun(fig_w_in, fig_h_in)
  grid.newpage()
  grid.rect(x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            width = unit(1, "npc"), height = unit(1, "npc"),
            gp = gpar(fill = "white", col = NA))

  lgd <- Legend(
    title = "log2FC",
    col_fun = cfun,
    at = c(-cap_abs, 0, cap_abs),
    labels = c(paste0("-", cap_abs), "0", paste0("+", cap_abs)),
    direction = "horizontal",
    legend_width = unit(46, "mm"),
    title_position = "topcenter",
    title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "sans"),
    labels_gp = gpar(fontsize = 9, fontfamily = "sans")
  )
  pushViewport(viewport(x = unit(0.5, "npc"), y = unit(0.93, "npc"), just = c("center", "top")))
  ComplexHeatmap::draw(lgd, just = c("center", "top"))
  popViewport()

  # 2x3 layout key: rows = RNA/Prot; cols = sh/C/OE
  key_x0 <- 0.18
  key_y0 <- 0.17
  key_w <- 0.64
  key_h <- 0.42

  grid.rect(x = unit(key_x0 + key_w / 2, "npc"), y = unit(key_y0 + key_h / 2, "npc"),
            width = unit(key_w, "npc"), height = unit(key_h, "npc"),
            gp = gpar(fill = adjustcolor("white", alpha.f = 0.98), col = "grey80", lwd = 0.35))

  cell_w <- key_w / 3
  cell_h <- key_h / 2
  for (ii in 1:2) {
    for (jj in 1:3) {
      grid.rect(
        x = unit(key_x0 + (jj - 0.5) * cell_w, "npc"),
        y = unit(key_y0 + (2 - ii + 0.5) * cell_h, "npc"),
        width = unit(cell_w, "npc"),
        height = unit(cell_h, "npc"),
        gp = gpar(fill = "white", col = "grey55", lwd = 0.3)
      )
    }
  }
  grid.text("RNA", x = unit(key_x0 - 0.02, "npc"), y = unit(key_y0 + cell_h * 1.5, "npc"),
            just = "right", gp = gpar(fontsize = 9, fontfamily = "sans", col = "grey20"))
  grid.text("Prot", x = unit(key_x0 - 0.02, "npc"), y = unit(key_y0 + cell_h * 0.5, "npc"),
            just = "right", gp = gpar(fontsize = 9, fontfamily = "sans", col = "grey20"))
  for (jj in 1:3) {
    grid.text(GENOTYPES[jj], x = unit(key_x0 + (jj - 0.5) * cell_w, "npc"),
              y = unit(key_y0 + key_h + 0.02, "npc"), just = "bottom",
              gp = gpar(fontsize = 9, fontface = "bold", fontfamily = "sans", col = "grey20"))
  }

  grid.text("H2O2 vs Normal", x = unit(0.5, "npc"), y = unit(0.64, "npc"),
            just = "center", gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "sans", col = "grey10"))
  grid.text("* p<0.05    ** p<0.01    *** p<0.001", x = unit(0.5, "npc"), y = unit(0.08, "npc"),
            just = "center", gp = gpar(fontsize = 9, fontfamily = "sans", col = "grey20"))

  dev.off()
}

overlay_legend_pdf <- file.path(out_dir, "pyrimidine_metabolism_pathway_overlay_legend.pdf")
overlay_legend_png <- file.path(out_dir, "pyrimidine_metabolism_pathway_overlay_legend.png")
draw_overlay_legend_only(function(w, h) pdf(overlay_legend_pdf, width = w, height = h, useDingbats = FALSE))
draw_overlay_legend_only(function(w, h) png(overlay_legend_png, width = w, height = h, units = "in", res = 600))
message("Saved: ", overlay_legend_pdf)
message("Saved: ", overlay_legend_png)
