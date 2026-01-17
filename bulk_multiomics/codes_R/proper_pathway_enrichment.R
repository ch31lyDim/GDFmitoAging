#!/usr/bin/env Rscript
# Proper Pathway Enrichment using KEGG_MapID from data
# Author: Automated Analysis
# Date: 2024

cat("==========================================================================\n")
cat("PROPER PATHWAY ENRICHMENT ANALYSIS - Using KEGG_MapID\n")
cat("==========================================================================\n\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(writexl)
  library(tidyr)
})

# Load results
results_dir <- "metabo_results_GDF11KO"
liver_results <- read.csv(file.path(results_dir, "liver/liver_differential_metabolites.csv"))
blood_results <- read.csv(file.path(results_dir, "blood/blood_differential_metabolites.csv"))
liver_sig <- liver_results[liver_results$significant == "Significant", ]
blood_sig <- blood_results[blood_results$significant == "Significant", ]

cat("Loaded results:\n")
cat(sprintf("  Liver: %d total, %d significant\n", nrow(liver_results), nrow(liver_sig)))
cat(sprintf("  Blood: %d total, %d significant\n", nrow(blood_results), nrow(blood_sig)))

# ==========================================================================
# Extract Metabolite-to-Pathway Mapping from KEGG_MapID
# ==========================================================================

cat("\nExtracting metabolite-to-pathway mappings...\n")

extract_pathway_mapping <- function(data) {
  # Filter rows with valid KEGG_MapID
  data_with_map <- data %>%
    filter(!is.na(KEGG_MapID) & KEGG_MapID != "-" & KEGG_MapID != "")
  
  cat(sprintf("  Metabolites with pathway info: %d / %d (%.1f%%)\n", 
              nrow(data_with_map), nrow(data), 100*nrow(data_with_map)/nrow(data)))
  
  # Expand KEGG_MapID (split by semicolon)
  pathway_map <- data_with_map %>%
    select(Compound_ID, Name, KEGG_MapID, KEGG_ID, HMDB_ID) %>%
    separate_rows(KEGG_MapID, sep = ";") %>%
    mutate(KEGG_MapID = trimws(KEGG_MapID))
  
  return(pathway_map)
}

liver_pathway_map <- extract_pathway_mapping(liver_results)
blood_pathway_map <- extract_pathway_mapping(blood_results)

cat(sprintf("\n  Liver: %d metabolite-pathway pairs\n", nrow(liver_pathway_map)))
cat(sprintf("  Blood: %d metabolite-pathway pairs\n", nrow(blood_pathway_map)))

# Get pathway statistics
liver_pathways <- liver_pathway_map %>% 
  group_by(KEGG_MapID) %>% 
  summarise(n_metabolites = n_distinct(Compound_ID)) %>%
  arrange(desc(n_metabolites))

cat(sprintf("\n  Unique pathways in liver: %d\n", nrow(liver_pathways)))
cat("\n  Top 10 pathways by metabolite count:\n")
print(head(liver_pathways, 10))

# ==========================================================================
# Load KEGG Pathway Names
# ==========================================================================

cat("\n\nLoading KEGG pathway database for names...\n")
kegg_db <- read.csv("kn_mouse_meta/analysis/scripts/kegg_pathway.csv")

# Extract pathway names (format: hsa00010 → map00010)
pathway_names <- kegg_db %>%
  mutate(KEGG_MapID = gsub("^hsa", "map", id)) %>%
  select(KEGG_MapID, name)

cat(sprintf("  Loaded %d pathway names\n", nrow(pathway_names)))

# ==========================================================================
# ORA Analysis using KEGG_MapID
# ==========================================================================

cat("\n\nPerforming ORA using KEGG_MapID...\n")

perform_kegg_ora <- function(sig_data, all_data, pathway_map, pathway_names, min_size = 3) {
  
  # Get sig and all compound IDs
  sig_compounds <- sig_data$Compound_ID
  all_compounds <- all_data$Compound_ID
  
  cat(sprintf("  Significant metabolites: %d\n", length(sig_compounds)))
  cat(sprintf("  Total metabolites: %d\n", length(all_compounds)))
  
  # Get pathways for significant and all metabolites
  sig_pathways <- pathway_map %>%
    filter(Compound_ID %in% sig_compounds) %>%
    group_by(KEGG_MapID) %>%
    summarise(sig_count = n_distinct(Compound_ID))
  
  all_pathways <- pathway_map %>%
    filter(Compound_ID %in% all_compounds) %>%
    group_by(KEGG_MapID) %>%
    summarise(all_count = n_distinct(Compound_ID))
  
  # Merge
  pathway_counts <- full_join(sig_pathways, all_pathways, by = "KEGG_MapID") %>%
    mutate(sig_count = ifelse(is.na(sig_count), 0, sig_count),
           all_count = ifelse(is.na(all_count), 0, all_count))
  
  # Filter by minimum size
  pathway_counts <- pathway_counts %>%
    filter(all_count >= min_size & sig_count >= min_size)
  
  cat(sprintf("  Pathways with ≥%d metabolites: %d\n", min_size, nrow(pathway_counts)))
  
  # Perform Fisher's exact test for each pathway
  results <- data.frame()
  
  n_sig_total <- length(unique(sig_compounds))
  n_all_total <- length(unique(all_compounds))
  
  for (i in 1:nrow(pathway_counts)) {
    pathway_id <- pathway_counts$KEGG_MapID[i]
    sig_in <- pathway_counts$sig_count[i]
    all_in <- pathway_counts$all_count[i]
    
    # Contingency table
    cont_table <- matrix(c(
      sig_in,                    # sig & in pathway
      n_sig_total - sig_in,      # sig & not in pathway
      all_in - sig_in,           # not sig & in pathway
      n_all_total - n_sig_total - (all_in - sig_in)  # not sig & not in pathway
    ), nrow = 2)
    
    # Fisher's exact test
    test <- fisher.test(cont_table, alternative = "greater")
    
    # Calculate fold enrichment
    expected <- (all_in / n_all_total) * n_sig_total
    fold_enrichment <- sig_in / expected
    
    results <- rbind(results, data.frame(
      KEGG_MapID = pathway_id,
      Significant_Count = sig_in,
      Pathway_Size = all_in,
      Total_Significant = n_sig_total,
      Expected = expected,
      Fold_Enrichment = fold_enrichment,
      Pvalue = test$p.value,
      Odds_Ratio = as.numeric(test$estimate),
      stringsAsFactors = FALSE
    ))
  }
  
  # Add FDR correction
  results$FDR <- p.adjust(results$Pvalue, method = "BH")
  
  # Add pathway names
  results <- left_join(results, pathway_names, by = "KEGG_MapID")
  
  # Sort by p-value
  results <- results %>%
    arrange(Pvalue) %>%
    select(KEGG_MapID, Pathway_Name = name, Significant_Count, Pathway_Size,
           Fold_Enrichment, Pvalue, FDR, Odds_Ratio, Expected)
  
  return(results)
}

# Perform ORA for liver
cat("\n=== LIVER ORA ===\n")
liver_ora <- perform_kegg_ora(
  sig_data = liver_sig,
  all_data = liver_results,
  pathway_map = liver_pathway_map,
  pathway_names = pathway_names,
  min_size = 3
)

write.csv(liver_ora,
          file.path(results_dir, "liver/pathway_enrichment", "liver_ORA_KEGG_proper.csv"),
          row.names = FALSE)

if (nrow(liver_ora) > 0) {
  cat(sprintf("\n  Total pathways tested: %d\n", nrow(liver_ora)))
  cat(sprintf("  Significant (p<0.05): %d\n", sum(liver_ora$Pvalue < 0.05)))
  cat(sprintf("  FDR significant (FDR<0.05): %d\n", sum(liver_ora$FDR < 0.05)))
  cat(sprintf("  FDR significant (FDR<0.1): %d\n", sum(liver_ora$FDR < 0.1)))
  cat(sprintf("  FDR significant (FDR<0.25): %d\n", sum(liver_ora$FDR < 0.25)))
  
  if (sum(liver_ora$Pvalue < 0.05) > 0) {
    cat("\n  Top enriched pathways (p<0.05):\n")
    top_liver <- liver_ora[liver_ora$Pvalue < 0.05, ]
    print(top_liver[, c("Pathway_Name", "Significant_Count", "Pathway_Size", 
                        "Fold_Enrichment", "Pvalue", "FDR")])
  }
}

# Perform ORA for blood
cat("\n\n=== BLOOD ORA ===\n")
blood_ora <- perform_kegg_ora(
  sig_data = blood_sig,
  all_data = blood_results,
  pathway_map = blood_pathway_map,
  pathway_names = pathway_names,
  min_size = 2
)

write.csv(blood_ora,
          file.path(results_dir, "blood/pathway_enrichment", "blood_ORA_KEGG_proper.csv"),
          row.names = FALSE)

if (nrow(blood_ora) > 0) {
  cat(sprintf("\n  Total pathways tested: %d\n", nrow(blood_ora)))
  cat(sprintf("  Significant (p<0.05): %d\n", sum(blood_ora$Pvalue < 0.05)))
  cat(sprintf("  FDR significant (FDR<0.1): %d\n", sum(blood_ora$FDR < 0.1)))
  
  if (sum(blood_ora$Pvalue < 0.05) > 0) {
    cat("\n  Top enriched pathways (p<0.05):\n")
    top_blood <- blood_ora[blood_ora$Pvalue < 0.05, ]
    print(top_blood[, c("Pathway_Name", "Significant_Count", "Pathway_Size",
                        "Fold_Enrichment", "Pvalue", "FDR")])
  } else if (sum(blood_ora$Pvalue < 0.1) > 0) {
    cat("\n  Nominally enriched pathways (p<0.1):\n")
    top_blood <- blood_ora[blood_ora$Pvalue < 0.1, ]
    print(top_blood[, c("Pathway_Name", "Significant_Count", "Pathway_Size",
                        "Fold_Enrichment", "Pvalue", "FDR")])
  }
}

# ==========================================================================
# Visualization
# ==========================================================================

cat("\n\nGenerating enrichment plots...\n")

create_ora_plot <- function(ora_results, title, output_path, top_n = 20, p_cutoff = 0.1) {
  
  # Filter and prepare data
  plot_data <- ora_results %>%
    filter(Pvalue < p_cutoff) %>%
    head(top_n) %>%
    mutate(neglog10p = -log10(Pvalue),
           FDR_sig = FDR < 0.05,
           Pathway_Name = ifelse(is.na(Pathway_Name), KEGG_MapID, Pathway_Name),
           Pathway_Name = factor(Pathway_Name, levels = rev(Pathway_Name)))
  
  if (nrow(plot_data) == 0) {
    cat(sprintf("  No pathways with p<%s for %s\n", p_cutoff, basename(output_path)))
    return(NULL)
  }
  
  p <- ggplot(plot_data, aes(x = Pathway_Name, y = neglog10p)) +
    geom_bar(aes(fill = FDR_sig), stat = "identity") +
    scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "#E31A1C"),
                     labels = c("FALSE" = "p<0.05", "TRUE" = "FDR<0.05")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) +
    coord_flip() +
    labs(title = title,
         subtitle = sprintf("Top %d enriched KEGG pathways", nrow(plot_data)),
         x = "KEGG Pathway",
         y = "-log10(p-value)",
         fill = "Significance") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 9))
  
  ggsave(output_path, p, width = 12, height = max(6, nrow(plot_data) * 0.35), dpi = 300)
  cat(sprintf("  Saved: %s\n", basename(output_path)))
}

if (nrow(liver_ora) > 0 && sum(liver_ora$Pvalue < 0.1) > 0) {
  create_ora_plot(liver_ora, 
                 "Liver - KEGG Pathway Enrichment (ORA)",
                 file.path(results_dir, "liver/figures", "liver_pathway_ORA_proper.png"),
                 top_n = 30, p_cutoff = 0.1)
}

if (nrow(blood_ora) > 0 && sum(blood_ora$Pvalue < 0.1) > 0) {
  create_ora_plot(blood_ora,
                 "Blood - KEGG Pathway Enrichment (ORA)",
                 file.path(results_dir, "blood/figures", "blood_pathway_ORA_proper.png"),
                 top_n = 20, p_cutoff = 0.1)
}

# ==========================================================================
# Summary
# ==========================================================================

cat("\n\n==========================================================================\n")
cat("PROPER PATHWAY ENRICHMENT COMPLETE\n")
cat("==========================================================================\n\n")

cat("Results Summary:\n")
cat(sprintf("  Liver ORA: %d pathways tested, %d significant (p<0.05), %d FDR<0.1\n",
            nrow(liver_ora), sum(liver_ora$Pvalue < 0.05), sum(liver_ora$FDR < 0.1)))
cat(sprintf("  Blood ORA: %d pathways tested, %d significant (p<0.05), %d FDR<0.1\n",
            nrow(blood_ora), sum(blood_ora$Pvalue < 0.05), sum(blood_ora$FDR < 0.1)))

cat("\nOutput files:\n")
cat("  - liver/pathway_enrichment/liver_ORA_KEGG_proper.csv\n")
cat("  - blood/pathway_enrichment/blood_ORA_KEGG_proper.csv\n")
cat("  - liver/figures/liver_pathway_ORA_proper.png\n")
cat("  - blood/figures/blood_pathway_ORA_proper.png\n")
cat("\n")
