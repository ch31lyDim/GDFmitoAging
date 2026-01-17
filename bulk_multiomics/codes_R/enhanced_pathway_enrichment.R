#!/usr/bin/env Rscript
# Enhanced Pathway Enrichment with Complete Database Annotations
# Part 1: Enhanced KEGG annotations using KEGG REST API
# Part 2: Alternative database enrichment (HMDB, Lipidmaps, SMPDB, Reactome)

library(dplyr)
library(httr)
library(jsonlite)

# Set up directories
base_dir <- "/data1/dyh/GDF11/metabo_results_GDF11KO"
output_dir <- file.path(base_dir, "enhanced_enrichment")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Enhanced Pathway Enrichment Analysis ===\n\n")

# ============================================================================
# PART 1: ENHANCED KEGG PATHWAY ANNOTATIONS
# ============================================================================

cat("PART 1: Downloading complete KEGG pathway information...\n")

# Function to query KEGG REST API with rate limiting
query_kegg_api <- function(endpoint, max_retries = 3) {
  url <- paste0("https://rest.kegg.jp/", endpoint)
  
  for (i in 1:max_retries) {
    tryCatch({
      Sys.sleep(0.5)  # Rate limiting
      response <- GET(url, timeout(30))
      
      if (status_code(response) == 200) {
        return(content(response, as = "text", encoding = "UTF-8"))
      } else {
        cat(sprintf("  Warning: Status %d for %s\n", status_code(response), endpoint))
      }
    }, error = function(e) {
      cat(sprintf("  Error on attempt %d: %s\n", i, e$message))
      if (i < max_retries) Sys.sleep(2)
    })
  }
  
  return(NULL)
}

# Get list of all KEGG pathways
cat("  Fetching pathway list...\n")
pathway_list_text <- query_kegg_api("list/pathway")

if (!is.null(pathway_list_text)) {
  pathway_lines <- strsplit(pathway_list_text, "\n")[[1]]
  pathway_lines <- pathway_lines[pathway_lines != ""]
  
  kegg_pathways <- data.frame(
    pathway_id = sapply(strsplit(pathway_lines, "\t"), `[`, 1),
    pathway_name = sapply(strsplit(pathway_lines, "\t"), `[`, 2),
    stringsAsFactors = FALSE
  )
  
  # Clean pathway IDs to match format in data (map format)
  kegg_pathways$map_id <- sub("^path:", "map", kegg_pathways$pathway_id)
  
  cat(sprintf("  Downloaded %d KEGG pathways\n", nrow(kegg_pathways)))
  
  # Get pathway hierarchies
  cat("  Fetching pathway hierarchies...\n")
  hierarchy_text <- query_kegg_api("get/br:br08901")
  
  if (!is.null(hierarchy_text)) {
    writeLines(hierarchy_text, file.path(output_dir, "kegg_pathway_hierarchy.txt"))
    cat("  Saved pathway hierarchy\n")
  }
  
  # Save complete pathway list
  write.csv(kegg_pathways, 
            file.path(output_dir, "kegg_complete_pathways.csv"),
            row.names = FALSE)
  cat("  Saved complete pathway list\n\n")
} else {
  cat("  Warning: Could not download KEGG pathway list\n\n")
  kegg_pathways <- NULL
}

# ============================================================================
# PART 2: ALTERNATIVE DATABASE ENRICHMENT
# ============================================================================

cat("PART 2: Performing alternative database enrichment...\n\n")

# Load significant metabolites
liver_sig <- read.csv(file.path(base_dir, "liver/liver_significant_metabolites.csv"))
blood_sig <- read.csv(file.path(base_dir, "blood/blood_significant_metabolites.csv"))

cat(sprintf("Loaded %d liver and %d blood significant metabolites\n\n", 
            nrow(liver_sig), nrow(blood_sig)))

# ----------------------------------------------------------------------------
# 2.1: HMDB Class Enrichment
# ----------------------------------------------------------------------------

cat("2.1: HMDB Class Enrichment Analysis\n")

perform_hmdb_enrichment <- function(sig_data, all_data, tissue_name) {
  cat(sprintf("  Processing %s...\n", tissue_name))
  
  # SuperClass enrichment
  sig_superclass <- sig_data %>%
    filter(!is.na(SuperClass.HMDB.) & SuperClass.HMDB. != "" & SuperClass.HMDB. != "-") %>%
    count(SuperClass.HMDB., name = "sig_count")
  
  all_superclass <- all_data %>%
    filter(!is.na(SuperClass.HMDB.) & SuperClass.HMDB. != "" & SuperClass.HMDB. != "-") %>%
    count(SuperClass.HMDB., name = "total_count")
  
  if (nrow(sig_superclass) > 0 && nrow(all_superclass) > 0) {
    superclass_results <- sig_superclass %>%
      left_join(all_superclass, by = "SuperClass.HMDB.") %>%
      mutate(
        sig_total = sum(sig_count),
        all_total = sum(total_count),
        expected = (total_count / all_total) * sig_total,
        fold_enrichment = sig_count / expected,
        pvalue = NA_real_
      )
    
    # Fisher's exact test for each class (use simulation for large tables)
    for (i in 1:nrow(superclass_results)) {
      contingency <- matrix(c(
        superclass_results$sig_count[i],
        superclass_results$sig_total - superclass_results$sig_count[i],
        superclass_results$total_count[i] - superclass_results$sig_count[i],
        superclass_results$all_total - superclass_results$total_count[i] - 
          (superclass_results$sig_total - superclass_results$sig_count[i])
      ), nrow = 2)
      
      # Use simulation for large tables to avoid memory issues
      tryCatch({
        superclass_results$pvalue[i] <- tryCatch(fisher.test(contingency, alternative = "greater")$p.value, error = function(e) fisher.test(contingency, alternative = "greater", simulate.p.value = TRUE, B = 2000)$p.value)
      }, error = function(e) {
        superclass_results$pvalue[i] <<- fisher.test(contingency, alternative = "greater", 
                                                     simulate.p.value = TRUE, B = 2000)$p.value
      })
    }
    
    superclass_results <- superclass_results %>%
      mutate(padj = p.adjust(pvalue, method = "BH")) %>%
      arrange(pvalue)
    
    write.csv(superclass_results,
              file.path(output_dir, sprintf("%s_HMDB_SuperClass_enrichment.csv", tolower(tissue_name))),
              row.names = FALSE)
    
    cat(sprintf("    SuperClass: %d classes tested, %d significant (p<0.05)\n",
                nrow(superclass_results), sum(superclass_results$pvalue < 0.05, na.rm = TRUE)))
  }
  
  # Class enrichment
  sig_class <- sig_data %>%
    filter(!is.na(Class.HMDB.) & Class.HMDB. != "" & Class.HMDB. != "-") %>%
    count(Class.HMDB., name = "sig_count")
  
  all_class <- all_data %>%
    filter(!is.na(Class.HMDB.) & Class.HMDB. != "" & Class.HMDB. != "-") %>%
    count(Class.HMDB., name = "total_count")
  
  if (nrow(sig_class) > 0 && nrow(all_class) > 0) {
    class_results <- sig_class %>%
      left_join(all_class, by = "Class.HMDB.") %>%
      mutate(
        sig_total = sum(sig_count),
        all_total = sum(total_count),
        expected = (total_count / all_total) * sig_total,
        fold_enrichment = sig_count / expected,
        pvalue = NA_real_
      )
    
    for (i in 1:nrow(class_results)) {
      contingency <- matrix(c(
        class_results$sig_count[i],
        class_results$sig_total - class_results$sig_count[i],
        class_results$total_count[i] - class_results$sig_count[i],
        class_results$all_total - class_results$total_count[i] - 
          (class_results$sig_total - class_results$sig_count[i])
      ), nrow = 2)
      
      tryCatch({
        class_results$pvalue[i] <- tryCatch(fisher.test(contingency, alternative = "greater")$p.value, error = function(e) fisher.test(contingency, alternative = "greater", simulate.p.value = TRUE, B = 2000)$p.value)
      }, error = function(e) {
        class_results$pvalue[i] <<- fisher.test(contingency, alternative = "greater",
                                                simulate.p.value = TRUE, B = 2000)$p.value
      })
    }
    
    class_results <- class_results %>%
      mutate(padj = p.adjust(pvalue, method = "BH")) %>%
      arrange(pvalue)
    
    write.csv(class_results,
              file.path(output_dir, sprintf("%s_HMDB_Class_enrichment.csv", tolower(tissue_name))),
              row.names = FALSE)
    
    cat(sprintf("    Class: %d classes tested, %d significant (p<0.05)\n",
                nrow(class_results), sum(class_results$pvalue < 0.05, na.rm = TRUE)))
  }
  
  # SubClass enrichment
  sig_subclass <- sig_data %>%
    filter(!is.na(SubClass.HMDB.) & SubClass.HMDB. != "" & SubClass.HMDB. != "-") %>%
    count(SubClass.HMDB., name = "sig_count")
  
  all_subclass <- all_data %>%
    filter(!is.na(SubClass.HMDB.) & SubClass.HMDB. != "" & SubClass.HMDB. != "-") %>%
    count(SubClass.HMDB., name = "total_count")
  
  if (nrow(sig_subclass) > 0 && nrow(all_subclass) > 0) {
    subclass_results <- sig_subclass %>%
      left_join(all_subclass, by = "SubClass.HMDB.") %>%
      mutate(
        sig_total = sum(sig_count),
        all_total = sum(total_count),
        expected = (total_count / all_total) * sig_total,
        fold_enrichment = sig_count / expected,
        pvalue = NA_real_
      )
    
    for (i in 1:nrow(subclass_results)) {
      contingency <- matrix(c(
        subclass_results$sig_count[i],
        subclass_results$sig_total - subclass_results$sig_count[i],
        subclass_results$total_count[i] - subclass_results$sig_count[i],
        subclass_results$all_total - subclass_results$total_count[i] - 
          (subclass_results$sig_total - subclass_results$sig_count[i])
      ), nrow = 2)
      
      subclass_results$pvalue[i] <- tryCatch(fisher.test(contingency, alternative = "greater")$p.value, error = function(e) fisher.test(contingency, alternative = "greater", simulate.p.value = TRUE, B = 2000)$p.value)
    }
    
    subclass_results <- subclass_results %>%
      mutate(padj = p.adjust(pvalue, method = "BH")) %>%
      arrange(pvalue)
    
    write.csv(subclass_results,
              file.path(output_dir, sprintf("%s_HMDB_SubClass_enrichment.csv", tolower(tissue_name))),
              row.names = FALSE)
    
    cat(sprintf("    SubClass: %d classes tested, %d significant (p<0.05)\n",
                nrow(subclass_results), sum(subclass_results$pvalue < 0.05, na.rm = TRUE)))
  }
  
  cat("\n")
}

# Load all metabolites for background
liver_all <- read.csv(file.path(base_dir, "liver/liver_differential_metabolites.csv"))
blood_all <- read.csv(file.path(base_dir, "blood/blood_differential_metabolites.csv"))

perform_hmdb_enrichment(liver_sig, liver_all, "Liver")
perform_hmdb_enrichment(blood_sig, blood_all, "Blood")

# ----------------------------------------------------------------------------
# 2.2: Lipidmaps Category Enrichment
# ----------------------------------------------------------------------------

cat("2.2: Lipidmaps Category Enrichment Analysis\n")

perform_lipidmaps_enrichment <- function(sig_data, all_data, tissue_name) {
  cat(sprintf("  Processing %s...\n", tissue_name))
  
  # Category enrichment
  sig_category <- sig_data %>%
    filter(!is.na(CATEGORY.Lipidmaps.) & CATEGORY.Lipidmaps. != "" & CATEGORY.Lipidmaps. != "-") %>%
    count(CATEGORY.Lipidmaps., name = "sig_count")
  
  all_category <- all_data %>%
    filter(!is.na(CATEGORY.Lipidmaps.) & CATEGORY.Lipidmaps. != "" & CATEGORY.Lipidmaps. != "-") %>%
    count(CATEGORY.Lipidmaps., name = "total_count")
  
  if (nrow(sig_category) > 0 && nrow(all_category) > 0) {
    category_results <- sig_category %>%
      left_join(all_category, by = "CATEGORY.Lipidmaps.") %>%
      mutate(
        sig_total = sum(sig_count),
        all_total = sum(total_count),
        expected = (total_count / all_total) * sig_total,
        fold_enrichment = sig_count / expected,
        pvalue = NA_real_
      )
    
    for (i in 1:nrow(category_results)) {
      contingency <- matrix(c(
        category_results$sig_count[i],
        category_results$sig_total - category_results$sig_count[i],
        category_results$total_count[i] - category_results$sig_count[i],
        category_results$all_total - category_results$total_count[i] - 
          (category_results$sig_total - category_results$sig_count[i])
      ), nrow = 2)
      
      category_results$pvalue[i] <- tryCatch(fisher.test(contingency, alternative = "greater")$p.value, error = function(e) fisher.test(contingency, alternative = "greater", simulate.p.value = TRUE, B = 2000)$p.value)
    }
    
    category_results <- category_results %>%
      mutate(padj = p.adjust(pvalue, method = "BH")) %>%
      arrange(pvalue)
    
    write.csv(category_results,
              file.path(output_dir, sprintf("%s_Lipidmaps_Category_enrichment.csv", tolower(tissue_name))),
              row.names = FALSE)
    
    cat(sprintf("    Category: %d categories tested, %d significant (p<0.05)\n",
                nrow(category_results), sum(category_results$pvalue < 0.05, na.rm = TRUE)))
  }
  
  # Main Class enrichment
  sig_mainclass <- sig_data %>%
    filter(!is.na(MAIN_CLASS.Lipidmaps.) & MAIN_CLASS.Lipidmaps. != "" & MAIN_CLASS.Lipidmaps. != "-") %>%
    count(MAIN_CLASS.Lipidmaps., name = "sig_count")
  
  all_mainclass <- all_data %>%
    filter(!is.na(MAIN_CLASS.Lipidmaps.) & MAIN_CLASS.Lipidmaps. != "" & MAIN_CLASS.Lipidmaps. != "-") %>%
    count(MAIN_CLASS.Lipidmaps., name = "total_count")
  
  if (nrow(sig_mainclass) > 0 && nrow(all_mainclass) > 0) {
    mainclass_results <- sig_mainclass %>%
      left_join(all_mainclass, by = "MAIN_CLASS.Lipidmaps.") %>%
      mutate(
        sig_total = sum(sig_count),
        all_total = sum(total_count),
        expected = (total_count / all_total) * sig_total,
        fold_enrichment = sig_count / expected,
        pvalue = NA_real_
      )
    
    for (i in 1:nrow(mainclass_results)) {
      contingency <- matrix(c(
        mainclass_results$sig_count[i],
        mainclass_results$sig_total - mainclass_results$sig_count[i],
        mainclass_results$total_count[i] - mainclass_results$sig_count[i],
        mainclass_results$all_total - mainclass_results$total_count[i] - 
          (mainclass_results$sig_total - mainclass_results$sig_count[i])
      ), nrow = 2)
      
      mainclass_results$pvalue[i] <- tryCatch(fisher.test(contingency, alternative = "greater")$p.value, error = function(e) fisher.test(contingency, alternative = "greater", simulate.p.value = TRUE, B = 2000)$p.value)
    }
    
    mainclass_results <- mainclass_results %>%
      mutate(padj = p.adjust(pvalue, method = "BH")) %>%
      arrange(pvalue)
    
    write.csv(mainclass_results,
              file.path(output_dir, sprintf("%s_Lipidmaps_MainClass_enrichment.csv", tolower(tissue_name))),
              row.names = FALSE)
    
    cat(sprintf("    Main Class: %d classes tested, %d significant (p<0.05)\n",
                nrow(mainclass_results), sum(mainclass_results$pvalue < 0.05, na.rm = TRUE)))
  }
  
  # Sub Class enrichment
  sig_subclass <- sig_data %>%
    filter(!is.na(SUB_CLASS.Lipidmaps.) & SUB_CLASS.Lipidmaps. != "" & SUB_CLASS.Lipidmaps. != "-") %>%
    count(SUB_CLASS.Lipidmaps., name = "sig_count")
  
  all_subclass <- all_data %>%
    filter(!is.na(SUB_CLASS.Lipidmaps.) & SUB_CLASS.Lipidmaps. != "" & SUB_CLASS.Lipidmaps. != "-") %>%
    count(SUB_CLASS.Lipidmaps., name = "total_count")
  
  if (nrow(sig_subclass) > 0 && nrow(all_subclass) > 0) {
    subclass_results <- sig_subclass %>%
      left_join(all_subclass, by = "SUB_CLASS.Lipidmaps.") %>%
      mutate(
        sig_total = sum(sig_count),
        all_total = sum(total_count),
        expected = (total_count / all_total) * sig_total,
        fold_enrichment = sig_count / expected,
        pvalue = NA_real_
      )
    
    for (i in 1:nrow(subclass_results)) {
      contingency <- matrix(c(
        subclass_results$sig_count[i],
        subclass_results$sig_total - subclass_results$sig_count[i],
        subclass_results$total_count[i] - subclass_results$sig_count[i],
        subclass_results$all_total - subclass_results$total_count[i] - 
          (subclass_results$sig_total - subclass_results$sig_count[i])
      ), nrow = 2)
      
      subclass_results$pvalue[i] <- tryCatch(fisher.test(contingency, alternative = "greater")$p.value, error = function(e) fisher.test(contingency, alternative = "greater", simulate.p.value = TRUE, B = 2000)$p.value)
    }
    
    subclass_results <- subclass_results %>%
      mutate(padj = p.adjust(pvalue, method = "BH")) %>%
      arrange(pvalue)
    
    write.csv(subclass_results,
              file.path(output_dir, sprintf("%s_Lipidmaps_SubClass_enrichment.csv", tolower(tissue_name))),
              row.names = FALSE)
    
    cat(sprintf("    Sub Class: %d classes tested, %d significant (p<0.05)\n",
                nrow(subclass_results), sum(subclass_results$pvalue < 0.05, na.rm = TRUE)))
  }
  
  cat("\n")
}

perform_lipidmaps_enrichment(liver_sig, liver_all, "Liver")
perform_lipidmaps_enrichment(blood_sig, blood_all, "Blood")

# ----------------------------------------------------------------------------
# 2.3: SMPDB Pathway Enrichment (if HMDB IDs available)
# ----------------------------------------------------------------------------

cat("2.3: SMPDB Pathway Enrichment\n")
cat("  Note: SMPDB requires HMDB IDs for metabolite matching\n")

# Function to query SMPDB pathways for a given HMDB ID
query_smpdb_for_metabolite <- function(hmdb_id) {
  if (is.na(hmdb_id) || hmdb_id == "" || hmdb_id == "-") return(NULL)
  
  # SMPDB API endpoint (if available)
  # Note: SMPDB may require direct database download
  # This is a placeholder for the logic
  
  return(NULL)
}

cat("  SMPDB enrichment requires database download or API access\n")
cat("  Skipping automated SMPDB enrichment - manual analysis recommended\n\n")

# Save HMDB IDs for manual SMPDB analysis
liver_hmdb <- liver_sig %>%
  filter(!is.na(HMDB_ID) & HMDB_ID != "" & HMDB_ID != "-") %>%
  select(Compound_ID, Name, HMDB_ID, log2FC, pvalue, padj)

blood_hmdb <- blood_sig %>%
  filter(!is.na(HMDB_ID) & HMDB_ID != "" & HMDB_ID != "-") %>%
  select(Compound_ID, Name, HMDB_ID, log2FC, pvalue, padj)

write.csv(liver_hmdb,
          file.path(output_dir, "liver_HMDB_IDs_for_SMPDB.csv"),
          row.names = FALSE)
write.csv(blood_hmdb,
          file.path(output_dir, "blood_HMDB_IDs_for_SMPDB.csv"),
          row.names = FALSE)

cat(sprintf("  Saved %d liver and %d blood HMDB IDs for manual SMPDB analysis\n\n",
            nrow(liver_hmdb), nrow(blood_hmdb)))

# ----------------------------------------------------------------------------
# 2.4: Reactome Pathway Enrichment (using PubChem IDs)
# ----------------------------------------------------------------------------

cat("2.4: Reactome Pathway Enrichment\n")
cat("  Note: Reactome primarily focuses on protein pathways\n")
cat("  Limited metabolite coverage - manual analysis recommended\n\n")

# Save PubChem IDs for Reactome analysis
liver_pubchem <- liver_sig %>%
  filter(!is.na(PubChemID) & PubChemID != "" & PubChemID != "-") %>%
  select(Compound_ID, Name, PubChemID, log2FC, pvalue, padj)

blood_pubchem <- blood_sig %>%
  filter(!is.na(PubChemID) & PubChemID != "" & PubChemID != "-") %>%
  select(Compound_ID, Name, PubChemID, log2FC, pvalue, padj)

write.csv(liver_pubchem,
          file.path(output_dir, "liver_PubChem_IDs_for_Reactome.csv"),
          row.names = FALSE)
write.csv(blood_pubchem,
          file.path(output_dir, "blood_PubChem_IDs_for_Reactome.csv"),
          row.names = FALSE)

cat(sprintf("  Saved %d liver and %d blood PubChem IDs for manual Reactome analysis\n\n",
            nrow(liver_pubchem), nrow(blood_pubchem)))

cat("=== Enhanced enrichment analysis complete! ===\n")
cat(sprintf("Results saved to: %s\n", output_dir))
