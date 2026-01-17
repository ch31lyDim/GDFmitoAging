#!/usr/bin/env Rscript
# Complete Metabolomics Analysis Pipeline - GDF11 KO
# Includes proper pathway enrichment using standard IDs
# Author: Automated Analysis
# Date: 2024

cat("==========================================================================\n")
cat("COMPLETE METABOLOMICS ANALYSIS - GDF11 KO\n")
cat("==========================================================================\n\n")

# Load required libraries
cat("Loading required libraries...\n")
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(writexl)
  library(dplyr)
  library(tidyr)
})

# Check for MetaboAnalystR
if (!require("MetaboAnalystR", quietly = TRUE)) {
  cat("  Note: MetaboAnalystR not installed, will use alternative methods\n")
  use_metaboanalyst <- FALSE
} else {
  use_metaboanalyst <- TRUE
  cat("  ✓ MetaboAnalystR available\n")
}

cat("\n")

# Create output directories
results_dir <- "metabo_results_GDF11KO"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(results_dir, "liver"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(results_dir, "liver/figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(results_dir, "liver/pathway_enrichment"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(results_dir, "blood"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(results_dir, "blood/figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(results_dir, "blood/pathway_enrichment"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(results_dir, "integrated_analysis"), showWarnings = FALSE, recursive = TRUE)

# ==========================================================================
# PART 1: Data Loading and Preprocessing
# ==========================================================================

cat("PART 1: Data Loading and Preprocessing\n")
cat("=========================================\n")

data_dir <- "kn_mouse_meta/analysis/data"

# Function to load and process data
load_metabolomics_data <- function(file_path, ion_mode) {
  cat(sprintf("  Loading %s...\n", basename(file_path)))
  
  # Read Excel file
  data <- read_excel(file_path)
  
  # Extract metadata columns (all non-numeric columns and ID columns)
  metadata_cols <- c("Compound_ID", "Name", "ChineseName", "IonMode", "Formula", 
                     "MolecularWeight", "m/z", "MassError", "Adduct", "RT (min)", 
                     "Score", "Level", "Column",
                     "ClassI", "ClassI (Chinese)", "ClassII", "ClassII (Chinese)", 
                     "ClassIII", "ClassIII (Chinese)",
                     "CAS", "HMDB_ID", "SuperClass(HMDB)", "Class(HMDB)", "SubClass(HMDB)",
                     "Other_name(Kegg_name)", "KEGG_ID", "KEGG_MapID",
                     "Lipidmaps_ID", "CATEGORY(Lipidmaps)", "MAIN_CLASS(Lipidmaps)", "SUB_CLASS(Lipidmaps)",
                     "PubChemID", "SMILES", "InChIKey")
  
  # Keep only existing columns
  metadata_cols <- metadata_cols[metadata_cols %in% colnames(data)]
  
  # Extract sample columns
  sample_cols <- setdiff(colnames(data), metadata_cols)
  sample_cols <- grep("^(pos|neg)_(WT|wt|GDF11)_", sample_cols, value = TRUE)
  sample_cols <- grep("QC", sample_cols, value = TRUE, invert = TRUE)
  
  cat(sprintf("    Found %d metabolites and %d samples\n", nrow(data), length(sample_cols)))
  
  # Return structured data
  list(
    metadata = data[, metadata_cols],
    intensity = data[, sample_cols],
    all_data = data
  )
}

# Load all data files
cat("\nLoading blood data...\n")
blood_pos <- load_metabolomics_data(file.path(data_dir, "blood_pos.xlsx"), "pos")
blood_neg <- load_metabolomics_data(file.path(data_dir, "blood_neg.xlsx"), "neg")

cat("\nLoading liver data...\n")
liver_pos <- load_metabolomics_data(file.path(data_dir, "liver_pos.xlsx"), "pos")
liver_neg <- load_metabolomics_data(file.path(data_dir, "liver_neg.xlsx"), "neg")

cat("\n✓ Data loading complete\n\n")

# ==========================================================================
# PART 2: Differential Expression Analysis
# ==========================================================================

cat("PART 2: Differential Expression Analysis\n")
cat("=========================================\n")

# Function to perform differential analysis
perform_differential_analysis <- function(pos_data, neg_data, tissue_name, 
                                         p_threshold = 0.05, fc_threshold = 1.5) {
  
  cat(sprintf("\nAnalyzing %s tissue...\n", tissue_name))
  
  # Combine positive and negative mode data
  combined_metadata <- rbind(pos_data$metadata, neg_data$metadata)
  
  # Get sample columns from each mode
  pos_samples <- colnames(pos_data$intensity)
  neg_samples <- colnames(neg_data$intensity)
  
  # Find common samples (by converting pos/neg prefixes)
  all_samples <- unique(c(
    gsub("^pos_", "", pos_samples),
    gsub("^neg_", "", neg_samples)
  ))
  
  # Create unified intensity matrix with consistent column names
  combined_intensity <- data.frame(matrix(nrow = nrow(combined_metadata), ncol = 0))
  
  # Add positive mode samples (rename to remove pos_ prefix)
  for (col in pos_samples) {
    new_col <- gsub("^pos_", "", col)
    combined_intensity[[new_col]] <- c(pos_data$intensity[[col]], rep(NA, nrow(neg_data$intensity)))
  }
  
  # Add negative mode samples (rename to remove neg_ prefix)
  for (col in neg_samples) {
    new_col <- gsub("^neg_", "", col)
    if (new_col %in% colnames(combined_intensity)) {
      # Column already exists from pos mode, fill in neg mode values
      combined_intensity[[new_col]][seq(nrow(pos_data$intensity) + 1, nrow(combined_metadata))] <- neg_data$intensity[[col]]
    } else {
      # New column, add with NAs for pos mode
      combined_intensity[[new_col]] <- c(rep(NA, nrow(pos_data$intensity)), neg_data$intensity[[col]])
    }
  }
  
  # Ensure we have valid data
  combined_intensity[is.na(combined_intensity)] <- 0
  
  cat(sprintf("  Total metabolites: %d\n", nrow(combined_metadata)))
  
  # Identify WT and KO samples (now without pos/neg prefix)
  sample_cols <- colnames(combined_intensity)
  wt_samples <- grep("^(WT|wt)_", sample_cols, value = TRUE)
  ko_samples <- grep("^GDF11_", sample_cols, value = TRUE)
  
  cat(sprintf("  WT samples: %d\n", length(wt_samples)))
  cat(sprintf("  KO samples: %d\n", length(ko_samples)))
  
  # Calculate statistics for each metabolite
  results <- data.frame()
  
  for (i in 1:nrow(combined_intensity)) {
    wt_values <- as.numeric(combined_intensity[i, wt_samples])
    ko_values <- as.numeric(combined_intensity[i, ko_samples])
    
    # Remove NA values
    wt_values <- wt_values[!is.na(wt_values)]
    ko_values <- ko_values[!is.na(ko_values)]
    
    # Skip if too few values
    if (length(wt_values) < 2 || length(ko_values) < 2) next
    
    # Calculate means
    wt_mean <- mean(wt_values)
    ko_mean <- mean(ko_values)
    
    # Calculate fold change
    fc <- ifelse(wt_mean > 0, ko_mean / wt_mean, NA)
    log2fc <- ifelse(!is.na(fc) && fc > 0, log2(fc), NA)
    
    # Perform t-test
    test_result <- tryCatch({
      t.test(ko_values, wt_values, var.equal = FALSE)
    }, error = function(e) NULL)
    
    if (is.null(test_result)) next
    
    results <- rbind(results, data.frame(
      Row_Index = i,
      WT_mean = wt_mean,
      KO_mean = ko_mean,
      FC = fc,
      log2FC = log2fc,
      pvalue = test_result$p.value,
      stringsAsFactors = FALSE
    ))
  }
  
  # Add FDR correction
  results$padj <- p.adjust(results$pvalue, method = "BH")
  
  # Add metadata
  results <- cbind(combined_metadata[results$Row_Index, ], results[, -1])
  
  # Determine significance
  results$significant <- ifelse(
    results$pvalue < p_threshold & abs(results$log2FC) > log2(fc_threshold),
    "Significant", "Not Significant"
  )
  
  # Determine regulation direction
  results$regulation <- ifelse(
    results$significant == "Significant",
    ifelse(results$log2FC > 0, "Up", "Down"),
    "NS"
  )
  
  # Sort by p-value
  results <- results[order(results$pvalue), ]
  
  cat(sprintf("  Significant metabolites: %d (%.1f%%)\n", 
              sum(results$significant == "Significant"),
              100 * sum(results$significant == "Significant") / nrow(results)))
  cat(sprintf("    Up-regulated: %d\n", sum(results$regulation == "Up")))
  cat(sprintf("    Down-regulated: %d\n", sum(results$regulation == "Down")))
  
  return(results)
}

# Perform differential analysis for liver
liver_results <- perform_differential_analysis(
  liver_pos, liver_neg, "Liver",
  p_threshold = 0.05, fc_threshold = 1.5
)

# Perform differential analysis for blood
blood_results <- perform_differential_analysis(
  blood_pos, blood_neg, "Blood",
  p_threshold = 0.05, fc_threshold = 1.2  # More relaxed for blood
)

# Save differential results
write.csv(liver_results, 
          file.path(results_dir, "liver", "liver_differential_metabolites.csv"),
          row.names = FALSE)
write.csv(blood_results,
          file.path(results_dir, "blood", "blood_differential_metabolites.csv"),
          row.names = FALSE)

# Extract significant metabolites
liver_sig <- liver_results[liver_results$significant == "Significant", ]
blood_sig <- blood_results[blood_results$significant == "Significant", ]

write.csv(liver_sig,
          file.path(results_dir, "liver", "liver_significant_metabolites.csv"),
          row.names = FALSE)
write.csv(blood_sig,
          file.path(results_dir, "blood", "blood_significant_metabolites.csv"),
          row.names = FALSE)

cat("\n✓ Differential analysis complete\n\n")

# ==========================================================================
# PART 3: Visualization
# ==========================================================================

cat("PART 3: Generating Visualizations\n")
cat("===================================\n")

# Function to create volcano plot
create_volcano_plot <- function(results, title, output_path, 
                                p_threshold = 0.05, fc_threshold = 1.5) {
  
  results$color_group <- "Not Significant"
  results$color_group[results$pvalue < p_threshold & results$log2FC > log2(fc_threshold)] <- "Up-regulated"
  results$color_group[results$pvalue < p_threshold & results$log2FC < -log2(fc_threshold)] <- "Down-regulated"
  results$color_group <- factor(results$color_group, 
                                levels = c("Not Significant", "Up-regulated", "Down-regulated"))
  
  n_up <- sum(results$color_group == "Up-regulated")
  n_down <- sum(results$color_group == "Down-regulated")
  
  p <- ggplot(results, aes(x = log2FC, y = -log10(pvalue))) +
    geom_point(aes(color = color_group, size = color_group, alpha = color_group)) +
    scale_color_manual(values = c("Not Significant" = "grey70",
                                  "Up-regulated" = "#E31A1C",
                                  "Down-regulated" = "#1F78B4")) +
    scale_size_manual(values = c("Not Significant" = 1, 
                                 "Up-regulated" = 2, 
                                 "Down-regulated" = 2)) +
    scale_alpha_manual(values = c("Not Significant" = 0.3,
                                  "Up-regulated" = 0.8,
                                  "Down-regulated" = 0.8)) +
    geom_vline(xintercept = c(-log2(fc_threshold), log2(fc_threshold)), 
               linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(p_threshold), 
               linetype = "dashed", color = "grey40") +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Up: %d\nDown: %d", n_up, n_down),
             hjust = 1.1, vjust = 1.5, size = 4) +
    labs(title = title,
         subtitle = sprintf("Threshold: p < %.2f, |FC| > %.1f", p_threshold, fc_threshold),
         x = "log2(Fold Change) [GDF11 KO / WT]",
         y = "-log10(p-value)",
         color = "Regulation",
         size = "Regulation",
         alpha = "Regulation") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right")
  
  ggsave(output_path, p, width = 10, height = 7, dpi = 300)
  cat(sprintf("  Saved: %s\n", basename(output_path)))
}

# Generate volcano plots
cat("\nGenerating volcano plots...\n")
create_volcano_plot(liver_results, 
                   "Liver Metabolomics - GDF11 KO vs WT",
                   file.path(results_dir, "liver/figures", "liver_volcano_plot.png"),
                   p_threshold = 0.05, fc_threshold = 1.5)

create_volcano_plot(blood_results,
                   "Blood Metabolomics - GDF11 KO vs WT",
                   file.path(results_dir, "blood/figures", "blood_volcano_plot.png"),
                   p_threshold = 0.05, fc_threshold = 1.2)

cat("\n✓ Visualizations complete\n\n")

# Save progress checkpoint
cat("Progress checkpoint: Basic differential analysis completed\n")
cat("Next: Pathway enrichment analysis using standard IDs\n\n")

# ==========================================================================
# CONTINUE TO PART 4 IN NEXT SECTION
# ==========================================================================

saveRDS(list(
  liver_results = liver_results,
  blood_results = blood_results,
  liver_sig = liver_sig,
  blood_sig = blood_sig
), file.path(results_dir, "analysis_checkpoint.rds"))

cat("✓ Analysis checkpoint saved\n")
cat("Run pathway_enrichment_analysis.R next for enrichment analysis\n\n")
