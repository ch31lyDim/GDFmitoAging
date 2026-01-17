#!/usr/bin/env Rscript
# Pathway Enrichment Analysis using Standard IDs
# Includes ORA and GSEA with multiple metabolite sets
# Author: Automated Analysis
# Date: 2024

cat("==========================================================================\n")
cat("PATHWAY ENRICHMENT ANALYSIS\n")
cat("==========================================================================\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(writexl)
})

# Load checkpoint data
results_dir <- "metabo_results_GDF11KO"
checkpoint <- readRDS(file.path(results_dir, "analysis_checkpoint.rds"))
liver_results <- checkpoint$liver_results
blood_results <- checkpoint$blood_results
liver_sig <- checkpoint$liver_sig
blood_sig <- checkpoint$blood_sig

cat("Loaded checkpoint data:\n")
cat(sprintf("  Liver: %d total, %d significant\n", nrow(liver_results), nrow(liver_sig)))
cat(sprintf("  Blood: %d total, %d significant\n", nrow(blood_results), nrow(blood_sig)))
cat("\n")

# ==========================================================================
# PART 1: Prepare ID Mapping
# ==========================================================================

cat("PART 1: Preparing ID Mapping\n")
cat("==============================\n")

# Function to clean and prepare IDs
prepare_id_mapping <- function(results, id_col) {
  # Extract non-missing IDs
  valid_ids <- results[[id_col]][!is.na(results[[id_col]]) & 
                                   results[[id_col]] != "-" & 
                                   results[[id_col]] != ""]
  return(valid_ids)
}

# Check ID coverage
check_id_coverage <- function(results, id_columns) {
  cat("\nID Coverage Analysis:\n")
  for (id_col in id_columns) {
    if (id_col %in% colnames(results)) {
      valid <- sum(!is.na(results[[id_col]]) & 
                   results[[id_col]] != "-" & 
                   results[[id_col]] != "")
      pct <- 100 * valid / nrow(results)
      cat(sprintf("  %s: %d / %d (%.1f%%)\n", id_col, valid, nrow(results), pct))
    }
  }
}

cat("\nLiver ID Coverage:\n")
check_id_coverage(liver_results, c("HMDB_ID", "KEGG_ID", "PubChemID", "KEGG_MapID"))

cat("\nBlood ID Coverage:\n")
check_id_coverage(blood_results, c("HMDB_ID", "KEGG_ID", "PubChemID", "KEGG_MapID"))

# ==========================================================================
# PART 2: Manual ORA using KEGG Pathways
# ==========================================================================

cat("\n\nPART 2: Over-Representation Analysis (ORA)\n")
cat("============================================\n")

# Load KEGG pathway database
kegg_db <- read.csv("kn_mouse_meta/analysis/scripts/kegg_pathway.csv")
cat(sprintf("\nLoaded KEGG pathway database: %d pathways\n", nrow(kegg_db)))

# Parse pathway database to create ID-to-pathway mapping
parse_kegg_for_ids <- function(kegg_df) {
  pathway_map <- list()
  
  for (i in 1:nrow(kegg_df)) {
    pathway_id <- kegg_df$id[i]
    pathway_name <- kegg_df$name[i]
    
    # Note: This is a simplified approach
    # The actual KEGG database would need proper compound ID mapping
    pathway_map[[pathway_id]] <- list(
      name = pathway_name,
      compounds = character(0)  # Would need proper KEGG compound IDs
    )
  }
  
  return(pathway_map)
}

# Function to perform ORA with compound names
perform_name_based_ora <- function(sig_compounds, all_compounds, kegg_df, min_hits = 3) {
  results <- data.frame()
  
  for (i in 1:nrow(kegg_df)) {
    pathway_id <- kegg_df$id[i]
    pathway_name <- kegg_df$name[i]
    pathway_members <- strsplit(kegg_df$member[i], "; ")[[1]]
    
    # Match by compound names
    sig_in_pathway <- sum(sig_compounds %in% pathway_members)
    all_in_pathway <- sum(all_compounds %in% pathway_members)
    
    # Skip if too few hits
    if (sig_in_pathway < min_hits || all_in_pathway < min_hits) next
    
    # Hypergeometric test (Fisher's exact test)
    contingency <- matrix(c(
      sig_in_pathway,
      length(sig_compounds) - sig_in_pathway,
      all_in_pathway - sig_in_pathway,
      length(all_compounds) - length(sig_compounds) - (all_in_pathway - sig_in_pathway)
    ), nrow = 2)
    
    test_result <- fisher.test(contingency, alternative = "greater")
    
    results <- rbind(results, data.frame(
      Pathway_ID = pathway_id,
      Pathway_Name = pathway_name,
      Hits_in_Pathway = sig_in_pathway,
      Pathway_Size = all_in_pathway,
      Total_Significant = length(sig_compounds),
      Total_Measured = length(all_compounds),
      Pvalue = test_result$p.value,
      Odds_Ratio = as.numeric(test_result$estimate),
      stringsAsFactors = FALSE
    ))
  }
  
  if (nrow(results) > 0) {
    results$FDR <- p.adjust(results$Pvalue, method = "BH")
    results <- results[order(results$Pvalue), ]
    results$Enrichment_Ratio <- (results$Hits_in_Pathway / results$Total_Significant) / 
                                 (results$Pathway_Size / results$Total_Measured)
  }
  
  return(results)
}

# Perform ORA for liver
cat("\nPerforming ORA for Liver significant metabolites...\n")
liver_ora <- perform_name_based_ora(
  sig_compounds = liver_sig$Name,
  all_compounds = liver_results$Name,
  kegg_df = kegg_db,
  min_hits = 3
)

if (nrow(liver_ora) > 0) {
  write.csv(liver_ora,
            file.path(results_dir, "liver/pathway_enrichment", "liver_ORA_KEGG.csv"),
            row.names = FALSE)
  
  cat(sprintf("  Found %d pathways with hits\n", nrow(liver_ora)))
  cat(sprintf("  Significant pathways (p<0.05): %d\n", sum(liver_ora$Pvalue < 0.05)))
  cat(sprintf("  FDR-significant (FDR<0.05): %d\n", sum(liver_ora$FDR < 0.05)))
  
  if (nrow(liver_ora) > 0) {
    cat("\n  Top 10 enriched pathways:\n")
    print(head(liver_ora[, c("Pathway_Name", "Hits_in_Pathway", "Pvalue", "FDR")], 10))
  }
} else {
  cat("  No enriched pathways found\n")
}

# Perform ORA for blood
cat("\nPerforming ORA for Blood significant metabolites...\n")
blood_ora <- perform_name_based_ora(
  sig_compounds = blood_sig$Name,
  all_compounds = blood_results$Name,
  kegg_df = kegg_db,
  min_hits = 2  # Lower threshold for blood
)

if (nrow(blood_ora) > 0) {
  write.csv(blood_ora,
            file.path(results_dir, "blood/pathway_enrichment", "blood_ORA_KEGG.csv"),
            row.names = FALSE)
  
  cat(sprintf("  Found %d pathways with hits\n", nrow(blood_ora)))
  cat(sprintf("  Significant pathways (p<0.05): %d\n", sum(blood_ora$Pvalue < 0.05)))
  
  if (sum(blood_ora$Pvalue < 0.05) > 0) {
    cat("\n  Top enriched pathways:\n")
    print(head(blood_ora[blood_ora$Pvalue < 0.05, 
                         c("Pathway_Name", "Hits_in_Pathway", "Pvalue", "FDR")], 10))
  }
} else {
  cat("  No enriched pathways found\n")
}

# ==========================================================================
# PART 3: GSEA (Gene Set Enrichment Analysis)
# ==========================================================================

cat("\n\nPART 3: Gene Set Enrichment Analysis (GSEA)\n")
cat("=============================================\n")

# Function to perform simplified GSEA
perform_simple_gsea <- function(ranked_data, pathway_sets, min_set_size = 5) {
  
  # Create ranked list by log2FC
  ranked_list <- setNames(ranked_data$log2FC, ranked_data$Name)
  ranked_list <- ranked_list[!is.na(ranked_list)]
  ranked_list <- sort(ranked_list, decreasing = TRUE)
  
  cat(sprintf("  Ranked list: %d metabolites\n", length(ranked_list)))
  
  results <- data.frame()
  
  for (i in 1:nrow(pathway_sets)) {
    pathway_name <- pathway_sets$name[i]
    pathway_members <- strsplit(pathway_sets$member[i], "; ")[[1]]
    
    # Find pathway members in ranked list
    hits <- names(ranked_list) %in% pathway_members
    
    if (sum(hits) < min_set_size) next
    
    # Calculate enrichment score (simplified running sum)
    hit_positions <- which(hits)
    n_hits <- length(hit_positions)
    n_total <- length(ranked_list)
    
    # Calculate normalized enrichment score
    running_sum <- cumsum(ifelse(hits, 1/n_hits, -1/(n_total - n_hits)))
    es <- max(abs(running_sum))
    
    # Get direction
    max_pos <- which.max(abs(running_sum))
    direction <- ifelse(running_sum[max_pos] > 0, "Up", "Down")
    
    # Simple permutation test (limited permutations for speed)
    n_perm <- 100
    null_es <- numeric(n_perm)
    for (j in 1:n_perm) {
      perm_hits <- sample(hits)
      perm_sum <- cumsum(ifelse(perm_hits, 1/n_hits, -1/(n_total - n_hits)))
      null_es[j] <- max(abs(perm_sum))
    }
    
    pvalue <- sum(null_es >= es) / n_perm
    
    results <- rbind(results, data.frame(
      Pathway_Name = pathway_name,
      Set_Size = n_hits,
      ES = es,
      NES = es / mean(null_es),  # Normalized ES
      Pvalue = pvalue,
      Direction = direction,
      stringsAsFactors = FALSE
    ))
  }
  
  if (nrow(results) > 0) {
    results$FDR <- p.adjust(results$Pvalue, method = "BH")
    results <- results[order(results$Pvalue), ]
  }
  
  return(results)
}

# Perform GSEA for liver
cat("\nPerforming GSEA for Liver (all metabolites)...\n")
liver_gsea <- perform_simple_gsea(
  ranked_data = liver_results,
  pathway_sets = kegg_db,
  min_set_size = 5
)

if (nrow(liver_gsea) > 0) {
  write.csv(liver_gsea,
            file.path(results_dir, "liver/pathway_enrichment", "liver_GSEA_KEGG.csv"),
            row.names = FALSE)
  
  cat(sprintf("  Analyzed %d pathways\n", nrow(liver_gsea)))
  cat(sprintf("  Nominally significant (p<0.05): %d\n", sum(liver_gsea$Pvalue < 0.05)))
  
  if (sum(liver_gsea$Pvalue < 0.05) > 0) {
    cat("\n  Top enriched pathways:\n")
    print(head(liver_gsea[liver_gsea$Pvalue < 0.05, 
                          c("Pathway_Name", "Set_Size", "NES", "Pvalue", "Direction")], 10))
  }
} else {
  cat("  No pathways analyzed\n")
}

# Perform GSEA for blood
cat("\nPerforming GSEA for Blood (all metabolites)...\n")
blood_gsea <- perform_simple_gsea(
  ranked_data = blood_results,
  pathway_sets = kegg_db,
  min_set_size = 3  # Lower threshold for blood
)

if (nrow(blood_gsea) > 0) {
  write.csv(blood_gsea,
            file.path(results_dir, "blood/pathway_enrichment", "blood_GSEA_KEGG.csv"),
            row.names = FALSE)
  
  cat(sprintf("  Analyzed %d pathways\n", nrow(blood_gsea)))
  cat(sprintf("  Nominally significant (p<0.05): %d\n", sum(blood_gsea$Pvalue < 0.05)))
  
  if (sum(blood_gsea$Pvalue < 0.05) > 0) {
    cat("\n  Top enriched pathways:\n")
    print(head(blood_gsea[blood_gsea$Pvalue < 0.05,
                          c("Pathway_Name", "Set_Size", "NES", "Pvalue", "Direction")], 10))
  }
}

# ==========================================================================
# PART 4: Visualization of Enrichment Results
# ==========================================================================

cat("\n\nPART 4: Visualizing Enrichment Results\n")
cat("========================================\n")

# Function to create enrichment bar plot
create_enrichment_plot <- function(ora_results, title, output_path, top_n = 20) {
  
  if (nrow(ora_results) == 0) {
    cat(sprintf("  Skipping %s - no results\n", basename(output_path)))
    return(NULL)
  }
  
  plot_data <- head(ora_results, top_n)
  plot_data$neglog10p <- -log10(plot_data$Pvalue)
  plot_data$Pathway_Name <- factor(plot_data$Pathway_Name,
                                   levels = rev(plot_data$Pathway_Name))
  plot_data$Significant <- plot_data$FDR < 0.05
  
  p <- ggplot(plot_data, aes(x = Pathway_Name, y = neglog10p, fill = Significant)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "#E31A1C"),
                     labels = c("FALSE" = "p < 0.05", "TRUE" = "FDR < 0.05")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) +
    coord_flip() +
    labs(title = title,
         subtitle = sprintf("Top %d pathways by p-value", min(top_n, nrow(plot_data))),
         x = "KEGG Pathway",
         y = "-log10(p-value)",
         fill = "Significance") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 8))
  
  ggsave(output_path, p, width = 12, height = max(8, top_n * 0.3), dpi = 300)
  cat(sprintf("  Saved: %s\n", basename(output_path)))
}

# Generate plots
cat("\nGenerating enrichment plots...\n")

if (exists("liver_ora") && nrow(liver_ora) > 0) {
  create_enrichment_plot(liver_ora,
                        "Liver - KEGG Pathway Enrichment (ORA)",
                        file.path(results_dir, "liver/figures", "liver_pathway_enrichment.png"))
}

if (exists("blood_ora") && nrow(blood_ora) > 0) {
  create_enrichment_plot(blood_ora,
                        "Blood - KEGG Pathway Enrichment (ORA)",
                        file.path(results_dir, "blood/figures", "blood_pathway_enrichment.png"))
}

# ==========================================================================
# PART 5: Summary Report
# ==========================================================================

cat("\n\nPART 5: Generating Summary Report\n")
cat("===================================\n")

summary_data <- list(
  liver_ora = if(exists("liver_ora")) liver_ora else data.frame(),
  blood_ora = if(exists("blood_ora")) blood_ora else data.frame(),
  liver_gsea = if(exists("liver_gsea")) liver_gsea else data.frame(),
  blood_gsea = if(exists("blood_gsea")) blood_gsea else data.frame()
)

saveRDS(summary_data, file.path(results_dir, "pathway_enrichment_summary.rds"))

cat("\n✓ Pathway enrichment analysis complete\n")
cat(sprintf("\nResults saved to: %s/\n", results_dir))
cat("\nRun integrated_analysis.R next for cross-tissue comparison\n\n")
