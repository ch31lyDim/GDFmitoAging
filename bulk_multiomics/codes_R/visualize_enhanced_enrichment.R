#!/usr/bin/env Rscript
# Visualization of Enhanced Enrichment Results
# Creates comprehensive figures for all enrichment analyses

library(dplyr)
library(ggplot2)
library(scales)

# Set up directories
base_dir <- "/data1/dyh/GDF11/metabo_results_GDF11KO"
enrich_dir <- file.path(base_dir, "enhanced_enrichment")
fig_dir <- file.path(enrich_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Visualizing Enhanced Enrichment Results ===\n\n")

# ============================================================================
# Function to create enrichment bar plot
# ============================================================================

plot_enrichment_barplot <- function(data, title, max_terms = 15, 
                                   x_col = "sig_count", 
                                   term_col = NULL,
                                   pval_col = "pvalue") {
  
  # Auto-detect term column if not specified
  if (is.null(term_col)) {
    possible_cols <- setdiff(names(data), 
                            c("sig_count", "total_count", "sig_total", 
                              "all_total", "expected", "fold_enrichment",
                              "pvalue", "padj"))
    term_col <- possible_cols[1]
  }
  
  # Select top terms
  plot_data <- data %>%
    arrange(!!sym(pval_col)) %>%
    head(max_terms) %>%
    mutate(
      term = factor(!!sym(term_col), levels = rev(!!sym(term_col))),
      neglog10p = -log10(!!sym(pval_col)),
      significant = !!sym(pval_col) < 0.05
    )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = term, y = neglog10p, fill = significant)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "#E41A1C"),
                     labels = c("FALSE" = "p >= 0.05", "TRUE" = "p < 0.05")) +
    labs(
      title = title,
      x = "",
      y = expression(-log[10](p-value)),
      fill = "Significance"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# ============================================================================
# Function to create fold enrichment dot plot
# ============================================================================

plot_enrichment_dotplot <- function(data, title, max_terms = 15,
                                   term_col = NULL) {
  
  if (is.null(term_col)) {
    possible_cols <- setdiff(names(data), 
                            c("sig_count", "total_count", "sig_total", 
                              "all_total", "expected", "fold_enrichment",
                              "pvalue", "padj"))
    term_col <- possible_cols[1]
  }
  
  plot_data <- data %>%
    arrange(pvalue) %>%
    head(max_terms) %>%
    mutate(
      term = factor(!!sym(term_col), levels = rev(!!sym(term_col))),
      neglog10p = -log10(pvalue),
      significant = pvalue < 0.05
    )
  
  p <- ggplot(plot_data, aes(x = fold_enrichment, y = term)) +
    geom_point(aes(size = sig_count, color = neglog10p)) +
    scale_color_gradient(low = "blue", high = "red",
                        name = expression(-log[10](p))) +
    scale_size_continuous(name = "Count", range = c(3, 10)) +
    labs(
      title = title,
      x = "Fold Enrichment",
      y = ""
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# ============================================================================
# HMDB Enrichment Visualization
# ============================================================================

cat("Visualizing HMDB enrichment results...\n")

# Liver HMDB
liver_hmdb_files <- list(
  SuperClass = "liver_HMDB_SuperClass_enrichment.csv",
  Class = "liver_HMDB_Class_enrichment.csv",
  SubClass = "liver_HMDB_SubClass_enrichment.csv"
)

for (level in names(liver_hmdb_files)) {
  file_path <- file.path(enrich_dir, liver_hmdb_files[[level]])
  if (file.exists(file_path)) {
    data <- read.csv(file_path)
    if (nrow(data) > 0) {
      # Bar plot
      p1 <- plot_enrichment_barplot(data, 
                                    sprintf("Liver HMDB %s Enrichment", level),
                                    max_terms = 20)
      ggsave(file.path(fig_dir, sprintf("liver_HMDB_%s_barplot.png", level)),
             p1, width = 10, height = 8, dpi = 300)
      
      # Dot plot
      p2 <- plot_enrichment_dotplot(data,
                                    sprintf("Liver HMDB %s Enrichment", level),
                                    max_terms = 20)
      ggsave(file.path(fig_dir, sprintf("liver_HMDB_%s_dotplot.png", level)),
             p2, width = 10, height = 8, dpi = 300)
      
      cat(sprintf("  Created plots for Liver HMDB %s\n", level))
    }
  }
}

# Blood HMDB
blood_hmdb_files <- list(
  SuperClass = "blood_HMDB_SuperClass_enrichment.csv",
  Class = "blood_HMDB_Class_enrichment.csv",
  SubClass = "blood_HMDB_SubClass_enrichment.csv"
)

for (level in names(blood_hmdb_files)) {
  file_path <- file.path(enrich_dir, blood_hmdb_files[[level]])
  if (file.exists(file_path)) {
    data <- read.csv(file_path)
    if (nrow(data) > 0) {
      p1 <- plot_enrichment_barplot(data, 
                                    sprintf("Blood HMDB %s Enrichment", level),
                                    max_terms = 20)
      ggsave(file.path(fig_dir, sprintf("blood_HMDB_%s_barplot.png", level)),
             p1, width = 10, height = 8, dpi = 300)
      
      p2 <- plot_enrichment_dotplot(data,
                                    sprintf("Blood HMDB %s Enrichment", level),
                                    max_terms = 20)
      ggsave(file.path(fig_dir, sprintf("blood_HMDB_%s_dotplot.png", level)),
             p2, width = 10, height = 8, dpi = 300)
      
      cat(sprintf("  Created plots for Blood HMDB %s\n", level))
    }
  }
}

# ============================================================================
# Lipidmaps Enrichment Visualization
# ============================================================================

cat("\nVisualizing Lipidmaps enrichment results...\n")

# Liver Lipidmaps
liver_lipid_files <- list(
  Category = "liver_Lipidmaps_Category_enrichment.csv",
  MainClass = "liver_Lipidmaps_MainClass_enrichment.csv",
  SubClass = "liver_Lipidmaps_SubClass_enrichment.csv"
)

for (level in names(liver_lipid_files)) {
  file_path <- file.path(enrich_dir, liver_lipid_files[[level]])
  if (file.exists(file_path)) {
    data <- read.csv(file_path)
    if (nrow(data) > 0) {
      p1 <- plot_enrichment_barplot(data, 
                                    sprintf("Liver Lipidmaps %s Enrichment", level),
                                    max_terms = 20)
      ggsave(file.path(fig_dir, sprintf("liver_Lipidmaps_%s_barplot.png", level)),
             p1, width = 10, height = 8, dpi = 300)
      
      p2 <- plot_enrichment_dotplot(data,
                                    sprintf("Liver Lipidmaps %s Enrichment", level),
                                    max_terms = 20)
      ggsave(file.path(fig_dir, sprintf("liver_Lipidmaps_%s_dotplot.png", level)),
             p2, width = 10, height = 8, dpi = 300)
      
      cat(sprintf("  Created plots for Liver Lipidmaps %s\n", level))
    }
  }
}

# Blood Lipidmaps
blood_lipid_files <- list(
  Category = "blood_Lipidmaps_Category_enrichment.csv",
  MainClass = "blood_Lipidmaps_MainClass_enrichment.csv",
  SubClass = "blood_Lipidmaps_SubClass_enrichment.csv"
)

for (level in names(blood_lipid_files)) {
  file_path <- file.path(enrich_dir, blood_lipid_files[[level]])
  if (file.exists(file_path)) {
    data <- read.csv(file_path)
    if (nrow(data) > 0) {
      p1 <- plot_enrichment_barplot(data, 
                                    sprintf("Blood Lipidmaps %s Enrichment", level),
                                    max_terms = 20)
      ggsave(file.path(fig_dir, sprintf("blood_Lipidmaps_%s_barplot.png", level)),
             p1, width = 10, height = 8, dpi = 300)
      
      p2 <- plot_enrichment_dotplot(data,
                                    sprintf("Blood Lipidmaps %s Enrichment", level),
                                    max_terms = 20)
      ggsave(file.path(fig_dir, sprintf("blood_Lipidmaps_%s_dotplot.png", level)),
             p2, width = 10, height = 8, dpi = 300)
      
      cat(sprintf("  Created plots for Blood Lipidmaps %s\n", level))
    }
  }
}

# ============================================================================
# Create Combined Summary Figure
# ============================================================================

cat("\nCreating combined summary figure...\n")

create_summary_comparison <- function() {
  # Collect top results from each analysis
  summaries <- list()
  
  # Liver HMDB SuperClass
  file_path <- file.path(enrich_dir, "liver_HMDB_SuperClass_enrichment.csv")
  if (file.exists(file_path)) {
    data <- read.csv(file_path) %>%
      arrange(pvalue) %>%
      head(5) %>%
      mutate(
        Tissue = "Liver",
        Database = "HMDB SuperClass",
        Term = SuperClass.HMDB.
      ) %>%
      select(Tissue, Database, Term, sig_count, fold_enrichment, pvalue)
    summaries[[length(summaries) + 1]] <- data
  }
  
  # Blood HMDB SuperClass
  file_path <- file.path(enrich_dir, "blood_HMDB_SuperClass_enrichment.csv")
  if (file.exists(file_path)) {
    data <- read.csv(file_path) %>%
      arrange(pvalue) %>%
      head(5) %>%
      mutate(
        Tissue = "Blood",
        Database = "HMDB SuperClass",
        Term = SuperClass.HMDB.
      ) %>%
      select(Tissue, Database, Term, sig_count, fold_enrichment, pvalue)
    summaries[[length(summaries) + 1]] <- data
  }
  
  # Liver Lipidmaps Category
  file_path <- file.path(enrich_dir, "liver_Lipidmaps_Category_enrichment.csv")
  if (file.exists(file_path)) {
    data <- read.csv(file_path) %>%
      arrange(pvalue) %>%
      head(5) %>%
      mutate(
        Tissue = "Liver",
        Database = "Lipidmaps Category",
        Term = CATEGORY.Lipidmaps.
      ) %>%
      select(Tissue, Database, Term, sig_count, fold_enrichment, pvalue)
    summaries[[length(summaries) + 1]] <- data
  }
  
  # Blood Lipidmaps Category
  file_path <- file.path(enrich_dir, "blood_Lipidmaps_Category_enrichment.csv")
  if (file.exists(file_path)) {
    data <- read.csv(file_path) %>%
      arrange(pvalue) %>%
      head(5) %>%
      mutate(
        Tissue = "Blood",
        Database = "Lipidmaps Category",
        Term = CATEGORY.Lipidmaps.
      ) %>%
      select(Tissue, Database, Term, sig_count, fold_enrichment, pvalue)
    summaries[[length(summaries) + 1]] <- data
  }
  
  if (length(summaries) > 0) {
    combined <- bind_rows(summaries) %>%
      mutate(
        neglog10p = -log10(pvalue),
        significant = pvalue < 0.05,
        label = sprintf("%s\n%s", Database, Term)
      )
    
    # Create faceted plot
    p <- ggplot(combined, aes(x = fold_enrichment, y = reorder(Term, -pvalue))) +
      geom_point(aes(size = sig_count, color = neglog10p)) +
      facet_grid(Tissue ~ Database, scales = "free_y", space = "free_y") +
      scale_color_gradient(low = "blue", high = "red",
                          name = expression(-log[10](p))) +
      scale_size_continuous(name = "Count", range = c(3, 8)) +
      labs(
        title = "Top Enriched Terms Across Databases",
        x = "Fold Enrichment",
        y = ""
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9)
      )
    
    ggsave(file.path(fig_dir, "combined_enrichment_summary.png"),
           p, width = 14, height = 10, dpi = 300)
    
    cat("  Created combined summary figure\n")
  }
}

create_summary_comparison()

# ============================================================================
# Generate Summary Report
# ============================================================================

cat("\nGenerating enrichment summary report...\n")

report_lines <- c(
  "# Enhanced Pathway Enrichment Analysis Report",
  "",
  sprintf("**Analysis Date:** %s", Sys.Date()),
  "",
  "## Overview",
  "",
  "This report summarizes the results of enhanced pathway enrichment analysis using multiple databases:",
  "1. **KEGG**: Complete pathway annotations from KEGG REST API",
  "2. **HMDB**: Chemical class enrichment (SuperClass, Class, SubClass)",
  "3. **Lipidmaps**: Lipid category enrichment (Category, Main Class, Sub Class)",
  "4. **SMPDB/Reactome**: IDs prepared for manual analysis",
  "",
  "## Methods",
  "",
  "- **Statistical Test**: Fisher's exact test (one-tailed, greater)",
  "- **Multiple Testing Correction**: Benjamini-Hochberg FDR",
  "- **Significance Threshold**: p < 0.05",
  "",
  "## Results Summary",
  ""
)

# Add HMDB results summary
for (tissue in c("liver", "blood")) {
  report_lines <- c(report_lines,
                   sprintf("### %s - HMDB Enrichment", tools::toTitleCase(tissue)),
                   "")
  
  for (level in c("SuperClass", "Class", "SubClass")) {
    file_path <- file.path(enrich_dir, 
                          sprintf("%s_HMDB_%s_enrichment.csv", tissue, level))
    if (file.exists(file_path)) {
      data <- read.csv(file_path)
      sig_count <- sum(data$pvalue < 0.05, na.rm = TRUE)
      fdr_count <- sum(data$padj < 0.05, na.rm = TRUE)
      
      report_lines <- c(report_lines,
                       sprintf("**%s:**", level),
                       sprintf("- Total terms tested: %d", nrow(data)),
                       sprintf("- Significant (p<0.05): %d", sig_count),
                       sprintf("- Significant (FDR<0.05): %d", fdr_count),
                       "")
      
      if (sig_count > 0) {
        top_terms <- data %>%
          filter(pvalue < 0.05) %>%
          arrange(pvalue) %>%
          head(5)
        
        report_lines <- c(report_lines,
                         "Top enriched terms:",
                         "")
        
        for (i in 1:nrow(top_terms)) {
          term_name <- top_terms[[1]][i]  # First column is the term
          report_lines <- c(report_lines,
                           sprintf("%d. **%s** (Count: %d, Fold: %.2f, p=%.2e)",
                                  i, term_name, top_terms$sig_count[i],
                                  top_terms$fold_enrichment[i], top_terms$pvalue[i]))
        }
        report_lines <- c(report_lines, "")
      }
    }
  }
}

# Add Lipidmaps results summary
for (tissue in c("liver", "blood")) {
  report_lines <- c(report_lines,
                   sprintf("### %s - Lipidmaps Enrichment", tools::toTitleCase(tissue)),
                   "")
  
  for (level in c("Category", "MainClass", "SubClass")) {
    file_path <- file.path(enrich_dir,
                          sprintf("%s_Lipidmaps_%s_enrichment.csv", tissue, level))
    if (file.exists(file_path)) {
      data <- read.csv(file_path)
      sig_count <- sum(data$pvalue < 0.05, na.rm = TRUE)
      fdr_count <- sum(data$padj < 0.05, na.rm = TRUE)
      
      report_lines <- c(report_lines,
                       sprintf("**%s:**", level),
                       sprintf("- Total terms tested: %d", nrow(data)),
                       sprintf("- Significant (p<0.05): %d", sig_count),
                       sprintf("- Significant (FDR<0.05): %d", fdr_count),
                       "")
      
      if (sig_count > 0) {
        top_terms <- data %>%
          filter(pvalue < 0.05) %>%
          arrange(pvalue) %>%
          head(5)
        
        report_lines <- c(report_lines,
                         "Top enriched terms:",
                         "")
        
        for (i in 1:nrow(top_terms)) {
          term_name <- top_terms[[1]][i]
          report_lines <- c(report_lines,
                           sprintf("%d. **%s** (Count: %d, Fold: %.2f, p=%.2e)",
                                  i, term_name, top_terms$sig_count[i],
                                  top_terms$fold_enrichment[i], top_terms$pvalue[i]))
        }
        report_lines <- c(report_lines, "")
      }
    }
  }
}

# Add files generated section
report_lines <- c(report_lines,
                 "## Files Generated",
                 "",
                 "### Enrichment Results (CSV)",
                 "")

all_files <- list.files(enrich_dir, pattern = "*.csv", full.names = FALSE)
for (f in all_files) {
  report_lines <- c(report_lines, sprintf("- `%s`", f))
}

report_lines <- c(report_lines,
                 "",
                 "### Figures (PNG)",
                 "")

all_figures <- list.files(fig_dir, pattern = "*.png", full.names = FALSE)
for (f in all_figures) {
  report_lines <- c(report_lines, sprintf("- `%s`", f))
}

report_lines <- c(report_lines,
                 "",
                 "## Next Steps",
                 "",
                 "1. Review top enriched terms for biological interpretation",
                 "2. Use HMDB IDs for SMPDB pathway analysis (web interface)",
                 "3. Use PubChem IDs for Reactome analysis (web interface)",
                 "4. Compare enrichment patterns between KEGG, HMDB, and Lipidmaps",
                 "5. Investigate metabolite-level details for key enriched categories",
                 "")

writeLines(report_lines, file.path(enrich_dir, "ENHANCED_ENRICHMENT_REPORT.md"))

cat(sprintf("\n=== Visualization complete! ===\n"))
cat(sprintf("Figures saved to: %s\n", fig_dir))
cat(sprintf("Report saved to: %s\n", 
            file.path(enrich_dir, "ENHANCED_ENRICHMENT_REPORT.md")))
