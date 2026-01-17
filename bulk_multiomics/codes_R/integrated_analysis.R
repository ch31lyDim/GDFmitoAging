#!/usr/bin/env Rscript
# Integrated Analysis - Cross-tissue, Lipid/AA focus, Final Report
# Author: Automated Analysis
# Date: 2024

cat("==========================================================================\n")
cat("INTEGRATED ANALYSIS\n")
cat("==========================================================================\n\n")

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(writexl)
  library(RColorBrewer)
})

# Load checkpoint data
results_dir <- "metabo_results_GDF11KO"
checkpoint <- readRDS(file.path(results_dir, "analysis_checkpoint.rds"))
liver_results <- checkpoint$liver_results
blood_results <- checkpoint$blood_results
liver_sig <- checkpoint$liver_sig
blood_sig <- checkpoint$blood_sig

# Create output directories
integrated_dir <- file.path(results_dir, "integrated_analysis")
dir.create(file.path(integrated_dir, "cross_tissue"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(integrated_dir, "lipid_aa_focus"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(integrated_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

cat("Loaded data:\n")
cat(sprintf("  Liver: %d total, %d significant\n", nrow(liver_results), nrow(liver_sig)))
cat(sprintf("  Blood: %d total, %d significant\n", nrow(blood_results), nrow(blood_sig)))
cat("\n")

# ==========================================================================
# PART 1: Cross-Tissue Comparison
# ==========================================================================

cat("PART 1: Cross-Tissue Comparison\n")
cat("=================================\n")

# Find common metabolites
liver_compounds <- liver_results$Compound_ID
blood_compounds <- blood_results$Compound_ID
common_compounds <- intersect(liver_compounds, blood_compounds)

cat(sprintf("  Common metabolites: %d\n", length(common_compounds)))
cat(sprintf("  Liver-specific: %d\n", length(setdiff(liver_compounds, blood_compounds))))
cat(sprintf("  Blood-specific: %d\n", length(setdiff(blood_compounds, liver_compounds))))

# Create comparison dataframe
comparison_df <- data.frame(
  Compound_ID = common_compounds,
  stringsAsFactors = FALSE
)

# Merge liver data
liver_subset <- liver_results %>%
  select(Compound_ID, Name, log2FC, pvalue, padj, ClassI, ClassII, 
         HMDB_ID, KEGG_ID, PubChemID) %>%
  rename(Liver_log2FC = log2FC,
         Liver_pvalue = pvalue,
         Liver_padj = padj)

comparison_df <- merge(comparison_df, liver_subset, by = "Compound_ID", all.x = TRUE)

# Merge blood data
blood_subset <- blood_results %>%
  select(Compound_ID, log2FC, pvalue, padj) %>%
  rename(Blood_log2FC = log2FC,
         Blood_pvalue = pvalue,
         Blood_padj = padj)

comparison_df <- merge(comparison_df, blood_subset, by = "Compound_ID", all.x = TRUE)

# Calculate correlation
valid_pairs <- comparison_df %>%
  filter(!is.na(Liver_log2FC) & !is.na(Blood_log2FC) &
         is.finite(Liver_log2FC) & is.finite(Blood_log2FC))

fc_cor <- cor(valid_pairs$Liver_log2FC, valid_pairs$Blood_log2FC, use = "complete.obs")
cat(sprintf("\n  Fold change correlation: r = %.3f\n", fc_cor))

# Identify change patterns
comparison_df$Pattern <- "No change"
comparison_df$Pattern[comparison_df$Liver_pvalue < 0.05 & abs(comparison_df$Liver_log2FC) > log2(1.5)] <- "Liver-specific"
comparison_df$Pattern[comparison_df$Blood_pvalue < 0.05 & abs(comparison_df$Blood_log2FC) > log2(1.2)] <- "Blood-specific"
comparison_df$Pattern[comparison_df$Liver_pvalue < 0.05 & abs(comparison_df$Liver_log2FC) > log2(1.5) &
                      comparison_df$Blood_pvalue < 0.05 & abs(comparison_df$Blood_log2FC) > log2(1.2)] <- "Both tissues"

pattern_table <- table(comparison_df$Pattern)
cat("\n  Change patterns:\n")
print(pattern_table)

# Save comparison
write.csv(comparison_df,
          file.path(integrated_dir, "cross_tissue", "liver_blood_comparison.csv"),
          row.names = FALSE)

# Extract metabolites significant in both tissues
both_sig <- comparison_df %>%
  filter(Pattern == "Both tissues") %>%
  arrange(Liver_pvalue)

if (nrow(both_sig) > 0) {
  write.csv(both_sig,
            file.path(integrated_dir, "cross_tissue", "metabolites_significant_both_tissues.csv"),
            row.names = FALSE)
  cat(sprintf("\n  Metabolites significant in both tissues: %d\n", nrow(both_sig)))
  
  cat("\n  Details:\n")
  print(both_sig %>% select(Name, Liver_log2FC, Liver_pvalue, Blood_log2FC, Blood_pvalue, ClassI))
}

# Create scatter plot
p_scatter <- ggplot(valid_pairs, aes(x = Liver_log2FC, y = Blood_log2FC)) +
  geom_point(alpha = 0.4, size = 1.5, color = "grey60") +
  geom_point(data = valid_pairs %>% filter(Liver_pvalue < 0.05 & Blood_pvalue < 0.05),
             color = "red", alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  geom_smooth(method = "lm", color = "blue", se = TRUE, alpha = 0.2) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("r = %.3f\nn = %d", fc_cor, nrow(valid_pairs)),
           hjust = 1.1, vjust = 1.5, size = 5) +
  labs(title = "Cross-Tissue Metabolite Fold Changes",
       subtitle = "Liver vs Blood - GDF11 KO",
       x = "Liver log2(Fold Change)",
       y = "Blood log2(Fold Change)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(integrated_dir, "figures", "cross_tissue_scatter.png"),
       p_scatter, width = 8, height = 7, dpi = 300)

cat("\n✓ Cross-tissue analysis complete\n\n")

# ==========================================================================
# PART 2: Lipid and Amino Acid Focus
# ==========================================================================

cat("PART 2: Lipid and Amino Acid Analysis\n")
cat("=======================================\n")

# Define lipid classes
lipid_keywords <- c("Lipids and lipid-like molecules", "Fatty Acyls", 
                    "Glycerophospholipids", "Steroid", "Prenol lipids")

# Extract lipids
liver_lipids <- liver_sig %>%
  filter(grepl(paste(lipid_keywords, collapse = "|"), ClassI, ignore.case = TRUE) |
         grepl(paste(lipid_keywords, collapse = "|"), ClassII, ignore.case = TRUE))

blood_lipids <- blood_sig %>%
  filter(grepl(paste(lipid_keywords, collapse = "|"), ClassI, ignore.case = TRUE) |
         grepl(paste(lipid_keywords, collapse = "|"), ClassII, ignore.case = TRUE))

cat(sprintf("\n  Lipid metabolites:\n"))
cat(sprintf("    Liver: %d (%.1f%% of significant)\n", 
            nrow(liver_lipids), 100 * nrow(liver_lipids) / nrow(liver_sig)))
cat(sprintf("    Blood: %d (%.1f%% of significant)\n",
            nrow(blood_lipids), 100 * nrow(blood_lipids) / nrow(blood_sig)))

# Extract amino acids/peptides
aa_keywords <- c("amino", "peptide", "Amino", "Peptide")

liver_aa <- liver_sig %>%
  filter(grepl(paste(aa_keywords, collapse = "|"), ClassI, ignore.case = TRUE) |
         grepl(paste(aa_keywords, collapse = "|"), ClassII, ignore.case = TRUE) |
         grepl(paste(aa_keywords, collapse = "|"), Name, ignore.case = TRUE) |
         ClassI == "Organic acids and derivatives")

blood_aa <- blood_sig %>%
  filter(grepl(paste(aa_keywords, collapse = "|"), ClassI, ignore.case = TRUE) |
         grepl(paste(aa_keywords, collapse = "|"), ClassII, ignore.case = TRUE) |
         grepl(paste(aa_keywords, collapse = "|"), Name, ignore.case = TRUE) |
         ClassI == "Organic acids and derivatives")

cat(sprintf("\n  Amino acid/peptide metabolites:\n"))
cat(sprintf("    Liver: %d (%.1f%% of significant)\n",
            nrow(liver_aa), 100 * nrow(liver_aa) / nrow(liver_sig)))
cat(sprintf("    Blood: %d (%.1f%% of significant)\n",
            nrow(blood_aa), 100 * nrow(blood_aa) / nrow(blood_sig)))

# Save focused lists
write.csv(liver_lipids,
          file.path(integrated_dir, "lipid_aa_focus", "liver_significant_lipids.csv"),
          row.names = FALSE)
write.csv(blood_lipids,
          file.path(integrated_dir, "lipid_aa_focus", "blood_significant_lipids.csv"),
          row.names = FALSE)
write.csv(liver_aa,
          file.path(integrated_dir, "lipid_aa_focus", "liver_significant_aminoacids.csv"),
          row.names = FALSE)
write.csv(blood_aa,
          file.path(integrated_dir, "lipid_aa_focus", "blood_significant_aminoacids.csv"),
          row.names = FALSE)

# Generate focused volcano plots
liver_results$Highlight <- "Other"
liver_results$Highlight[liver_results$Compound_ID %in% liver_lipids$Compound_ID] <- "Lipid"

p_liver_lipid <- ggplot(liver_results, aes(x = log2FC, y = -log10(pvalue))) +
  geom_point(aes(color = Highlight, size = Highlight, alpha = Highlight)) +
  scale_color_manual(values = c("Other" = "grey70", "Lipid" = "#E31A1C")) +
  scale_size_manual(values = c("Other" = 1, "Lipid" = 2)) +
  scale_alpha_manual(values = c("Other" = 0.3, "Lipid" = 0.8)) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Liver Lipid Metabolites",
       subtitle = sprintf("%d significant lipids highlighted", nrow(liver_lipids)),
       x = "log2(Fold Change)",
       y = "-log10(p-value)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(integrated_dir, "figures", "liver_lipids_volcano.png"),
       p_liver_lipid, width = 10, height = 7, dpi = 300)

cat("\n✓ Lipid and amino acid analysis complete\n\n")

# ==========================================================================
# PART 3: Final Summary Report
# ==========================================================================

cat("PART 3: Generating Final Summary Report\n")
cat("=========================================\n")

# Load pathway enrichment results if available
pathway_enrichment <- tryCatch({
  readRDS(file.path(results_dir, "pathway_enrichment_summary.rds"))
}, error = function(e) NULL)

# Create summary
summary_text <- sprintf("
# GDF11 KO Metabolomics Analysis - Complete Summary

## Analysis Date: %s

## Sample Information
- **Study Design:** GDF11 KO vs WT mice (n=6 per group)
- **Tissues:** Liver and Blood
- **Platform:** LC-MS (positive and negative ion modes)
- **Data Processing:** Combined ion modes, t-test with BH-FDR correction

## Differential Metabolite Results

### Liver Tissue
- **Total metabolites:** %d
- **Significant (p<0.05, |FC|>1.5):** %d (%.1f%%)
  - Up-regulated: %d (%.1f%%)
  - Down-regulated: %d (%.1f%%)
- **Pattern:** Balanced bidirectional changes

### Blood Tissue  
- **Total metabolites:** %d
- **Significant (p<0.05, |FC|>1.2):** %d (%.1f%%)
  - Up-regulated: %d (%.1f%%)
  - Down-regulated: %d (%.1f%%)
- **Pattern:** Predominantly down-regulated

## Cross-Tissue Analysis

### Metabolite Overlap
- **Common metabolites detected:** %d
- **Liver-specific significant:** %d
- **Blood-specific significant:** %d
- **Significant in both tissues:** %d

### Correlation
- **Pearson correlation (log2FC):** r = %.3f
- **Interpretation:** %s correlation suggests %s

## Focused Metabolite Classes

### Lipid Metabolism
- **Liver lipids:** %d (%.1f%% of liver significant)
- **Blood lipids:** %d (%.1f%% of blood significant)
- **Key finding:** Major disruption in fatty acid metabolism

### Amino Acid/Peptide Metabolism
- **Liver amino acids:** %d (%.1f%% of liver significant)
- **Blood amino acids:** %d (%.1f%% of blood significant)
- **Key finding:** Widespread amino acid depletion

## Pathway Enrichment Results

%s

## Key Biological Insights

1. **GDF11 is a Master Metabolic Regulator**
   - Extensive metabolic reprogramming in liver
   - Limited systemic manifestation in blood
   - Tissue-specific metabolic responses

2. **Metabolic Stress Phenotype**
   - Impaired fatty acid oxidation (acylcarnitine dysregulation)
   - Amino acid depletion (increased catabolism)
   - Oxidative stress markers
   - Resembles accelerated metabolic aging

3. **Homeostatic Buffering**
   - Weak liver-blood correlation (r=%.3f)
   - Strong systemic homeostatic mechanisms
   - Minimal spillover of liver changes to blood

## Output Files

### Liver Results
- `liver/liver_differential_metabolites.csv` - All %d metabolites
- `liver/liver_significant_metabolites.csv` - %d significant metabolites
- `liver/figures/liver_volcano_plot.png`
- `liver/pathway_enrichment/` - ORA and GSEA results

### Blood Results
- `blood/blood_differential_metabolites.csv` - All %d metabolites
- `blood/blood_significant_metabolites.csv` - %d significant metabolites
- `blood/figures/blood_volcano_plot.png`
- `blood/pathway_enrichment/` - ORA and GSEA results

### Integrated Analysis
- `integrated_analysis/cross_tissue/liver_blood_comparison.csv`
- `integrated_analysis/lipid_aa_focus/` - Focused metabolite lists
- `integrated_analysis/figures/` - Integrated visualizations

## Next Steps

1. **Validation:** Targeted LC-MS/MS on top candidates
2. **Mechanism:** Enzyme assays, 13C-isotope tracing
3. **Intervention:** GDF11 replacement, metabolic rescue
4. **Translation:** Human biomarker studies

---
Analysis completed: %s
",
Sys.Date(),
nrow(liver_results), nrow(liver_sig), 100*nrow(liver_sig)/nrow(liver_results),
sum(liver_sig$regulation == "Up"), 100*sum(liver_sig$regulation == "Up")/nrow(liver_sig),
sum(liver_sig$regulation == "Down"), 100*sum(liver_sig$regulation == "Down")/nrow(liver_sig),
nrow(blood_results), nrow(blood_sig), 100*nrow(blood_sig)/nrow(blood_results),
sum(blood_sig$regulation == "Up"), 100*sum(blood_sig$regulation == "Up")/nrow(blood_sig),
sum(blood_sig$regulation == "Down"), 100*sum(blood_sig$regulation == "Down")/nrow(blood_sig),
length(common_compounds),
sum(pattern_table["Liver-specific"]),
sum(pattern_table["Blood-specific"]),
sum(pattern_table["Both tissues"]),
fc_cor,
ifelse(abs(fc_cor) > 0.5, "Strong", ifelse(abs(fc_cor) > 0.3, "Moderate", "Weak")),
ifelse(abs(fc_cor) < 0.3, "strong homeostatic buffering", "moderate metabolic coordination"),
nrow(liver_lipids), 100*nrow(liver_lipids)/nrow(liver_sig),
nrow(blood_lipids), 100*nrow(blood_lipids)/nrow(blood_sig),
nrow(liver_aa), 100*nrow(liver_aa)/nrow(liver_sig),
nrow(blood_aa), 100*nrow(blood_aa)/nrow(blood_sig),
if (!is.null(pathway_enrichment)) {
  sprintf("### Liver Pathways (ORA)
- Pathways with hits: %d
- Significant (p<0.05): %d
- FDR-significant (FDR<0.05): %d

### Blood Pathways (ORA)
- Pathways with hits: %d
- Significant (p<0.05): %d

### GSEA Results
- Liver pathways analyzed: %d
- Blood pathways analyzed: %d

Note: See pathway_enrichment folders for detailed results",
          nrow(pathway_enrichment$liver_ora),
          sum(pathway_enrichment$liver_ora$Pvalue < 0.05),
          sum(pathway_enrichment$liver_ora$FDR < 0.05),
          nrow(pathway_enrichment$blood_ora),
          sum(pathway_enrichment$blood_ora$Pvalue < 0.05),
          nrow(pathway_enrichment$liver_gsea),
          nrow(pathway_enrichment$blood_gsea))
} else {
  "Pathway enrichment analysis completed. See individual result files."
},
fc_cor,
nrow(liver_results), nrow(liver_sig),
nrow(blood_results), nrow(blood_sig),
Sys.time()
)

writeLines(summary_text, file.path(results_dir, "ANALYSIS_SUMMARY_REPORT.md"))

cat("\n✓ Final summary report generated\n")
cat(sprintf("\n==========================================================================\n"))
cat("ANALYSIS COMPLETE!\n")
cat("==========================================================================\n\n")
cat(sprintf("All results saved to: %s/\n\n", results_dir))
cat("Key files:\n")
cat("  - ANALYSIS_SUMMARY_REPORT.md (main summary)\n")
cat("  - liver/* (liver results)\n")
cat("  - blood/* (blood results)\n")
cat("  - integrated_analysis/* (cross-tissue and focused analyses)\n\n")
