#!/usr/bin/env Rscript
# Master script to run enhanced enrichment analysis
# Executes both enrichment and visualization scripts

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║     Enhanced Pathway Enrichment Analysis - Complete Suite     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

start_time <- Sys.time()

# Step 1: Enhanced enrichment analysis
cat("STEP 1: Running enhanced enrichment analysis...\n")
cat("────────────────────────────────────────────────────────────────\n")
source("/data1/dyh/GDF11/codes/enhanced_pathway_enrichment.R")

cat("\n")

# Step 2: Visualization
cat("STEP 2: Creating visualizations...\n")
cat("────────────────────────────────────────────────────────────────\n")
source("/data1/dyh/GDF11/codes/visualize_enhanced_enrichment.R")

# Calculate runtime
end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║                    Analysis Complete!                          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat(sprintf("\nTotal runtime: %.2f minutes\n", runtime))
cat(sprintf("Results directory: /data1/dyh/GDF11/metabo_results_GDF11KO/enhanced_enrichment/\n\n"))

cat("Key outputs:\n")
cat("  1. KEGG complete pathway annotations\n")
cat("  2. HMDB class enrichment (SuperClass, Class, SubClass)\n")
cat("  3. Lipidmaps enrichment (Category, Main Class, Sub Class)\n")
cat("  4. Comprehensive visualization suite\n")
cat("  5. IDs prepared for SMPDB and Reactome analysis\n")
cat("  6. Detailed summary report\n\n")
