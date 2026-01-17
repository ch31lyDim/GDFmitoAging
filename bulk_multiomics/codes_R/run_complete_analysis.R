#!/usr/bin/env Rscript
# Master Script - Run Complete Metabolomics Analysis Pipeline
# Author: Automated Analysis
# Date: 2024

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║     GDF11 KO METABOLOMICS ANALYSIS - COMPLETE PIPELINE          ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")
cat("\n")

start_time <- Sys.time()

# Set working directory
setwd("/data1/dyh/GDF11")

# Define scripts in order
scripts <- c(
  "codes/complete_metabolomics_analysis.R",
  "codes/pathway_enrichment_analysis.R",
  "codes/integrated_analysis.R"
)

# Run each script
for (i in seq_along(scripts)) {
  cat(sprintf("\n>>> STEP %d/%d: %s <<<\n\n", i, length(scripts), basename(scripts[i])))
  
  script_start <- Sys.time()
  
  tryCatch({
    source(scripts[i], echo = FALSE)
    script_end <- Sys.time()
    script_duration <- difftime(script_end, script_start, units = "mins")
    cat(sprintf("\n✓ Step %d completed in %.2f minutes\n", i, as.numeric(script_duration)))
  }, error = function(e) {
    cat(sprintf("\n✗ Error in Step %d: %s\n", i, e$message))
    stop(sprintf("Pipeline failed at step %d", i))
  })
}

end_time <- Sys.time()
total_duration <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                   PIPELINE COMPLETED SUCCESSFULLY                ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat(sprintf("Total analysis time: %.2f minutes\n", as.numeric(total_duration)))
cat(sprintf("Results directory: metabo_results_GDF11KO/\n"))
cat(sprintf("\nCheck ANALYSIS_SUMMARY_REPORT.md for complete results\n\n"))
