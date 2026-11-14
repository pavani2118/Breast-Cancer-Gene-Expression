
# ==========================================
# 06 - FINAL RESULTS SUMMARY
# Breast Cancer Gene Expression Analysis
# ==========================================

# Set project directory
project_dir <- file.path(
  Sys.getenv("USERPROFILE"),
  "Downloads",
  "Breast-Cancer-Gene-Expression"
)

setwd(project_dir)

# Load results
results <- read.csv(
  "results/tables/all_differential_expression_results.csv"
)

significant_genes <- read.csv(
  "results/tables/significant_genes.csv"
)

up_go <- read.csv(
  "results/tables/upregulated_GO_enrichment.csv"
)

down_go <- read.csv(
  "results/tables/downregulated_GO_enrichment.csv"
)

sample_counts <- read.csv(
  "results/tables/sample_counts.csv"
)

# Calculate differential expression counts
total_probes <- nrow(results)
significant_count <- nrow(significant_genes)

upregulated_count <- sum(
  significant_genes$logFC > 1,
  na.rm = TRUE
)

downregulated_count <- sum(
  significant_genes$logFC < -1,
  na.rm = TRUE
)

# Create overall summary
summary_table <- data.frame(
  Analysis = c(
    "Total probes",
    "Significant probes",
    "Upregulated probes",
    "Downregulated probes",
    "Upregulated GO terms",
    "Downregulated GO terms"
  ),
  Result = c(
    total_probes,
    significant_count,
    upregulated_count,
    downregulated_count,
    nrow(up_go),
    nrow(down_go)
  )
)

# Display summary
print(summary_table)

# Save summary
write.csv(
  summary_table,
  "results/tables/final_results_summary.csv",
  row.names = FALSE
)

cat("\nFinal results summary saved successfully!\n")
