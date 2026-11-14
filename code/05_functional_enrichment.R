
# ============================================================
# Project: Breast Cancer Gene Expression Analysis
# Dataset: GSE42568
# File: 05_functional_enrichment.R
# Purpose: GO Functional Enrichment Analysis
# ============================================================

# 1. Load Packages

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ggplot2)


# 2. Create Output Folders

dir.create("results/figures", recursive = TRUE,
           showWarnings = FALSE)

dir.create("results/tables", recursive = TRUE,
           showWarnings = FALSE)


# 3. Load Differential Expression Results

cat("Loading differential expression results...\n")

results <- read.csv(
  "results/tables/all_differential_expression_results.csv"
)


# 4. Map Probe IDs to Entrez Gene IDs

cat("Mapping probe IDs to Entrez IDs...\n")

entrez_ids <- mapIds(
  hgu133plus2.db,
  keys = as.character(results$Probe_ID),
  column = "ENTREZID",
  keytype = "PROBEID",
  multiVals = "first"
)

results$ENTREZID <- unname(
  entrez_ids[as.character(results$Probe_ID)]
)

# Save mapped results

write.csv(
  results,
  "results/tables/gene_annotation_results.csv",
  row.names = FALSE
)


# 5. Define Background Genes

background_genes <- unique(
  na.omit(results$ENTREZID)
)

cat("Total mapped background genes:",
    length(background_genes), "\n")


# 6. Select Upregulated and Downregulated Genes

significant <- results[
  !is.na(results$adj.P.Val) &
    results$adj.P.Val < 0.05 &
    abs(results$logFC) > 1,
]

upregulated <- significant[
  significant$logFC > 1,
]

downregulated <- significant[
  significant$logFC < -1,
]

up_genes <- unique(
  na.omit(upregulated$ENTREZID)
)

down_genes <- unique(
  na.omit(downregulated$ENTREZID)
)

cat("Upregulated mapped genes:",
    length(up_genes), "\n")

cat("Downregulated mapped genes:",
    length(down_genes), "\n")


# 7. GO Enrichment for Upregulated Genes

cat("\nAnalyzing upregulated genes...\n")

if (length(up_genes) >= 5) {
  
  up_go <- enrichGO(
    gene = up_genes,
    universe = background_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 5,
    readable = TRUE
  )
  
  up_go_table <- as.data.frame(up_go)
  
  write.csv(
    up_go_table,
    "results/tables/upregulated_GO_enrichment.csv",
    row.names = FALSE
  )
  
  if (nrow(up_go_table) > 0) {
    
    up_plot <- dotplot(
      up_go,
      showCategory = 15
    ) +
      ggtitle("GO Enrichment - Upregulated Genes")
    
    ggsave(
      "results/figures/upregulated_GO_plot.png",
      plot = up_plot,
      width = 10,
      height = 7,
      dpi = 300
    )
  }
  
} else {
  cat("Not enough mapped upregulated genes.\n")
}


# 8. GO Enrichment for Downregulated Genes

cat("\nAnalyzing downregulated genes...\n")

if (length(down_genes) >= 5) {
  
  down_go <- enrichGO(
    gene = down_genes,
    universe = background_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 5,
    readable = TRUE
  )
  
  down_go_table <- as.data.frame(down_go)
  
  write.csv(
    down_go_table,
    "results/tables/downregulated_GO_enrichment.csv",
    row.names = FALSE
  )
  
  if (nrow(down_go_table) > 0) {
    
    down_plot <- dotplot(
      down_go,
      showCategory = 15
    ) +
      ggtitle("GO Enrichment - Downregulated Genes")
    
    ggsave(
      "results/figures/downregulated_GO_plot.png",
      plot = down_plot,
      width = 10,
      height = 7,
      dpi = 300
    )
  }
  
} else {
  cat("Not enough mapped downregulated genes.\n")
}


# 9. Final Summary

cat("\n====================================\n")
cat("GO ENRICHMENT ANALYSIS COMPLETED!\n")
cat("====================================\n")

cat("Background genes:", length(background_genes), "\n")
cat("Upregulated genes:", length(up_genes), "\n")
cat("Downregulated genes:", length(down_genes), "\n")

cat("\nGenerated files:\n")

print(
  list.files(
    "results",
    recursive = TRUE
  )
)
