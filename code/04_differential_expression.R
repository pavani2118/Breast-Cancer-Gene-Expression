
# ============================================================
# Project: Breast Cancer Gene Expression Analysis
# Dataset: GSE42568
# File: 04_differential_expression.R
# Purpose: Differential Gene Expression Analysis
# ============================================================


# 1. Load Required Packages

library(limma)
library(ggplot2)
library(pheatmap)
library(AnnotationDbi)
library(hgu133plus2.db)


# 2. Create Output Folders

dir.create("results/figures", recursive = TRUE,
           showWarnings = FALSE)

dir.create("results/tables", recursive = TRUE,
           showWarnings = FALSE)


# 3. Load Expression Data and Sample Information

cat("Loading expression data...\n")

expr <- readRDS("data/expression_data.rds")

sample_info <- read.csv(
  "data/sample_groups.csv",
  stringsAsFactors = FALSE
)

# Match sample information to expression columns
sample_info <- sample_info[
  match(colnames(expr), sample_info$Sample_ID),
]

if (any(is.na(sample_info$Sample_ID))) {
  stop("Sample information does not match expression data.")
}

group <- factor(
  sample_info$Group,
  levels = c("Normal", "Cancer")
)

if (any(is.na(group))) {
  stop("Missing or invalid sample groups found.")
}

cat("Expression data loaded successfully!\n")
cat("Total samples:", ncol(expr), "\n")
print(table(group))


# 4. Create Design Matrix

cat("\nCreating design matrix...\n")

design <- model.matrix(~ 0 + group)

colnames(design) <- c("Normal", "Cancer")

print(design[1:6, ])


# 5. Fit Linear Model

cat("\nFitting linear model...\n")

fit <- lmFit(expr, design)


# 6. Compare Cancer and Normal Groups

cat("\nCalculating differential expression...\n")

contrast_matrix <- makeContrasts(
  Cancer - Normal,
  levels = design
)

fit2 <- contrasts.fit(
  fit,
  contrast_matrix
)

fit2 <- eBayes(fit2)


# 7. Extract Differential Expression Results

results <- topTable(
  fit2,
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

results$Probe_ID <- rownames(results)


# 8. Add Gene Symbols

cat("\nAdding gene annotations...\n")

gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = results$Probe_ID,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

results$Gene_Symbol <- unname(
  gene_symbols[results$Probe_ID]
)


# 9. Select Significant Genes

significant_genes <- results[
  !is.na(results$adj.P.Val) &
    results$adj.P.Val < 0.05 &
    abs(results$logFC) > 1,
]

cat("\nDifferential expression analysis completed!\n")

cat("Total probes:", nrow(results), "\n")
cat("Significant probes:", nrow(significant_genes), "\n")


# 10. Save All Differential Expression Results

write.csv(
  results,
  "results/tables/all_differential_expression_results.csv",
  row.names = FALSE
)


# 11. Save Significant Genes

write.csv(
  significant_genes,
  "results/tables/significant_genes.csv",
  row.names = FALSE
)


# 12. Create Volcano Plot

results$Significance <- "Not Significant"

results$Significance[
  !is.na(results$adj.P.Val) &
    results$adj.P.Val < 0.05 &
    results$logFC > 1
] <- "Upregulated"

results$Significance[
  !is.na(results$adj.P.Val) &
    results$adj.P.Val < 0.05 &
    results$logFC < -1
] <- "Downregulated"

results$Significance <- factor(
  results$Significance,
  levels = c(
    "Upregulated",
    "Downregulated",
    "Not Significant"
  )
)

volcano_plot <- ggplot(
  results,
  aes(
    x = logFC,
    y = -log10(pmax(adj.P.Val, 1e-300)),
    color = Significance
  )
) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Volcano Plot: Cancer vs Normal",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

print(volcano_plot)

ggsave(
  "results/figures/volcano_plot.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# 13. Create Heatmap of Top 30 Significant Probes

if (nrow(significant_genes) >= 2) {
  
  top_genes <- head(
    significant_genes,
    30
  )
  
  heatmap_matrix <- expr[
    top_genes$Probe_ID,
    ,
    drop = FALSE
  ]
  
  # Row-wise scaling
  heatmap_matrix <- t(
    scale(t(heatmap_matrix))
  )
  
  # Replace missing values if any
  heatmap_matrix[!is.finite(heatmap_matrix)] <- 0
  
  annotation_col <- data.frame(
    Group = group
  )
  
  rownames(annotation_col) <- colnames(expr)
  
  png(
    "results/figures/top_genes_heatmap.png",
    width = 1400,
    height = 1200,
    res = 150
  )
  
  pheatmap(
    heatmap_matrix,
    annotation_col = annotation_col,
    show_colnames = FALSE,
    show_rownames = TRUE,
    main = "Top Differentially Expressed Probes"
  )
  
  dev.off()
  
} else {
  cat("Not enough significant probes for heatmap.\n")
}


# 14. Display Top Results

cat("\nTop 10 Differentially Expressed Probes:\n")

print(
  head(
    results[
      ,
      c(
        "Probe_ID",
        "Gene_Symbol",
        "logFC",
        "P.Value",
        "adj.P.Val"
      )
    ],
    10
  )
)


# 15. Final Summary

cat("\n====================================\n")
cat("DIFFERENTIAL EXPRESSION COMPLETED!\n")
cat("====================================\n")

cat("Total probes:", nrow(results), "\n")
cat("Significant probes:", nrow(significant_genes), "\n")

cat("\nGenerated files:\n")

print(
  list.files(
    "results",
    recursive = TRUE
  )
)
