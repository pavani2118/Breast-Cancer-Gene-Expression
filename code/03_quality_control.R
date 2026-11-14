
# ============================================================
# Project: Breast Cancer Gene Expression Analysis
# Dataset: GSE42568
# File: 03_quality_control.R
# Purpose: Quality Control Analysis
# ============================================================


# 1. Load Required Packages

library(ggplot2)
library(pheatmap)


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

cat("Data loaded successfully!\n")

cat("Expression dimensions:\n")
print(dim(expr))


# 4. Match Sample Information with Expression Data

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
  stop("Some samples have missing group information.")
}


# 5. Check Missing Values

cat("\nChecking missing values...\n")

missing_values <- sum(is.na(expr))

cat("Total missing values:", missing_values, "\n")

write.csv(
  data.frame(Missing_Values = missing_values),
  "results/tables/missing_values.csv",
  row.names = FALSE
)


# 6. Boxplot - Expression Distribution

cat("\nCreating boxplot...\n")

png(
  "results/figures/expression_boxplot.png",
  width = 1600,
  height = 900,
  res = 150
)

boxplot(
  expr,
  outline = FALSE,
  las = 2,
  cex.axis = 0.5,
  main = "Gene Expression Distribution Across Samples",
  xlab = "Samples",
  ylab = "Expression Values"
)

dev.off()


# 7. Density Plot - Expression Distribution

cat("Creating density plot...\n")

png(
  "results/figures/expression_density.png",
  width = 1200,
  height = 800,
  res = 150
)

plot(
  density(expr[, 1], na.rm = TRUE),
  main = "Gene Expression Density Plot",
  xlab = "Expression Values",
  ylab = "Density",
  col = 1,
  lwd = 2
)

for (i in 2:ncol(expr)) {
  lines(
    density(expr[, i], na.rm = TRUE),
    lwd = 1
  )
}

dev.off()


# 8. Sample Correlation Analysis

cat("Calculating sample correlations...\n")

sample_correlation <- cor(
  expr,
  use = "pairwise.complete.obs",
  method = "pearson"
)

write.csv(
  sample_correlation,
  "results/tables/sample_correlation.csv"
)


# 9. Correlation Heatmap

cat("Creating correlation heatmap...\n")

annotation <- data.frame(
  Group = group
)

rownames(annotation) <- colnames(expr)

png(
  "results/figures/sample_correlation_heatmap.png",
  width = 1400,
  height = 1200,
  res = 150
)

pheatmap(
  sample_correlation,
  annotation_col = annotation,
  annotation_row = annotation,
  main = "Sample Correlation Heatmap",
  show_rownames = FALSE,
  show_colnames = FALSE
)

dev.off()


# 10. Principal Component Analysis (PCA)

cat("Performing PCA...\n")

# Select 500 most variable probes

probe_variance <- apply(
  expr,
  1,
  var,
  na.rm = TRUE
)

valid_probes <- which(
  is.finite(probe_variance) & probe_variance > 0
)

top_probes <- valid_probes[
  order(
    probe_variance[valid_probes],
    decreasing = TRUE
  )
]

top_probes <- head(top_probes, 500)

pca_data <- t(expr[top_probes, , drop = FALSE])

# Perform PCA
pca_result <- prcomp(
  pca_data,
  center = TRUE,
  scale. = FALSE
)

# Calculate percentage of variance
variance_percent <- (
  pca_result$sdev^2 /
    sum(pca_result$sdev^2)
) * 100

pca_plot_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = group,
  Sample = colnames(expr)
)


# Create PCA Plot

pca_plot <- ggplot(
  pca_plot_data,
  aes(
    x = PC1,
    y = PC2,
    color = Group
  )
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA of Breast Cancer Gene Expression Data",
    x = paste0(
      "PC1 (",
      round(variance_percent[1], 2),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_percent[2], 2),
      "%)"
    )
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

print(pca_plot)

ggsave(
  "results/figures/pca_plot.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# 11. Save PCA Results

write.csv(
  pca_plot_data,
  "results/tables/pca_results.csv",
  row.names = FALSE
)


# 12. Final Output

cat("\n====================================\n")
cat("QUALITY CONTROL COMPLETED!\n")
cat("====================================\n")

cat("\nTotal samples:", ncol(expr), "\n")
cat("Total probes:", nrow(expr), "\n")
cat("Missing values:", missing_values, "\n")

cat("\nGenerated files:\n")

print(
  list.files(
    "results",
    recursive = TRUE
  )
)
