
# ============================================================
# Project: Breast Cancer Gene Expression Analysis
# Dataset: GSE42568
# File: 02_data_preprocessing.R
# Purpose: Data preprocessing and sample classification
# ============================================================


# 1. Load Required Packages
library(GEOquery)
library(limma)
library(ggplot2)


# 2. Create Required Folders
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)


# 3. Load Previously Extracted Data
cat("Loading dataset...\n")

eset <- readRDS("data/GSE42568_eset.rds")

expr <- exprs(eset)
pdata <- pData(eset)

cat("Dataset loaded successfully!\n")

cat("Expression data dimensions:\n")
print(dim(expr))

cat("Sample information dimensions:\n")
print(dim(pdata))


# 4. Extract Tissue Information
cat("\nExtracting tissue information...\n")

tissue <- tolower(trimws(as.character(pdata[["tissue:ch1"]])))

# Classify samples into Normal and Cancer groups
group <- ifelse(
  grepl("normal", tissue), "Normal",
  ifelse(grepl("cancer", tissue), "Cancer", NA)
)

group <- factor(group, levels = c("Normal", "Cancer"))


# 5. Check Sample Classification
cat("\nSample classification:\n")

print(table(group, useNA = "ifany"))

cat("\nNumber of unclassified samples:\n")
print(sum(is.na(group)))

if (any(is.na(group))) {
  stop("Some samples could not be classified. Check tissue information.")
}


# 6. Create Sample Information Table
sample_info <- data.frame(
  Sample_ID = rownames(pdata),
  Tissue = tissue,
  Group = group
)

head(sample_info)

# Save sample information
write.csv(
  sample_info,
  "data/sample_groups.csv",
  row.names = FALSE
)

cat("\nSample information saved successfully!\n")


# 7. Count Samples in Each Group
sample_counts <- as.data.frame(table(group))

colnames(sample_counts) <- c("Group", "Count")

print(sample_counts)

# Save sample counts
write.csv(
  sample_counts,
  "results/tables/sample_counts.csv",
  row.names = FALSE
)


# 8. Create Sample Distribution Graph
sample_plot <- ggplot(
  sample_counts,
  aes(x = Group, y = Count, fill = Group)
) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(
    title = "Distribution of Breast Tissue Samples",
    x = "Sample Group",
    y = "Number of Samples"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    legend.position = "none"
  )

# Display graph in RStudio
print(sample_plot)


# 9. Save Graph as PNG Image
ggsave(
  filename = "results/figures/sample_distribution.png",
  plot = sample_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cat("\nGraph saved successfully!\n")


# 10. Check Expression and Sample Alignment
cat("\nChecking expression and sample alignment...\n")

if (!all(colnames(expr) %in% rownames(pdata))) {
  stop("Expression samples and sample information do not match.")
}

cat("Expression and sample information are aligned.\n")


# 11. Save Processed Data
saveRDS(
  expr,
  "data/expression_data.rds"
)

saveRDS(
  sample_info,
  "data/sample_information_processed.rds"
)

cat("\nProcessed data saved successfully!\n")


# 12. Display Final Results
cat("\n====================================\n")
cat("DATA PREPROCESSING COMPLETED!\n")
cat("====================================\n")

cat("\nTotal number of samples:", ncol(expr), "\n")

cat("Normal samples:",
    sum(group == "Normal"), "\n")

cat("Cancer samples:",
    sum(group == "Cancer"), "\n")

cat("\nGenerated files:\n")

print(list.files(
  "results",
  recursive = TRUE
))
