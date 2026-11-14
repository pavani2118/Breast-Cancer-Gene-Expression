
# ==========================================
# Project: Breast Cancer Gene Expression Analysis
# Dataset: GSE42568
# Step 01: Data Collection
# ==========================================

# 1. Load required package
library(GEOquery)

# 2. Create data directory if it does not exist
if (!dir.exists("data")) {
  dir.create("data", recursive = TRUE)
}

# 3. Download GEO dataset
cat("Downloading GEO dataset...\n")

gse <- getGEO(
  "GSE42568",
  GSEMatrix = TRUE,
  getGPL = FALSE,
  destdir = file.path(getwd(), "data")
)

# 4. Check number of datasets
cat("Number of datasets:", length(gse), "\n")

# 5. Extract the first dataset
eset <- gse[[1]]

# 6. Save dataset locally
saveRDS(eset, "data/GSE42568_eset.rds")

# 7. Display dataset information
cat("Dataset downloaded successfully!\n")
cat("Dataset name:", annotation(eset), "\n")
cat("Number of probes and samples:", dim(eset), "\n")

# 8. Check downloaded files
print(list.files("data"))

