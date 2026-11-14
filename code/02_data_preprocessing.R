
# ==========================================
# Project: Breast Cancer Gene Expression
# Step 02: Data Preprocessing
# Dataset: GSE42568
# ==========================================

# 1. Load required library
library(GEOquery)
library(limma)

# 2. Load downloaded dataset
eset <- readRDS("data/GSE42568_eset.rds")

# 3. Extract gene expression data
expr <- exprs(eset)

# 4. Extract sample information
pdata <- pData(eset)

# 5. Extract probe information
fdata <- fData(eset)

# 6. Display dimensions
cat("Expression Data Dimensions:\n")
print(dim(expr))

cat("Sample Information Dimensions:\n")
print(dim(pdata))

cat("Probe Information Dimensions:\n")
print(dim(fdata))

# 7. Check first six expression values
head(expr[, 1:6])

# 8. Check sample information
head(pdata)

# 9. Check probe information
head(fdata)

# 10. Save extracted data
saveRDS(expr, "data/expression_data.rds")
saveRDS(pdata, "data/sample_information.rds")
saveRDS(fdata, "data/probe_information.rds")

cat("Data extraction completed successfully!\n")
