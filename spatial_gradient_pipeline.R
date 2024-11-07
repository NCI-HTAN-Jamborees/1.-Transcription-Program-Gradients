#----------------------------------------
# Project: Spatial Gradient Analysis 
# This script processes spatial transcriptomics data using Seurat, 
# applies non-negative matrix factorization (NMF) for dimensionality reduction, 
# performs LSGI analysis for local trajectory inference, and 
# conducts functional annotation using hallmark gene sets. 
# The goal is to analyze spatial patterns in breast cancer tissue samples and understand functional modules.
#----------------------------------------
# Requiring Libraries
library(Seurat)
library(Matrix)
library(RcppML)  # Ref: https://github.com/zdebruine/RcppML
library(ggplot2)
library(dplyr)
library(LSGI) # Ref: https://github.com/qingnanl/LSGI
library(hypeR)
library(msigdbr)

# Set seed for reproducibility
set.seed(2023)  

#----------------------------------------
# Section 1: Helper Functions
#----------------------------------------

# Function to clean gene list by removing certain genes like ribosomal and mitochondrial genes
clean_glist <- function(glist) {
  idx <- grep("^RPS|XIST|RP11|RPL|^MT|\\.", glist)
  glist <- glist[-idx]
  return(glist)
}

# Function to process Seurat object by filtering and subsetting features
sr.process <- function(object) {
  counts <- GetAssayData(object)
  counts <- counts[clean_glist(rownames(counts)), ]
  object <- subset(object, features = rownames(counts))
  return(object)
}

#----------------------------------------
# Section 2: Load and Prepare Data
#----------------------------------------

# Set input directory
InDir = '~/Projects/Data/HTAN_breast_cancer/HT305B1/HT305B1-S1H1Fc2U1Z1Bs1/'

# Load 10X data from directory
Xdata <- Seurat::Read10X(data.dir = paste(InDir, "filtered_feature_bc_matrix", sep = ""))

# Create Seurat object for spatial data
XF <- CreateSeuratObject(counts = Xdata, project = 'sample', assay = "Spatial")

# Read image files
Ximage <- Read10X_Image(image.dir = paste(InDir, "spatial", sep = ""))
Seurat::DefaultAssay(Ximage) <- "Spatial"

# Link matrix and image file
Ximage <- Ximage[colnames(XF)]
XF[["image"]] <- Ximage

# Get the barcodes that have image coordinates and subset Seurat object accordingly
barcodes_in_image <- rownames(XF@images$image@coordinates)
XF_clean <- subset(XF, cells = barcodes_in_image)

# Filter cells to keep those with more than 200 features (genes)
TumorST <- XF_clean[, XF_clean$nFeature_Spatial > 200]

#----------------------------------------
# Section 3: Preprocessing
#----------------------------------------

# Preprocess the Seurat object
TumorST <- sr.process(TumorST)
TumorST <- NormalizeData(TumorST)

#----------------------------------------
# Section 4: Non-negative Matrix Factorization (NMF)
#----------------------------------------

# Function to scan NMF ranks and calculate Mean Squared Error (MSE) for each rank
scan.nmf.mse <- function(obj, ranks = seq(1, 30, 2), tol = 1e-04) {
  dat <- obj@assays$Spatial$data
  errors <- c()
  for (i in ranks) {
    mod <- RcppML::nmf(dat, i, tol = tol, verbose = F)
    mse_i <- mse(dat, mod$w, mod$d, mod$h)
    errors <- c(errors, mse_i)
  }
  results <- data.frame(rank = ranks, MSE = errors)
  return(results)
}

# Function to perform NMF and add embeddings to Seurat object
sr.nmf <- function(obj, k = 10, tol = 1e-06, assay = "RNA") {
  dat <- obj@assays$Spatial$data
  nmf_model <- RcppML::nmf(dat, k = k, tol = tol, verbose = F)
  embeddings <- t(nmf_model$h)
  rownames(embeddings) <- colnames(obj)
  colnames(embeddings) <- paste0("nmf_", 1:k)
  loadings <- nmf_model$w
  rownames(loadings) <- rownames(obj)
  obj@reductions$nmf <- CreateDimReducObject(embeddings = embeddings,
                                             loadings = loadings, key = "nmf_", assay = assay)
  return(obj)
}

# Run NMF with k = 10
TumorST <- sr.nmf(obj = TumorST, k = 10, tol = 1e-05)

#----------------------------------------
# Section 5: Local Spatial Gradient Inference (LSGI) Analysis
#----------------------------------------

# Extract spatial coordinates and embeddings for LSGI
spatial_coords <- TumorST@images$image@coordinates[, c(4, 5)]
colnames(spatial_coords) <- c("X", "Y")
embeddings <- TumorST@reductions$nmf@cell.embeddings

# Run LSGI analysis
lsgi.res <- local.traj.preprocessing(spatial_coords = spatial_coords,
                                     n.grids.scale = 5, embeddings = embeddings, n.cells.per.meta = 25)

# Save processed data and LSGI results
saveRDS(TumorST, file = "HT206B1-S1Fc1U2Z1B1_processed_data.rds")
saveRDS(lsgi.res, file = "HT206B1-S1Fc1U2Z1B1_lsgi.rds")

#----------------------------------------
# Section 6: Visualization of NMF Gradients
#----------------------------------------

# Plot selected gradients with r_squared threshold
plot <- plt.factors.gradient.ind(info = lsgi.res, r_squared_thresh = 0.6,
                                 sel.factors = c("nmf_3", "nmf_4"), minimum.fctr = 10)

# Calculate average distance between NMF gradients
dist.mat <- avg.dist.calc(info = lsgi.res, r_squared_thresh = 0.6, minimum.fctr = 10)
plt.dist.heat(dist.mat)  # Plot distance heatmap

#----------------------------------------
# Section 7: Functional Annotation using HypeR
#----------------------------------------

# Perform functional annotation using HypeR
feature.loadings <- as.data.frame(TumorST@reductions$nmf@feature.loadings)
top.gene.list <- list()
n.genes.factor = 50
for (i in 1:ncol(feature.loadings)) {
  o <- order(feature.loadings[, i], decreasing = T)[1:n.genes.factor]
  features <- rownames(feature.loadings)[o]
  top.gene.list[[colnames(feature.loadings)[i]]] <- features
}

# Obtain hallmark gene sets
genesets <- msigdb_gsets("Homo sapiens", "H", clean=TRUE)

# Run hypergeometric test using HypeR
mhyp <- hypeR(top.gene.list, genesets, test="hypergeometric", background=30000)

# Save results to PDF
pdf(file="Hallmark_co_expression_modules_merge.pdf", h=8, w=12)
hyp_dots(mhyp, top = 50, merge=T, fdr=0.01, val="fdr", title="Co-expression Modules")
dev.off()

# Save results to Excel
hyp_to_excel(mhyp, file_path="MESO_fov19_Hallmark_hypeR.xlsx")
