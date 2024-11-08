# === Load Necessary Libraries ===
# These libraries are essential for data handling, visualization, and distance computation
library("Seurat")        # For handling Seurat objects
library("viridis")       # Color palette for visualizations
library(ggplot2)          # For generating visual plots
library(ComplexHeatmap)   # For generating heatmaps
library(proxy)            # For computing distance metrics
library(reshape2)         # For reshaping data frames

# === Load Data ===
# Load pre-processed Seurat objects for different samples
HT206B1 <- readRDS("./output/HT206B1-S1Fc1U2Z1B1_processed_data.rds")
HT268B1 <- readRDS("./output/HT268B1-Th1H3Fc2U12Z1Bs1_processed_data.rds")
HT271B1 <- readRDS("./output/HT271B1-S1H3Fc2U1Z1Bs1_processed_data.rds")
HT305B1 <- readRDS("./output/HT305B1-S1H1Fc2U1Z1Bs1_processed_data.rds")
HT323B1 <- readRDS("./output/HT323B1-S1H1Fc2U1Z1Bs1_processed_data.rds")
HT243B1H3A2 <- readRDS("./output/HT243B1H3A2-S1Fc1U1Z1B1_processed_data.rds")

# === Function Definitions ===
# Define functions used for analysis and feature extraction

# Function to extract the top N genes based on NMF feature loadings
get.nmf.info <- function(obj, top.n = top.n) {
  # Extract NMF feature loadings from the Seurat object
  feature.loadings <- as.data.frame(obj@reductions$nmf@feature.loadings)
  
  # Initialize an empty list to store the top genes
  top.gene.list <- list()
  
  # Loop through each component of the feature loadings
  for (i in 1:ncol(feature.loadings)) {
    # Order genes by their loading values (descending) and select the top N
    o <- order(feature.loadings[, i], decreasing = T)[1:top.n]
    tmp_mat <- feature.loadings[o, i, drop = F]
    idx <- which(tmp_mat > 0.0000001)  # Filter genes based on threshold
    features <- rownames(tmp_mat)[idx]  # Extract gene names
    top.gene.list[[colnames(feature.loadings)[i]]] <- features  # Store genes
  }
  
  # Return the list of top genes
  return(top.gene.list)
}

# === Sample Grouping ===
# Define sample lists for basal and luminal subtypes
list_of_samples_basal <- c(HT206B1 = HT206B1, HT268B1 = HT268B1, HT271B1 = HT271B1)
list_of_samples_luminal <- c(HT305B1 = HT305B1, HT323B1 = HT323B1, HT243B1H3A2 = HT243B1H3A2)

# === Feature Extraction ===
# Extract top genes for basal and luminal samples

# Basal Subtype: Extract top genes
list_of_vectors_basal <- unlist(lapply(list_of_samples_basal, get.nmf.info, 50), recursive = F)
names(list_of_vectors_basal) <- paste0(names(list_of_vectors_basal), paste0("_", "Basal", sep = ''))

# Luminal Subtype: Extract top genes
list_of_vectors_luminal <- unlist(lapply(list_of_samples_luminal, get.nmf.info, 50), recursive = F)
names(list_of_vectors_luminal) <- paste0(names(list_of_vectors_luminal), paste0("_", "Luminal", sep = ''))

# === Jaccard Distance Computation ===
# Function to compute Jaccard distance between two sets
jaccard_dist <- function(x, y) {
  1 - (length(intersect(x, y)) / length(union(x, y)))
}

# Compute Jaccard distance matrices for basal and luminal subtypes
dist_matrix_B <- proxy::dist(list_of_vectors_basal, method = jaccard_dist)
dist_matrix_L <- proxy::dist(list_of_vectors_luminal, method = jaccard_dist)

# Convert distance objects to matrices for visualization
dist_matrix_B <- as.matrix(dist_matrix_B)
dist_matrix_L <- as.matrix(dist_matrix_L)

# === Visualization: Heatmap Generation ===
# Generate heatmaps for the Jaccard distance matrices
pdf(file = "./output/heatmap.pdf", h = 20, w = 20)
Heatmap(dist_matrix_B, col = colorRampPalette(c("red", "white", "blue"))(100))
Heatmap(dist_matrix_L, col = colorRampPalette(c("red", "white", "blue"))(100))
dev.off()

# === Cross-Subtype Jaccard Distance ===
# Compute Jaccard distances between basal and luminal subtype vectors
result.df.jacc <- expand.grid(names(list_of_vectors_basal), names(list_of_vectors_luminal))
colnames(result.df.jacc) <- c("Basal", "Luminal")
results_hyp <- apply(result.df.jacc, 1, function(idx) {
  jaccard_dist(list_of_vectors_basal[[idx[1]]], list_of_vectors_luminal[[idx[2]]])
}) 
result.df.jacc$jacc <- results_hyp

# Reshape the data frame to a matrix for heatmap visualization
result.mat.jacc <- reshape2::acast(data = result.df.jacc, Basal ~ Luminal)

# Generate a heatmap for the comparison between basal and luminal subtypes
pdf(file = "./output/heatmap_BasalvsLuminal.pdf", h = 20, w = 20)
Heatmap(result.mat.jacc, col = colorRampPalette(c("red", "white", "blue"))(100))
dev.off()
