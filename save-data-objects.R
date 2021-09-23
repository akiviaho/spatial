library(Seurat)
library(hdf5r)
library(tidyr)
library(ggplot2)

# Paths and variables
data_folder_path <- "~/Dropbox (Compbio)/prostate_spatial/results/spaceranger-outs/data"
results_folder_path <- "~/Dropbox (Compbio)/prostate_spatial/"
setwd(results_folder_path)

sample_names <- list("BPH_651","BPH_688","CRPC_278","CRPC_489","CRPC_697","PC_03_6712","PC_15420OIK","PC_7875OIK","PC_4980")
resolution_param = 0.4
n_neighbor_UMAP = 30

# Implementation

samples <- lapply(file.path(data_folder_path,sample_names),Load10X_Spatial)
for(idx in 1:length(samples)){ Idents(samples[[idx]])<- sample_names[[idx]]}

# Normalize 
samples <- lapply(samples,SCTransform,assay="Spatial",variable.features.n=500)
var_features <- lapply(samples,VariableFeatures)
var_features_combined <- unique(unlist(var_features))

# Aggregate samples to one bulk
samples_aggregated <- merge(samples[[1]], y = unlist(samples)[-1], add.cell.ids = unlist(sample_names), project = "aggregated")
samples_aggregated$orig.sample <- Idents(samples_aggregated)
names(samples_aggregated@images) <- sample_names
VariableFeatures(samples_aggregated) <- var_features_combined

# Run dim reduction on samples
samples_aggregated <- RunPCA(samples_aggregated, assay = "SCT", verbose = FALSE)
samples_aggregated <- FindNeighbors(samples_aggregated, reduction = "pca", dims = 1:PCA_dims)
samples_aggregated <- FindClusters(samples_aggregated, verbose = FALSE, resolution = resolution_param)
samples_aggregated <- RunUMAP(samples_aggregated, reduction = "pca", dims = 1:UMAP_dims, n.neighbors = n_neighbor_UMAP)

for(i in 1:length(samples)){DefaultAssay(samples[[i]]) <- "Spatial"}

saveRDS(samples,file="samples-seurat-object.rds")
saveRDS(samples_aggregated,file="aggregated_samples_seurat_obj.rds")

