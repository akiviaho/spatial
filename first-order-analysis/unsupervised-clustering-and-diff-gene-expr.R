# Antti Kiviaho 20.7.2021
# This script is intended for the unsupervised clustering of a Spatial Transcriptomics
# dataset generated with 10X Visium technology.
# The script outputs the following plots:
#
# 1. UMAP of the unfiltered data on a 2D graph with labeled clusters
# 2. Spatial plots of the samples where the clusters are visible
# 3. UMAP of the DE-genes FILTERED data on a 2D graph with labeled clusters
# 4. Spatial plots of the samples where the DE-genes filtered clusters are visible
#
# All of the plots above are also produced in such way that the clustering is done
# sample-independently after normalizing across all samples.
#
# Additionally the script outputs an txt file with the DE-genes (padj < 0.05), abs(lFC) > 1
# Associated with each cluster

## SETUP ####

# Packages
library(Seurat)
library(hdf5r)
library(tidyr)
library(ggplot2)

# Paths
data_folder_path <- "~/Dropbox (Compbio)/prostate_spatial/results/spaceranger-outs"
results_folder_path <- "~/Dropbox (Compbio)/prostate_spatial/results/secondary-analysis-in-r"
gene_lengths_file<- "~/Dropbox (Compbio)/prostate_spatial/data/reference/gene_lengths.txt"
setwd(results_folder_path)

# Variables

sample_names <- list("BPH_651","BPH_688","CRPC_278","CRPC_489","CRPC_697","PC_03_6712","PC_15420OIK","PC_7875OIK","PC_4980")
sample_types <- list("BPH","BPH", "CRPC","CRPC","CRPC","PC","PC","PC","PC")
resolution_param = 0.4
n_neighbor_UMAP = 30

# Functions

do.dimensionality.reduction <- function(data, PCA_dims = 30, UMAP_dims = 30)
{
  data <- RunPCA(data, assay = "SCT", verbose = FALSE)
  data <- FindNeighbors(data, reduction = "pca", dims = 1:PCA_dims)
  data <- FindClusters(data, verbose = FALSE, resolution = resolution_param)
  data <- RunUMAP(data, reduction = "pca", dims = 1:UMAP_dims, n.neighbors = n_neighbor_UMAP)
  return(data)
}

calculate.variable.genes <- function(data)
{
  df_result <- data.frame()
  for(idx in 1:length(levels(Idents(data))))
  {
    df <- FindMarkers(data, ident.1 = levels(Idents(data))[idx] , min.pct = 0.25)
    df <- subset(df,avg_log2FC > 1 | avg_log2FC < -1 & p_val_adj < 1e-02)
    df$cluster <- rep(levels(Idents(data))[idx],nrow(df))
    df <- tibble::rownames_to_column(df,var = "SYMBOL")
    df_result <- rbind(df_result,df)
  }
  return(df_result)
}

plot.UMAP.and.spatial.aggregated <- function(data, names = sample_names)
{
  plot(DimPlot(data,label = T))
  plot(DimPlot(data,group.by = 'orig.sample', label = T))
  for(idx in 1:length(sample_names))
  {
    plot(SpatialPlot(data,images = names[[idx]],label=T) + ggtitle(names[[idx]]))
  }
}

plot.UMAP.and.spatial.individual <- function(data, data_aggregated = samples_aggregated, do_filter = F, 
                                             names = sample_names, genelist = DE_genes_per_cluster$SYMBOL)
{
  for(idx in 1:length(names))
  {
    data[[idx]] <- subset(data_aggregated,subset = orig.sample == names[[idx]])
    if(do_filter)
    {
      data[[idx]] <- subset(data[[idx]],features = genelist)
    }
    data[[idx]] <- do.dimensionality.reduction(data[[idx]])
   plot(DimPlot(data[[idx]]) + ggtitle(names[[idx]]))
   plot(SpatialPlot(data[[idx]],images = names[[idx]],label=T) + ggtitle(names[[idx]]))
  }
  return(data)
}

## RUN ####

# Import data
samples <- lapply(file.path(data_folder_path,sample_names),Load10X_Spatial)
for(idx in 1:length(samples)){ Idents(samples[[idx]])<- sample_names[[idx]]}

# Normalize and get common variable features 
samples <- lapply(samples,SCTransform,assay="Spatial",variable.features.n=500)
var_features <- lapply(samples,VariableFeatures)
var_features_combined <- unique(unlist(var_features))
# 
# features <- SelectIntegrationFeatures(object.list = samples, nfeatures = 3000)
# samples <- PrepSCTIntegration(object.list = samples, anchor.features = features)
# 
# anchors <- FindIntegrationAnchors(object.list = samples, normalization.method = "SCT",
#                                          anchor.features = features)
# samples_aggregated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# samples_aggregated$orig.sample <- Idents(samples_aggregated)
# names(samples_aggregated@images) <- sample_names


# Aggregate and set identifiers
samples_aggregated <- merge(samples[[1]], y = unlist(samples)[-1], add.cell.ids = unlist(sample_names), project = "aggregated")
samples_aggregated$orig.sample <- Idents(samples_aggregated)
names(samples_aggregated@images) <- sample_names
VariableFeatures(samples_aggregated) <- var_features_combined


# Normalize, reduce and calculate variable features between clusters
samples_aggregated <- do.dimensionality.reduction(samples_aggregated)

DE_genes_per_cluster <- calculate.variable.genes(samples_aggregated) # This takes time
write.table(DE_genes_per_cluster,file = "DE-genes-per-cluster.txt") # Save the table
write.csv(unique(DE_genes_per_cluster$SYMBOL),file = "DE-genes-per-cluster-names.csv")


# Dimensionality reduction and plotting individually 

pdf(file="UMAP-and-spatial-individually-unfiltered.pdf")
plot.UMAP.and.spatial.individual(samples)
dev.off()


pdf(file="UMAP-and-spatial-individually-DE-genes-filtered.pdf")
plot.UMAP.and.spatial.individual(samples, samples_aggregated, do_filter = TRUE)
dev.off()

# Dimensionality reduction and plotting together
pdf(file="UMAP-and-Spatial-cluster-unfiltered.pdf")
plot.UMAP.and.spatial.aggregated(samples_aggregated)
dev.off()


# Filter samples and plot them
filtered_samples_aggregated <- subset(samples_aggregated,features = DE_genes_per_cluster$SYMBOL)
filtered_samples_aggregated <- do.dimensionality.reduction(filtered_samples_aggregated)

pdf(file="UMAP-and-Spatial-cluster-DE-genes-filtered.pdf")
plot.UMAP.and.spatial.aggregated(filtered_samples_aggregated)
dev.off()



