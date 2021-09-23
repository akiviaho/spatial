# Antti Kiviaho 26.7.2021
#
# Script 2/3 in cell-type identification series of scripts
#
# This script is meant for unsupervised clustering of Mutual Nearest Neighbor (MNN)
# transformed counts generated from log-normalized count matrix 
# published in "Single-cell analysis reveals transcriptomic remodellings
# in distinct cell types that contribute to human prostate cancer progression" by
# Chen et al. in Nature: https://www.nature.com/articles/s41556-020-00613-6
#
# The script writes cluster, pca and tsne reductions per cell for all cells in the 
# dataset. These can be further used in cell-type identification:
# script 3/3, cluster-recognition.R
#
# The script follows the methods used in the original article as closely as possible.
# Uses R 3.5.1

library(Seurat)

setwd("~/Dropbox (Compbio)/prostate_spatial/data/chen-sc-reference/runs/")
results_path <- "~/Dropbox (Compbio)/prostate_spatial/results/chen-sc-reference/"

data <- read.table(file = "~/Dropbox (Compbio)/prostate_spatial/data/chen-sc-reference/MNN-corrected-counts.txt",header = TRUE)
data <- CreateSeuratObject(data,project="single-cell")

# Run clustering
data <- ScaleData(object = data)
data <- FindVariableGenes(object= data,  mean.function = ExpMean, dispersion.function = LogVMR)
data <- RunPCA(object =data, pcs.compute = 8, do.print = F, seed.use = 42)
data <- ProjectPCA(object = data, do.print = FALSE)
data <- FindClusters(object = data, reduction.type = "pca",dims.use = 1:8, resolution = 0.1)
data <- RunTSNE(object = data, reduction.use = "pca",  dims.use = 1:8, perplexity = 30 ,do.fast = TRUE)
# 
saveRDS(data,"tSNE-reduced-MNN-seurat-obj.rds")

# Plotting ####
pdf(file="chen-sc-tsne-plot_v3.pdf")
DimPlot(data, reduction.use = 'tsne', pt.size = 1, vector.friendly = TRUE, label.size = 12, no.legend = FALSE, no.axes = FALSE,do.label = T)
dev.off()

idents <- data.frame(data@ident)
colnames(idents) <- c("barcode","cluster")
write.table(idents,file="MNN-cluster-idents.txt")

pca <- data.frame(GetCellEmbeddings(data, reduction.type = "pca"))
write.table(pca,file="MNN-cluster-pca-reduction.txt")

tsne <- data.frame(GetCellEmbeddings(data, reduction.type = "tsne"))
write.table(tsne,file="MNN-cluster-tsne-reduction.txt")
