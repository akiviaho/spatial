# 9.8.2021 Antti Kiviaho
#
# Script 3/3 in cell-type identification series of scripts
#
# This script is meant for cell-type identification through dimensionality reduction and 
# expression of marker genes in clusters identified through dr. The cell-type identities
# are inspected manually. The data has been published in
# "Single-cell analysis reveals transcriptomic remodellings in distinct cell types that 
# contribute to human prostate cancer progression" by Chen et al. in Nature: 
# https://www.nature.com/articles/s41556-020-00613-6.
#
# Cell-barcode/Cell-type matrix is outputted as plain text in order to make it easily
# accessible. This information can be used together with the expression data.
#
# The script follows the methods used in the original article as closely as possible.
# Uses R 3.5.1

# SETUP ####
library(Seurat)


setwd("~/Dropbox (Compbio)/prostate_spatial/data/single-cell-integration/")
results_path <- "~/Dropbox (Compbio)/prostate_spatial/results/single-cell-integration/"

data <- read.table("GSM4203181_data.matrix.txt")
data <- CreateSeuratObject(data,project="single-cell")

# Read Cell cluster identities and add annotation
clusters <- read.table("MNN-cluster-idents.txt")
data <- AddMetaData(object =  data, metadata =  clusters,col.name = "cluster")

# Read tSNE reduction
tsne <- read.table("MNN-cluster-tsne-reduction.txt")
data <- SetDimReduction(object = data, reduction.type = "tsne", slot = "cell.embeddings", 
                        new.data = as.matrix(tsne))
data <- SetDimReduction(object = data, reduction.type = "tsne", slot = "key", 
                        new.data = "tsne")

# Add original sample annotation
original <- data.frame(matrix(ncol=1,nrow=ncol(data@data)))
rownames(original) <- names(data@data)
colnames(original) <- c("original")
original$original <- gsub('.*\\.(.*)','\\1',rownames(original))
data <- AddMetaData(object = data, metadata = original, col.name = "original")

# Set marker genes to use (same as in articles fig 1b.)
marker_genes <- c("CD3D","CD3E","CD3G","CD2","CD7","DPP4","SLAMF1","CD27","CD28",
                  "CCR7","IL7R","IL2RA","BTLA","CD4","PTPRC","SELL","CD8A","CTLA4",
                  "PDCD1","HAVCR2","UCHL1","CD14","CD163","CD68","CSF1R","FCGR3A",
                  "LYZ","KRT5","KRT14","TP63","KRT8","KRT18","KRT19","AR","TPSB2",
                  "TPSAB1","MS4A2","CMA1","ENG","VWF","PECAM1","ACTA2")

## PLOTTING ####

# Plot initial clusters and gene expression per clusters
png(paste0(results_path,"gene-dotplot.png"),width = 1000, height = 400)
DotPlot(data,genes.plot = marker_genes,group.by = "cluster",x.lab.rot = T)
dev.off()

pdf(file=paste0(results_path,"tsne-clusters.pdf"))
DimPlot(data, reduction.use = 'tsne', pt.size = 1, vector.friendly = TRUE,  no.legend = FALSE, no.axes = FALSE,do.label = T,group.by = "cluster")
dev.off()

pdf(file=paste0(results_path,"tsne-original-samples.pdf"))
DimPlot(data, reduction.use = 'tsne', pt.size = 1, vector.friendly = TRUE,  no.legend = FALSE, no.axes = FALSE,do.label = F,group.by = "original")
dev.off()
# Check clustering in respect to the original sample


# Review the dotplot and filter accordingly. Edit these if necessary
clusters[which(clusters$cluster %in% c("0")),"cluster"] <- "T"
clusters[which(clusters$cluster %in% c("1")),"cluster"] <- "Endothelial"
clusters[which(clusters$cluster %in% c("2","3","4","5","6","7","8","9")),"cluster"] <- "Luminal"
clusters[which(clusters$cluster %in% c("10")),"cluster"] <- "Fibroblast"
clusters[which(clusters$cluster %in% c("11")),"cluster"] <- "Monolytic"
clusters[which(clusters$cluster %in% c("12")),"cluster"] <- "Mast"
clusters[which(clusters$cluster %in% c("13")),"cluster"] <- "Basal/intermediate"
clusters <- subset(clusters, cluster != c("14")) # Ended up removing cluster 14, since it is fibrobalsts/chondrocytes judging from DE-genes enrichment

data <- AddMetaData(object =  data, metadata =  clusters,col.name = "cluster")
data <- SubsetData(data,cells.use = rownames(clusters))

png(paste0(results_path,"gene-dotplot-tsne-annotated_final.png"),width = 1000, height = 400)
DotPlot(data,genes.plot = marker_genes,group.by = "cluster",x.lab.rot = T)
dev.off()

pdf(file=paste0(results_path,"tsne-clusters-annotated_final.pdf"))
DimPlot(data, reduction.use = 'tsne', pt.size = 1, vector.friendly = TRUE, no.legend = FALSE, no.axes = FALSE,do.label = T,group.by = "cluster")
dev.off()

write.table(clusters,file="chen-2021-cell-type-annotations.txt")
