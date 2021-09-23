library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(GSVA)


samples <- readRDS("~/Dropbox (Compbio)/prostate_spatial/data/samples-seurat-object.rds")
de_genes_down <- read.csv("~/ctcf_depletion_downreg_genes.csv",header = T, row.names = 1)
de_genes_up <- read.csv("~/ctcf_depletion_upreg_genes.csv",header = T, row.names = 1)

sample <- samples$CRPC_278

sample <- SCTransform(sample,assay = "Spatial")

sample <- RunPCA(sample, assay = "SCT", verbose = FALSE)
sample <- FindNeighbors(sample, reduction = "pca", dims = 1:30)
sample <- FindClusters(sample, verbose = FALSE, resolution = 0.7)
sample <- RunUMAP(sample, reduction = "pca", dims = 1:30, n.neighbors = 20)

p1 <- DimPlot(sample,reduction = "umap")
p2 <- SpatialPlot(sample,pt.size.factor = 2.2, alpha = 0.8)

my_gene_sets <- list(CTCF_up_markers = rownames(de_genes_up), CTCF_down_markers = rownames(de_genes_down))

sct_data <- sample@assays$SCT@counts

############

res <- gsva(as.matrix(sct_data),my_gene_sets)
res <- as.data.frame(t(res))

res$top_upreg <- "no"
res$top_downreg <- "no"

res[which(res$CTCF_up_markers > quantile(res$CTCF_up_markers,0.9)), "top_upreg"] <- "yes"
res[which(res$CTCF_down_markers < quantile(res$CTCF_down_markers,0.1)), "top_downreg"] <- "yes"

cluster <- as.data.frame(sample$seurat_clusters)
res <- merge(res, cluster,by=0)
colnames(res)[6] <- "cluster"

#######
sample <- samples$CRPC_278
res <- gsva_crpc_278

pdf("CTCF_depletion_enriched_CRPC_278.pdf")
DimPlot(sample, reduction = "umap") + ggtitle("CRPC_278")
SpatialPlot(sample, images = "CRPC_278", label = T)  + ggtitle("CRPC_278")
DimPlot(sample,cells.highlight = paste0("CRPC_278_",res[which(res$top_upreg=="yes" & res$top_downreg == "yes"),"Row.names"])) + ggtitle("CRPC_278")
SpatialPlot(sample,cells.highlight = paste0("CRPC_278_",res[which(res$top_upreg=="yes" & res$top_downreg == "yes"),"Row.names"]),images = "CRPC_278") + ggtitle("CRPC_278")
dev.off()
# try to get barplot of categories

