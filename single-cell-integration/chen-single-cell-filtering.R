# Antti Kiviaho 26.7.2021
# This script is meant for QC, filtering and cell-type signature generation from 
# the dataset published in "Single-cell analysis reveals transcriptomic remodellings
# in distinct cell types that contribute to human prostate cancer progression" by
# Chen et al. in Nature: https://www.nature.com/articles/s41556-020-00613-6
#
# The script follows the methods used in the original article as closely as possible.

library(scran)
library(Seurat)
library(qusage)
library(infercnv)

results_path <- "~/Dropbox (Compbio)/prostate_spatial/results/single-cell-integration/"

data <- read.table(file = "~/Dropbox (Compbio)/prostate_spatial/data/single-cell-integration/MNN-corrected-counts.txt",header = TRUE)
data <- CreateSeuratObject(data,project="single-cell")

# Run clustering
data <- ScaleData(object = data)
data <- FindVariableGenes(object= data,  mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = -1, x.high.cutoff = Inf, y.cutoff = 1, do.plot = F)

data <- RunPCA(object =data, pcs.compute = 20, do.print = F, seed.use = 42)
data <- ProjectPCA(object = data, do.print = FALSE)
data <- FindClusters(object = data, reduction.type = "pca", 
                           dims.use = 1:8, print.output = FALSE, save.SNN = TRUE, force.recalc = TRUE);
data <- RunTSNE(object = data, reduction.use = "pca",  dims.use = 1:8, perplexity = 30 ,do.fast = TRUE);



# Plotting ####
fontsize <- theme(axis.text=element_text(size=48, colour = 'black'), axis.title=element_text(size=48), legend.text = element_text(size = 24), legend.title = element_text(size = 0), 
                  plot.margin=unit(c(.5,.5,1,.5),"cm"), axis.text.x = element_text(angle = 0, hjust = 1), 
                  panel.border = element_blank(), axis.line = element_line(colour = "black", size = rel(1)));

pdf(file="~/Dropbox (Compbio)/prostate_spatial/results/single-cell-integration/chen-sc-tsne-plot_v2.pdf")
DimPlot(data, reduction.use = 'tsne', pt.size = 1, vector.friendly = TRUE, label.size = 24, no.legend = FALSE, no.axes = FALSE) + 
 + fontsize + labs(x = 'Dimension 1', y = 'Dimension 2'); # group.by = "type"
dev.off()

# saveRDS(data,file="~/Dropbox (Compbio)/prostate_spatial/results/single-cell-integration/tsne-seurat-object.rds")

marker_genes <- c("CD3D","CD3E","CD3G","CD2","CD7","DPP4","SLAMF1","CD27","CD28",
                 "CCR7","IL7R","IL2RA","BTLA","CD4","PTPRC","SELL","CD8A","CTLA4",
                 "PDCD1","HAVCR2","UCHL1","CD14","CD163","CD68","CSF1R","FCGR3A",
                 "LYZ","KRT5","KRT14","TP63","KRT8","KRT18","KRT19","AR","TPSB2",
                 "TPSAB1","MS4A2","CMA1","ENG","VWF","PECAM1","ACTA2")
