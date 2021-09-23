
library(Seurat)
library(ggplot2)

setwd("~/Dropbox (Compbio)/prostate_spatial/data/chen-sc-reference/")

data = read.table("GSM4203181_data.matrix.txt",sep='\t', header = T, row.names = 1)
chen = CreateSeuratObject(counts = data, project = "chen_et_al", min.cells = 3, min.features = 200)
rm(data)
chen = ScaleData(object = chen)
chen = FindVariableFeatures(object= chen,  mean.function = ExpMean, dispersion.function = LogVMR)
chen = RunPCA(object = chen, npcs = 20)
chen = RunTSNE(object = chen,dims=1:8)
chen = FindNeighbors(chen)
chen = FindClusters(object = chen, resolution = 0.1)
P = DimPlot(chen,label = T)

marker_genes_chen = c("CD3D","CD3E","CD3G","CD2","CD7","DPP4","SLAMF1","CD27","CD28","CCR7","IL7R","IL2RA","BTLA","CD4","PTPRC","SELL","CD8A","CTLA4","PDCD1","HAVCR2","UCHL1","CD14","CD163","CD68","CSF1R","FCGR3A","LYZ","KRT5","KRT14","TP63","KRT8","KRT18","KRT19","AR","TPSB2","TPSAB1","MS4A2","CMA1","ENG","VWF","PECAM1","ACTA2")
P = DotPlot(chen,features = marker_genes_chen, dot.scale = 8) + RotatedAxis() + labs(y="Cell cluster",x="Genes")
ggsave("dotplot_chen_etal_clusters.pdf",P,device="pdf",width=15,height = 4)

P = DimPlot(object = chen, reduction = "tsne")
ggsave("tsneplot_chen_etal_clusters.pdf",P,device="pdf",width=15,height = 15)

#Clusters 0, 1, 2, 3, 6 = Luminal
#Cluster 4 = Endothelial
#Cluster 5 = T cell
#Cluster 7 = Monocytic
#Cluster 8 = Fibroblast
#Cluster 9 = Mast
#Cluster 10 = Basal/intermediate

chen <- RenameIdents(
  object = chen,
  '0' = 'Luminal',
  '1' = 'Luminal',
  '2' = 'Luminal',
  '3' = 'Luminal',
  '4' = 'Endothelial',
  '5' = 'T cell',
  '6' = 'Luminal',
  '7' = 'Monocytic',
  '8' = 'Fibroblast',
  '9' = 'Mast',
  '10' = 'Basal/intermediate',
  '11' = 'Unknown',
  '12' = 'Unknown'
)

P = DotPlot(chen,features = marker_genes_chen, dot.scale = 8) + RotatedAxis() + labs(y="Cell cluster",x="Genes")
ggsave("dotplot_chen_etal_annotated_clusters.pdf",P,device="pdf",width=15,height = 4)

P = DimPlot(object = chen, reduction = "tsne")
ggsave("tsneplot_chen_etal_annotated_clusters.pdf",P,device="pdf",width=15,height = 15)

# Unknowns removed
chen = subset(chen,cells = WhichCells(chen,idents = "Unknown",invert = T))

P = DotPlot(chen,features = marker_genes_chen, dot.scale = 8) + RotatedAxis() + labs(y="Cell cluster",x="Genes")
ggsave("dotplot_chen_etal_annotated_clusters_final.pdf",P,device="pdf",width=15,height = 4)

P = DimPlot(object = chen, reduction = "tsne")
ggsave("tsneplot_chen_etal_annotated_clusters_final.pdf",P,device="pdf",width=15,height = 15)

write.table(Idents(chen),file = "chen-et-al-revised-clusters.txt")
