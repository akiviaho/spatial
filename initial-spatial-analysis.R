# Antti Kiviaho 8.7.2021
# This script is meant for initial analysis of a Spatial Transcriptomics dataset implemented with the Visium technology. Three plots are produced:
# 1. UMI count plot by spot + violin spots for all samples
# 2. Sequencing saturation plot, with UMI-count and Gene-count as axes
# 3. PCA plot with two principal components

library(Seurat)
library(hdf5r)
library(tidyr)
library(ggplot2)
library(patchwork)

# SETUP  ####

#Functions
create.nUMI.data.frame <- function(object_list, use.genes, names = sample_names, max_n_spots = MAX_N_SPOTS)
{
  df <- c(1:max_n_spots)
  
  if(use.genes)
  {
    for(idx in 1:length(object_list))
    {
      sample_nFeature <- c(object_list[[idx]]$nFeature_RNA, rep(NA,max_n_spots-length(object_list[[idx]]$nFeature_RNA)))
      sample_nFeature <- sort(sample_nFeature,decreasing = T,na.last = T)
      df <- cbind(df,sample_nFeature)
    }
  }
  
  else
  {
    for(idx in 1:length(object_list))
    {
      sample_nUMI <- c(object_list[[idx]]$nCount_RNA, rep(NA,max_n_spots-length(object_list[[idx]]$nCount_RNA)))
      sample_nUMI <- sort(sample_nUMI,decreasing = T,na.last = T)
      df <- cbind(df,sample_nUMI)
    }
  }
  rownames(df) <- NULL
  colnames(df) <- c("x",names)
  df <- data.frame(df)
  df$x <- factor(df$x)
  df <- gather(df, sample, count, names[[1]]:names[[length(names)]],factor_key = T)
  df$count <- log10(df$count)
  return(df)
}


plot.nUMI.dotplot <- function(data, names = sample_names, use.genes = F)
{
  data <- create.nUMI.data.frame(data, use.genes)
  
  if(use.genes)
  {
    ggplot(data, mapping = aes(data, x=x, y=count, col = sample)) +
      geom_point() +
      labs(y = "log10 Gene Count", x = "Spots") +
      coord_cartesian(xlim=c(-100,nrow(data)/length(names)+100)) +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
  }
  else
  {
    ggplot(data, mapping = aes(data, x=x, y=count, col = sample)) +
      geom_point() +
      labs(y = "log10 UMI Count", x = "Spots") +
      coord_cartesian(xlim=c(-100,nrow(data)/length(names)+100)) +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
  }
}

# Variables
data_folder_path <- "~/Dropbox (Compbio)/prostate_spatial/results/spaceranger-outs"
results_folder_path <- "~/Dropbox (Compbio)/prostate_spatial/results/secondary-analysis-in-r"
setwd(results_folder_path)

sample_names <- list("BPH_651","BPH_688","CRPC_278","CRPC_489","CRPC_697","PC_03_6712","PC_15420OIK","PC_7875OIK")


# CODE ####

# Import data and assign names
samples <- lapply(file.path(data_folder_path,paste0(sample_names,"_filtered_feature_bc_matrix")),Read10X)
samples <- lapply(samples,CreateSeuratObject)

for(idx in 1:length(samples)){ Idents(samples[[idx]])<- sample_names[[idx]]}

samples_aggregated <- merge(samples[[1]], y = unlist(samples)[-1], add.cell.ids = unlist(sample_names), project = "aggregated")
samples_aggregated$log_nCount_RNA <- log10(samples_aggregated$nCount_RNA)
samples_aggregated$log_nFeature_RNA <- log10(samples_aggregated$nFeature_RNA)

n_of_spots <-  lapply(samples,function(samples)
                      {
                      n_of_spots <- WhichCells(samples)
                      n_of_spots <- length(n_of_spots)
                      n_of_spots
                      })
MAX_N_SPOTS <- max(unlist(n_of_spots))

nUMI_data_frame <- create.nUMI.data.frame(samples)

pdf(file = "UMI-count-plots.pdf")
plot.nUMI.dotplot(samples)
VlnPlot(samples_aggregated,features = "nCount_RNA",y.max = quantile(samples_aggregated$nCount_RNA, .999))
VlnPlot(samples_aggregated,features = "log_nCount_RNA",y.max = quantile(samples_aggregated$log_nCount_RNA, .999))
dev.off()


pdf(file = "gene-count-plots.pdf")
plot.nUMI.dotplot(samples,use.genes=T)
VlnPlot(samples_aggregated,features = "nFeature_RNA",y.max = quantile(samples_aggregated$nFeature_RNA, .999))
VlnPlot(samples_aggregated,features = "log_nFeature_RNA",y.max = quantile(samples_aggregated$log_nFeature_RNA, .999))
dev.off()
