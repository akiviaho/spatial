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
create.nUMI.data.frame <- function(object_list, lower_quantile, upper_quantile, names = sample_names, use.genes =F, max_n_spots = MAX_N_SPOTS)
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
  df <- gather(df, sample, count, names[[1]]:names[[length(names)]],factor_key = T)
  df <- na.omit(df)
  df <- df[which(df$count < quantile(df$count,upper_quantile) & df$count
                 > quantile(df$count,lower_quantile)),]
  return(df)
}

plot.nUMI.dotplot <- function(data, lower_quantile = 0, upper_quantile = 1, names = sample_names, use.genes = F)
{
  data <- create.nUMI.data.frame(data, lower_quantile, upper_quantile, names, use.genes)
  
  if(use.genes)
  {
    ggplot(data, mapping = aes(data, x=x, y=count, col = sample)) +
      geom_point() +
      labs(title = paste0("Filtered with: lower quantile: ",lower_quantile," upper quantile: ", upper_quantile), y = "Gene Count", x = "Spots") +
      scale_y_log10() +
      scale_x_continuous(breaks=seq(0,max(data$x),200))
  }
  else
  {
    ggplot(data, mapping = aes(data, x=x, y=count, col = sample)) +
      geom_point() +
      labs(title = paste0("Filtered with: lower quantile: ",lower_quantile," upper quantile: ", upper_quantile), y = "UMI Count", x = "Spots") +
      scale_y_log10() +
      scale_x_continuous(breaks=seq(0,max(data$x),200))
  }
}

create.saturation.data.frame <- function(object_list, lower_quantile, upper_quantile, names = sample_names)
{
  df <- matrix(nrow = 0,ncol = 4)
  for(idx in 1:length(object_list))
  {
    UMI_count <- c(object_list[[idx]]$nCount_RNA)
    Gene_count <- c(object_list[[idx]]$nFeature_RNA)
    sample <- c(rep(names[[idx]],length(UMI_count)))
    x <- c(rep(NA,length(UMI_count)))
    sample.df <- data.frame(cbind(x,sample,UMI_count,Gene_count))
    sample.df$UMI_count <- as.integer(sample.df$UMI_count)
    sample.df$Gene_count <- as.integer(sample.df$Gene_count)
    sample.df <- sample.df[order(sample.df$UMI_count,decreasing = F),]
    df <- rbind(df,sample.df)
  }
  rownames(df) <- NULL
  df <- data.frame(df)
  df$x <- c(1:nrow(df))
  df$x <- factor(df$x)
  df <- df[which(df$UMI_count < quantile(df$UMI_count,upper_quantile) & df$UMI_count
                 > quantile(df$UMI_count,lower_quantile)),]
  return(df)
}

plot.saturation <- function(data, lower_quantile = 0, upper_quantile = 1, names = sample_names)
{
  data <- create.saturation.data.frame(data, lower_quantile, upper_quantile)
  
  ggplot(data, mapping = aes(data, x=UMI_count, y=Gene_count, col = sample)) +
    geom_smooth() +
    geom_point(size = 0.2, alpha = 0.5) +
    labs(title = paste0("Filtered with: lower quantile: ",lower_quantile," upper quantile: ", upper_quantile),y = "Gene count", x = "UMI count") +
    scale_x_continuous(breaks=seq(0,max(data$UMI_count),5000))
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

n_of_spots <-  lapply(samples,function(samples)
                      {
                      n_of_spots <- WhichCells(samples)
                      n_of_spots <- length(n_of_spots)
                      n_of_spots
                      })
MAX_N_SPOTS <- max(unlist(n_of_spots))

# Generate plots

pdf(file = "UMI-count-plots.pdf")
plot.nUMI.dotplot(samples)
VlnPlot(samples_aggregated,features = "nCount_RNA", log = T)
plot.nUMI.dotplot(samples, lower_quantile = 0.05, upper_quantile = 0.95)
VlnPlot(subset(samples_aggregated, subset = 
                 nCount_RNA > quantile(samples_aggregated$nCount_RNA, 0.05) & 
                 nCount_RNA < quantile(samples_aggregated$nCount_RNA, 0.95)),features = "nCount_RNA", log = T)
dev.off()


pdf(file = "gene-count-plots.pdf")
plot.nUMI.dotplot(samples,use.genes=T)
VlnPlot(samples_aggregated,features = "nFeature_RNA", log = T)
plot.nUMI.dotplot(samples, lower_quantile = 0.05, upper_quantile = 0.95, use.genes = T)
VlnPlot(subset(samples_aggregated, subset = 
                 nFeature_RNA > quantile(samples_aggregated$nFeature_RNA, 0.05) & 
                 nFeature_RNA < quantile(samples_aggregated$nFeature_RNA, 0.95)),features = "nFeature_RNA", log = T)
dev.off()


pdf(file = "saturation-plot.pdf")
plot.saturation(samples)
plot.saturation(samples,lower_quantile = 0.01, upper_quantile = 0.99)
plot.saturation(samples,lower_quantile = 0.05, upper_quantile = 0.95)
dev.off()
