# Antti Kiviaho 8.7.2021
# This script is meant for initial analysis of a Spatial Transcriptomics dataset implemented with the Visium technology. The plots produced are:
# 1. UMI count plot by spot + violin plots
# 2. Gene count plot by spot + violin plots
# 2. Sequencing saturation plot, with UMI-count and Gene-count as axes
# 3. PCA plot with two principal components
# The dataset can be trimmed by setting percentile thresholds for UMI-counts. Plot for both
# the whole data and the trimmed data are produced.

library(Seurat)
library(hdf5r)
library(tidyr)
library(ggplot2)


# SETUP  ####

# Variables
data_folder_path <- "~/Dropbox (Compbio)/prostate_spatial/results/spaceranger-outs"
results_folder_path <- "~/Dropbox (Compbio)/prostate_spatial/results/secondary-analysis-in-r"
gene_lengths_file<- "~/Dropbox (Compbio)/prostate_spatial/data/reference/gene_lengths.txt"
setwd(results_folder_path)

sample_names <- list("BPH_651","BPH_688","CRPC_278","CRPC_489","CRPC_697","PC_03_6712","PC_15420OIK","PC_7875OIK","PC_4980")
sample_types <- list("BPH","BPH", "CRPC","CRPC","CRPC","PC","PC","PC","PC")

lower_percentile_threshold = 0.05
upper_percentile_threshold = 0.95

#Functions
create.nUMI.data.frame <- function(object_list, lower_quantile, upper_quantile, names = sample_names, use.genes =F, max_n_spots = MAX_N_SPOTS)
{
  df <- c(1:max_n_spots)
  
  if(use.genes)
  {
    for(idx in 1:length(object_list))
    {
      sample_nFeature <- c(object_list[[idx]]$nFeature_Spatial, rep(NA,max_n_spots-length(object_list[[idx]]$nFeature_Spatial)))
      sample_nFeature <- sort(sample_nFeature,decreasing = T,na.last = T)
      df <- cbind(df,sample_nFeature)
    }
  }
  
  else
  {
    for(idx in 1:length(object_list))
    {
      sample_nUMI <- c(object_list[[idx]]$nCount_Spatial, rep(NA,max_n_spots-length(object_list[[idx]]$nCount_Spatial)))
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
    UMI_count <- c(object_list[[idx]]$nCount_Spatial)
    Gene_count <- c(object_list[[idx]]$nFeature_Spatial)
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

create.bulk.sample <- function(data, ncount_lower_limit , ncount_upper_limit , names = sample_names)
{
  result_df <- matrix(ncol = 0,nrow = nrow(data[[1]]@assays$Spatial)) # All samples should contain same n of features
  for(idx in 1:length(data))
  {
    sample <- subset(data[[idx]],nCount_Spatial > ncount_lower_limit & nCount_Spatial < ncount_upper_limit)
    sample_df <- sample@assays$Spatial[]
    bulk_df <- rowSums(sample_df)
    result_df <- cbind(result_df,bulk_df)
  }
  colnames(result_df) <- unlist(names)
  return(data.frame(result_df))
}

TPM.normalize <- function(data, gene_lengths = gene_lengths_file, names = sample_names)
{
  
  gene_lengths <- read.csv(gene_lengths, row.names=1)
  gene_lengths <- gene_lengths[,c("SYMBOL","LENGTH")]
  gene_lengths <- gene_lengths[!duplicated(gene_lengths[,c("SYMBOL")]),]
  
  data <- data[-which(rowSums(data) == 0),]
  data <- tibble::rownames_to_column(data,var="SYMBOL")
  
  merged_df <- merge(data,gene_lengths,by="SYMBOL",all=FALSE)
  merged_df <- tibble::column_to_rownames(merged_df,"SYMBOL")
  merged_df <- as.data.frame(sweep(as.matrix(merged_df),1,merged_df$LENGTH/1e+03,'/'))
  merged_df <- merged_df[,unlist(names)]
  
  merged_df <- as.data.frame(sweep(as.matrix(merged_df),2,colSums(merged_df)/1e+06,'/'))
  return(merged_df)

}

plot.PCA <- function(data, ncount_lower_limit = 0, ncount_upper_limit = Inf, names = sample_names, types = sample_types)
{
  bulk_sample <- create.bulk.sample(data, ncount_lower_limit, ncount_upper_limit)
  normalized_bulk_sample <- TPM.normalize(bulk_sample)
  
  pca <- prcomp(t(normalized_bulk_sample),scale. = FALSE)
  var.percentage <- round(pca$sdev^2/sum(pca$sdev^2)*100,2)
  pca_data <- data.frame(sample=colnames(normalized_bulk_sample),
                         type = unlist(types),
                         X = pca$x[,1],
                         Y = pca$x[,2])
  
  ggplot(data=pca_data, aes(x=X, y=Y, col =sample, shape =  type)) +
    geom_point(size = 3) +
    xlab(paste("PC1 - ", var.percentage[1], "%", sep="")) +
    ylab(paste("PC2 - ", var.percentage[2], "%", sep="")) +
    ggtitle(paste0("Pseudobulk PCA, ", ncount_lower_limit," < nUMI < ",ncount_upper_limit))
}

# CODE ####

# Import data and assign names
samples <- lapply(file.path(data_folder_path,sample_names),Load10X_Spatial)

for(idx in 1:length(samples)){ Idents(samples[[idx]])<- sample_names[[idx]]}

samples_aggregated <- merge(samples[[1]], y = unlist(samples)[-1], add.cell.ids = unlist(sample_names), project = "aggregated")
n_lower <- quantile(samples_aggregated$nCount_Spatial,lower_percentile_threshold)
n_upper <- quantile(samples_aggregated$nCount_Spatial,upper_percentile_threshold)


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

VlnPlot(samples_aggregated,features = "nCount_Spatial", log = T) +
  stat_summary(fun = median, geom='point', size = 20, colour = "black",shape=95) +
  theme(legend.position = "none")

plot.nUMI.dotplot(samples, lower_quantile = lower_percentile_threshold, upper_quantile = upper_percentile_threshold)

VlnPlot(subset(samples_aggregated, subset = 
                 nCount_Spatial > quantile(samples_aggregated$nCount_Spatial, lower_percentile_threshold) & 
                 nCount_Spatial < quantile(samples_aggregated$nCount_Spatial, upper_percentile_threshold)),features = "nCount_Spatial", log = T) +
  stat_summary(fun = median, geom='point', size = 20, colour = "black",shape=95) +
  theme(legend.position = "none")
dev.off()

####

pdf(file = "gene-count-plots.pdf")
plot.nUMI.dotplot(samples,use.genes=T)

VlnPlot(samples_aggregated,features = "nFeature_Spatial", log = T,pt.size = 0.005) +
  stat_summary(fun = median, geom='point', size = 20, colour = "black",shape=95) +
  theme(legend.position = "none")

plot.nUMI.dotplot(samples, lower_quantile = lower_percentile_threshold, upper_quantile = upper_percentile_threshold, use.genes = T)

VlnPlot(subset(samples_aggregated, subset = 
                 nCount_Spatial > quantile(samples_aggregated$nCount_Spatial, lower_percentile_threshold) & 
                 nCount_Spatial < quantile(samples_aggregated$nCount_Spatial, upper_percentile_threshold)),features = "nFeature_Spatial", log = T,pt.size = 0.005) +
  stat_summary(fun = median, geom='point', size = 20, colour = "black",shape=95) +
  theme(legend.position = "none")
dev.off()

####

pdf(file = "saturation-plot.pdf")
plot.saturation(samples)
plot.saturation(samples,lower_quantile = lower_percentile_threshold, upper_quantile = upper_percentile_threshold)
dev.off()

####

pdf(file ="PCA-plots.pdf")

plot.PCA(samples)
plot.PCA(samples,ncount_lower_limit = n_lower,ncount_upper_limit = n_upper)

dev.off()
