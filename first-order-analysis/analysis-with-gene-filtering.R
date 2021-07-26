# Antti Kiviaho 8.7.2021
# This script is meant for producing plots of a Spatial dataset filtered with
# a geneset of interest. The genelist must be by in SYMBOl format.
# Outputted plots are:
# 1. PCA-plot with 2 PC's
# 2. Similarity matrix with all PC's

library(Seurat)
library(hdf5r)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)



# SETUP ####
# Variables
data_folder_path <- "~/Dropbox (Compbio)/prostate_spatial/results/spaceranger-outs"
results_folder_path <- "~/Dropbox (Compbio)/prostate_spatial/results/secondary-analysis-in-r"
gene_lengths_file<- "~/Dropbox (Compbio)/prostate_spatial/data/reference/gene_lengths.txt"
setwd(results_folder_path)

sample_names <- list("BPH_651","BPH_688","CRPC_278","CRPC_489","CRPC_697","PC_03_6712","PC_15420OIK","PC_7875OIK","PC_4980")
sample_types <- list("BPH","BPH", "CRPC","CRPC","CRPC","PC","PC","PC","PC")

genes_of_interest_file <- "DE-genes-per-cluster-names.csv"
genelist_name <- "DE-genes from Seurat clusters"

col <- magma(50,direction = -1)
# Functions
create.bulk.sample <- function(data, names = sample_names)
{
  result_df <- matrix(ncol = 0,nrow = nrow(data[[1]]@assays$Spatial)) # All samples should contain same n of features
  for(idx in 1:length(data))
  {
    sample_df <- data[[idx]]@assays$Spatial
    bulk_df <- rowSums(sample_df)
    result_df <- cbind(result_df,bulk_df)
  }
  colnames(result_df) <- unlist(names)
  return(data.frame(result_df))
}

TPM.normalize <- function(data, filter_genes = TRUE, gene_lengths = gene_lengths_file, names = sample_names)
{
  data <- create.bulk.sample(data)
  
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
  
  # if (filter_genes)
  # {
  #  merged_df <- merged_df[rownames(merged_df) %in% genelist,]
  # }
  return(merged_df)
  
}

plot.PCA <- function(data, names = sample_names, types = sample_types)
{
  pca <- prcomp(t(data),scale. = FALSE)
  var.percentage <- round(pca$sdev^2/sum(pca$sdev^2)*100,3)
  pca_data <- data.frame(sample=colnames(data),
                         type = unlist(types),
                         X = pca$x[,1],
                         Y = pca$x[,2])
  
  ggplot(data=pca_data, aes(x=X, y=Y, col =sample, shape =  type)) +
    geom_point(size = 3) +
    xlab(paste("PC1 - ", var.percentage[1], "%", sep="")) +
    ylab(paste("PC2 - ", var.percentage[2], "%", sep="")) +
    ggtitle(paste0("Pseudobulk PCA, using genes from ", genelist_name))
}

calculate.similarity.matrix <- function(data, names = sample_names, types = sample_types)
{
  pca <- prcomp(t(data),scale. = FALSE)
  var.percentage <- pca$sdev^2/sum(pca$sdev^2)
  pca_scaled <- sweep(pca$x,2,var.percentage,'*')
  result <- matrix( nrow = nrow(pca_scaled))
  for(idx1 in 1:nrow(pca_scaled))
  {
    temp_vector2 <- c()
    for(idx2 in 1:ncol(pca_scaled))
    {
      temp_vector3 <- c()
      for(idx3 in 1:ncol(pca_scaled))
      {
        temp_vector3[idx3] <- sqrt((pca_scaled[idx1,idx3]-pca_scaled[idx2,idx3])^2)
      }
      temp_vector2[idx2] <- sum(temp_vector3)
    }
    result <- cbind(result,temp_vector2)
  }
 result <- as.data.frame(result[,-1])
 rownames(result) <- unlist(names)
 colnames(result) <- unlist(names)
 return(result)
}


# EXECUTE ####

# Import data and assign names
samples <- lapply(file.path(data_folder_path,sample_names),Load10X_Spatial)

genelist <- read.table(genes_of_interest_file, quote="\"", comment.char="")$V1

for(idx in 1:length(samples)){ Idents(samples[[idx]])<- sample_names[[idx]]} 

pdf(file="gene-filtered-PCA-and-similarity-matrix.pdf")
normalized_bulk_sample <- TPM.normalize(samples)
plot.PCA(normalized_bulk_sample)

result <- calculate.similarity.matrix(normalized_bulk_sample)
pheatmap(-log(result+1),color = col,angle_col = c("315"),legend = F)
dev.off()
