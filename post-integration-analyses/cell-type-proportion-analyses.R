# Antti Kiviaho 18.8.2021
#
# 
# 
library(Seurat)
library(patchwork)
library(data.table)

setwd("~/Dropbox (Compbio)/prostate_spatial/results/secondary-analysis-in-r/") # Where to save plots
path_for_data <- "~/Dropbox (Compbio)/prostate_spatial/data/inferCNV-subsets/" # where to save data created for inferCNV
sample_names <- list("BPH_651","BPH_688","CRPC_278","CRPC_489","CRPC_697","PC_03_6712","PC_15420OIK","PC_7875OIK","PC_4980")
# sample_quantiles <- list(0.38,0.36,1.00,0.37,0.50,0.81,0.49,0.46,0.51)

# Use this if the goal is to just find the most representative spots in BPH samples and get all spots for cancer samples.
sample_quantiles <- list(0.38,0.36,1,1,1,1,1,1,1) 

USE_QUANTILES = TRUE
per_threshold <- 0.50
qnt_threshold <- 0.70 # What quantile of the purity metric should be included per sample? This is a hard-coded value not accounting for
# differences in tissue composition. Some samples might have a larger share of epithelial tissue and thus require different 
# quantiles for thresholding. Alternative method is a vector of sample-specific thresholds evaluated by examining the share of 
# epithelial tissue.

###############

# Create Seurat object
data <- readRDS(file="~/Dropbox (Compbio)/prostate_spatial/data/samples-seurat-object.rds")
data <- merge(data[[1]], y = unlist(data)[-1], add.cell.ids = unlist(sample_names), project = "aggregated")
names(data@images) <- unlist(sample_names)

# Add Cell-type proportions to Seurat object
proportions <- read.table("~/Dropbox (Compbio)/prostate_spatial/data/single-cell-mapping/W_mRNA_count_q05.csv",
                          sep = ",",header = T) # Read in the cell-associated mRNA-counts

# Modify matrix into a spot-proportion matrix
cell_type_proportions <- data.frame(matrix(nrow=nrow(proportions),ncol=ncol(proportions)))
colnames(cell_type_proportions) <- gsub("q05_nUMI_factors",replacement = "",colnames(proportions))
cell_type_proportions$spot_id <- proportions$spot_id
for(idx in 2:ncol(proportions))
{
  cell_type_proportions[,idx] <- proportions[,idx]/rowSums(proportions[,-1])
}

cell_type_proportions$epithelial <- cell_type_proportions$Basal.intermediate + cell_type_proportions$Luminal
cell_type_proportions <- tibble::column_to_rownames(cell_type_proportions,var="spot_id")
cell_type_proportions$representative_epithelial <- NA

for(idx in 1:length(sample_names))
{
  if(USE_QUANTILES)
  {
    threshold <- round(quantile(cell_type_proportions[which(rownames(cell_type_proportions) %like% 
                                                            sample_names[[idx]]),"epithelial"],1-sample_quantiles[[idx]]),2)
  } else { threshold <- per_threshold }
  
  cell_type_proportions[which(rownames(cell_type_proportions) %like% sample_names[[idx]] & 
                                cell_type_proportions$epithelial > threshold),"representative_epithelial"] <- paste0("epithelial") 
}
data <- AddMetaData(data,metadata = cell_type_proportions)

##############

# Plot proportions on the tissue sections to evaluate the accuracy
if(USE_QUANTILES)
{
  pdf(paste0("epithelial_proportion_on_spatial_samplewise_sample_specific_threshold.pdf"))
  
  for(idx in 1:length(sample_names))
  {
    p1 <- SpatialFeaturePlot(data,features = "epithelial", images = sample_names[[idx]], crop = T, alpha = 1) + 
      ggplot2::ggtitle(paste0(sample_names[[idx]],", qnt thr ",sample_quantiles[[idx]]))
    p2 <- SpatialPlot(data, group.by = "representative_epithelial", images = sample_names[[idx]], crop = T, alpha = 1, cols = "blue")
    plot(p1 + p2) 
  }
  
} else { pdf(paste0("epithelial_proportion_on_spatial_samplewise_sample_specific_threshold.pdf")) }

for(idx in 1:length(sample_names))
{
  p1 <- SpatialFeaturePlot(data,features = "epithelial", images = sample_names[[idx]], crop = T, alpha = 1) + 
    ggplot2::ggtitle(paste0(sample_names[[idx]],", per threshold ",per_threshold))
  p2 <- SpatialPlot(data, group.by = "representative_epithelial", images = sample_names[[idx]], crop = T, alpha = 1, cols = "blue")
  plot(p1 + p2) 
}

dev.off()


###############

# Subset data to only contain epithelial representative spots

epithelial_data <- subset(data, representative_epithelial %like% "epithelial")
epithelial_data <- data.frame(epithelial_data@assays$Spatial@counts)
epithelial_data <- epithelial_data[-which(rowSums(epithelial_data) == 0),] # Remove genes with no counts


# Construct infercnv reference file
gene_locations <- read.table("~/Dropbox (Compbio)/prostate_spatial/data/reference/all_genes_locations_for_icnv.txt")
colnames(gene_locations) <- c("chr","start","end","SYMBOL")
gene_locations <- gene_locations[which(!duplicated(gene_locations$SYMBOL)),]
rownames(gene_locations) <- NULL
gene_locations <- tibble::column_to_rownames(gene_locations,var = "SYMBOL")

# Create annotations
annotations <- data.frame(matrix(nrow = ncol(epithelial_data), ncol = 1))
colnames(annotations) <- c("sample")
rownames(annotations) <- colnames(epithelial_data)
for(idx in 1:length(sample_names))
{
  annotations[which(rownames(annotations) %like% sample_names[[idx]]),"sample"] <- sample_names[[idx]]
}

data_dir <- paste0(path_for_data,"inferCNV_analysis_data_BPH_specific_threshold_epithelial_all_cancer_sample_spots")
dir.create(data_dir)
setwd(data_dir)

write.table(epithelial_data,file = paste0("raw_counts_",ncol(epithelial_data),"_spots.txt"), sep = "\t",quote = F)
write.table(x = gene_locations,file = "gene_locations.txt",col.names = F,sep = "\t", quote =F)
write.table(annotations, file = "spot_annotations.txt",col.names = F, sep = "\t", quote = F)
