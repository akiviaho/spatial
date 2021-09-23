# Antti Kiviaho 17.8.2021
#
# The purpose of this script is to explore the distribution of each cell type
# and find the most representative spots for them all.

library(ggplot2)
library(data.table)

setwd("~/Dropbox (Compbio)/prostate_spatial/results/secondary-analysis-in-r/")
sample_names <- list("BPH_651","BPH_688","CRPC_278","CRPC_489","CRPC_697","PC_03_6712","PC_15420OIK","PC_7875OIK","PC_4980")
sample_types <- list("BPH","BPH", "CRPC","CRPC","CRPC","PC","PC","PC","PC")

# Cell-type proportions matrix ####
data <- read.table("~/Dropbox (Compbio)/prostate_spatial/data/single-cell-mapping/W_mRNA_count_q05.csv",sep = ",",header = T)

cell_type_proportions <- data.frame(matrix(nrow=nrow(data),ncol=ncol(data)))
colnames(cell_type_proportions) <- gsub("q05_nUMI_factors",replacement = "",colnames(data))
cell_type_proportions$spot_id <- data$spot_id

for(idx in 2:ncol(data))
{
  cell_type_proportions[,idx] <- data[,idx]/rowSums(data[,-1])
}
cell_type_proportions$Epithelial <- cell_type_proportions$Basal.intermediate + cell_type_proportions$Luminal

# Functions ####
plot.histogram <- function(proportions,type,y.max)
{
  p1 <- ggplot(data = proportions, mapping = aes_(x = as.name(type))) +
    geom_histogram(binwidth = 0.04,stat = "count") +
    scale_x_binned(n.breaks = 10,limits = c(0,1)) +
    scale_y_continuous(limits = c(0,y.max)) +
    labs(title = paste0(type,"-cell type proportion distribution") 
         ,x = paste0("Cell type proportion - ",type)) +
    theme_bw()
  plot(p1)
}

plot.epithelial.density <- function(data)
{
  tmp <- data[,c("spot_id","Epithelial")]
  tmp$sample <- NA
  for(sample_idx in 1:length(sample_names))
  {
    tmp[tmp$spot_id %like% sample_names[[sample_idx]],"sample"] <- sample_names[[sample_idx]]
  }
  
  ggplot(data = tmp, aes(x=Epithelial,fill=sample)) +
    geom_density(alpha = 0.5) +
    theme_bw() +
    labs(title = "Proportion density for Epithelial-cell type", x = "Proportion of epithelial cells")
}
# Plotting ####

# Plot Epithelial cell-type density
pdf(file="Epithelial_proportion_density_sample_wise.pdf")
plot.epithelial.density(cell_type_proportions)
dev.off()


# Distributions of cell-proportions for the whole dataset
pdf(file="all-cell-type-proportion-distributions.pdf")
for(idx in 2:ncol(cell_type_proportions))
{
  plot.histogram(cell_type_proportions,colnames(cell_type_proportions)[idx], nrow(cell_type_proportions))
}
dev.off()


# Distributions of cell-proportions for each independent sample
for(sample_idx in 1:length(sample_names))
{
  sample <- cell_type_proportions[cell_type_proportions$spot_id %like% sample_names[[sample_idx]],]
  pdf(file=paste0(sample_names[[sample_idx]],"-cell-type-proportion-distributions.pdf"))
  for(idx in 2:ncol(sample))
  {
    plot.histogram(sample,colnames(sample)[idx],nrow(sample))
  }
  dev.off()
}


# Subsetting test ####
epithelial_spots <- subset(cell_type_proportions,Epithelial > 0.5)

epithelial_share_per_sample <- data.frame(matrix(ncol=2,nrow=length(sample_names)))
colnames(epithelial_share_per_sample) <- c("sample","epithelial_per")
epithelial_share_per_sample$sample <- as.vector(sample_names)
  
for(idx in 1:length(sample_names))
{
  epithelial_spots_per_sample <- nrow(epithelial_spots[epithelial_spots$spot_id %like% sample_names[[idx]],])
  total_spots_per_sample <- nrow(cell_type_proportions[cell_type_proportions$spot_id %like% sample_names[[idx]],])
  epithelial_share_per_sample[idx,2] <- round(epithelial_spots_per_sample/total_spots_per_sample*100,2)
}
View(epithelial_share_per_sample)


