# Antti Kiviaho 4.8.2021
#
# This script is for the integration of cell-count data inferred through nuclei-
# segmentation and barcode mapping. Saves a Seurat object with necessary data


library(Seurat)
library(hdf5r)

# SETUP  ####

# Download data
data_folder_path <- "~/Dropbox (Compbio)/prostate_spatial/results/spaceranger-outs"
sample_names <- list("BPH_651","BPH_688","CRPC_278","CRPC_489","CRPC_697","PC_03_6712","PC_15420OIK","PC_7875OIK","PC_4980")
samples <- lapply(file.path(data_folder_path,sample_names),Load10X_Spatial)

## Download cell counts ###

cells_per_spot <- lapply(paste0("~/Dropbox (Compbio)/prostate_spatial/results/segmentation/",sample_names,"-cells-per-spot.csv"),read.csv)
names(cells_per_spot) <- sample_names

# Integrate cell counts into Seurat object
for(idx in 1:length(samples))
{
  sample_cell_count <- as.vector(cells_per_spot[[idx]]$n_of_cells)
  names(sample_cell_count) <- cells_per_spot[[idx]]$barcode
  samples[[idx]]$nCells <- sample_cell_count
}
names(samples) <- sample_names

saveRDS(samples,file = "~/Dropbox (Compbio)/prostate_spatial/data/samples-seurat-object.rds")