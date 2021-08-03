# Cell segmentation matching
# Antti Kiviaho 3.8.2021
#
#
# The purpose of this script is to match the pixel-coordinates of a Visium experiment
# and a cell-location matrix acquired through nuclei-segmentation. The result
# is a cell-count-by-barcode matrix, which can be exported and integrated with
# gene-expression data.

setwd("~/Dropbox (Compbio)/prostate_spatial/results/segmentation/")
sample_names <- list("BPH_651","BPH_688","CRPC_278","CRPC_489","CRPC_697","PC_03_6712","PC_15420OIK","PC_7875OIK" )#,"PC_4980")
segmentation_path <- paste0("~/Dropbox (Compbio)/prostate_spatial/results/segmentation/run_2/",sample_names,"-seg-results/",sample_names,"-cell-locations.csv")
spot_location_path <- paste0("~/Dropbox (Compbio)/prostate_spatial/results/spaceranger-outs/",sample_names,"/spatial/tissue_positions_list.csv")
spot_radius <- 65

create.n.of.cell.by.bc.matrix <- function(spot_location_df, cell_location_df, spot_size = spot_radius)
{
  result_df <- data.frame(matrix(ncol = 2,nrow = 0))
  colnames(result_df) <- c("barcode","n_of_cells")
  result_df$n_of_cells <- as.numeric(result_df$n_of_cells)
  for(bc_idx in 1:nrow(spot_location_df))
  {
    result_df[bc_idx,] <- c(spot_location_df[bc_idx,"barcode"],0)
    bc_x <- spot_location_df[bc_idx,"x"]
    bc_y <- spot_location_df[bc_idx,"y"]
    for(seg_idx in 1:nrow(cell_location_df))
    {
      cell_x <- cell_location_df[seg_idx,"x"]
      cell_y <- cell_location_df[seg_idx,"y"]
      if(sqrt((bc_x-cell_x)^2+(bc_y-cell_y)^2) < spot_size)
      {
        result_df[bc_idx,"n_of_cells"] <- as.numeric(result_df[bc_idx,"n_of_cells"]) + 1
      }
    }
  }
  result_df$n_of_cells <- as.numeric(result_df$n_of_cells)
  return(result_df)
}

run.mapping.in.pieces <- function(loc_matrix, seg_matrix, nrows = 200)
{
  loc_matrix_split <- split(loc_matrix[order(loc_matrix$y),], (seq(nrow(loc_matrix))-1) %/% nrows)
  result_df <- data.frame(matrix(ncol = 2,nrow = 0))
  colnames(result_df) <- c("barcode","n_of_cells")
  for(idx in 1:length(loc_matrix_split))
  {
    loc_matrix_part <- loc_matrix_split[[idx]]
    seg_matrix_part <- subset(seg_matrix, y >= min(loc_matrix_part$y) - spot_radius & y <= max(loc_matrix_part$y) + spot_radius)
    result_df <- rbind(result_df,create.n.of.cell.by.bc.matrix(loc_matrix_part,seg_matrix_part))
  }
  return(result_df)
}

cell_location_matrices <- lapply(segmentation_path, read.csv, header = T)
spot_location_matrices <- lapply(spot_location_path,read.csv, header = F,col.names = c("barcode","in_tissue","array_row","array_col","y","x"))
spot_location_matrices <- lapply(spot_location_matrices, subset,in_tissue == 1)

result_matrices <- list()

ptm <- proc.time()
for(idx in 1:length(cell_location_matrices))
{
  cell_location_matrix <- cell_location_matrices[[idx]]
  spot_location_matrix <- spot_location_matrices[[idx]]
  result_matrices[[idx]] <- run.mapping.in.pieces(spot_location_matrix, cell_location_matrix)
}
proc.time() - ptm

for(idx in 1:length(result_matrices))
{
  write.csv(result_matrices[[idx]],file=paste0(sample_names[idx],"-cells-per-spot.csv"))
}
