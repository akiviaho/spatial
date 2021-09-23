# Antti Kiviaho 19.8.2021
#
# Script for running inferCNV with ready-subsetted data
#
#

library(infercnv)
library(data.table)

setwd("~/Dropbox (Compbio)/prostate_spatial/data/inferCNV-subsets/inferCNV_analysis_data_0.3_epithelial_all_samples/")
results_dir = "~/Dropbox (Compbio)/prostate_spatial/results/secondary-analysis-in-r/inferCNV_analysis_results_9585_spots_all_samples_ref_group_BPH_cutoff_0.01"
counts_file = "raw_counts_9585_spots.txt"
annotations = "spot_annotations.txt"
gene_locations = "gene_locations.txt"
ref_group_id = "BPH"
own_cutoff = 0.01

sample_names <- list("BPH_651","BPH_688","CRPC_278","CRPC_489","CRPC_697","PC_03_6712","PC_15420OIK","PC_7875OIK","PC_4980")

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_file,
                                    annotations_file=annotations,
                                    gene_order_file=gene_locations,
                                    ref_group_names=unlist(sample_names[which(sample_names %like% ref_group_id)]))
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=own_cutoff,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=results_dir,  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)

