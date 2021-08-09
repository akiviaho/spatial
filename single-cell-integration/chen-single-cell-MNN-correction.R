# 3.8.2021 Antti Kiviaho
#
# Script 1/3 in cell-type identification series of scripts
#
# This script is meant for Mutual Nearest Neighbor (MNN) normalization of data published in
# "Single-cell analysis reveals transcriptomic remodellings in distinct cell types that 
# contribute to human prostate cancer progression" by Chen et al. in Nature: 
# https://www.nature.com/articles/s41556-020-00613-6. 
# MNN-normalized-matrix is outputted in plain txt in order to make it easily accessible
# Further analysis is done with script 2/3, chen-single-cell-clustering
#
# Normalized data is available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4203181
# The script follows the methods used in the original article as closely as possible.
# Uses R 4.0.4

library(batchelor)

setwd("~/Dropbox (Compbio)/prostate_spatial/results/single-cell-integration/")

mycount <- read.table(file = "~/Dropbox (Compbio)/prostate_spatial/data/single-cell-integration/GSM4203181_data.matrix.txt",header = TRUE)
mycount <- SingleCellExperiment(assays=list(logcounts = as.matrix(mycount)));


### block batch effect
original <- list(logcounts(mycount)[, endsWith(colnames(mycount),'.1')],
                 logcounts(mycount)[, endsWith(colnames(mycount),'.2')],
                 logcounts(mycount)[, endsWith(colnames(mycount),'.3')],
                 logcounts(mycount)[, endsWith(colnames(mycount),'.4')],
                 logcounts(mycount)[, endsWith(colnames(mycount),'.5')],
                 logcounts(mycount)[, endsWith(colnames(mycount),'.6')],
                 logcounts(mycount)[, endsWith(colnames(mycount),'.7')],
                 logcounts(mycount)[, endsWith(colnames(mycount),'.8')],
                 logcounts(mycount)[, endsWith(colnames(mycount),'.9')],
                 logcounts(mycount)[, endsWith(colnames(mycount),'.10')],
                 logcounts(mycount)[, endsWith(colnames(mycount),'.11')],
                 logcounts(mycount)[, endsWith(colnames(mycount),'.12')],
                 logcounts(mycount)[, endsWith(colnames(mycount),'.13')])

names(original) <- 1:13

####run fastMNN
out <- do.call(fastMNN, c(original, list(k = 5, d = 50)));



normalized <- data.frame(out@assays@data$reconstructed)
write.table(normalized,file = "MNN-corrected-counts.txt")

# Good practice to save a smaller size file, since filesize is very large
write.table(round(normalized,10),file = "rounded-MNN-corrected-counts.txt")
