
library(batchelor)


setwd("~/Dropbox (Compbio)/prostate_spatial/results/single-cell-integration/")
mycount <- read.table(file = "~/Dropbox (Compbio)/prostate_spatial/data/single-cell-integration/GSM4203181_data.lite.version.txt",header = TRUE)

mycount <- SingleCellExperiment(assays=list(logcounts = as.matrix(mycount)));


### block batch effect
original <- list(logcounts(mycount)[, grepl('1$',colnames(mycount))],
logcounts(mycount)[, grepl('2',colnames(mycount))],
logcounts(mycount)[, grepl('3',colnames(mycount))],
logcounts(mycount)[, grepl('4',colnames(mycount))],
logcounts(mycount)[, grepl('5',colnames(mycount))],
logcounts(mycount)[, grepl('6',colnames(mycount))],
logcounts(mycount)[, grepl('7',colnames(mycount))],
logcounts(mycount)[, grepl('8',colnames(mycount))],
logcounts(mycount)[, grepl('9',colnames(mycount))],
logcounts(mycount)[, grepl('10',colnames(mycount))],
logcounts(mycount)[, grepl('11',colnames(mycount))],
logcounts(mycount)[, grepl('12',colnames(mycount))],
logcounts(mycount)[, grepl('13',colnames(mycount))])


names(original) <- 1:13
####run fastMNN
out <- do.call(fastMNN, c(original, list(k = 20, d = 50)));
####
# combined <- do.call(cbind, original);
# sce <- SingleCellExperiment(list(logcounts = combined));
# reducedDim(sce, 'MNN') <- out$corrected;
# sce$Batch <- out$batch

saveRDS(out, file = "MNN-out.rds");
write.table(as.data.frame(out@assays@data$reconstructed),file = "MNN-corrected-counts.txt");
# saveRDS(combined, file ="combined.rds");
# saveRDS(sce, file ="sce.rds");

####run pca
# set.seed(100);
# osce <- runPCA(sce, ntop=Inf, method="irlba");
# osce <- runTSNE(osce, use_dimred="PCA")
# set.seed(100)
# csce <- runTSNE(sce, use_dimred="MNN");
# 
# ####
# ot <- plotTSNE(osce, colour_by="Batch") + ggtitle("Original") 
# ct <- plotTSNE(csce, colour_by="Batch") + ggtitle("Corrected")
# pdf(generate.filename('plottsne', paste0('batch_compare_', name), 'pdf'), width = 15)
# multiplot(ot, ct, cols=2);
# dev.off();
###
# csce$group <- epi.ne@meta.data[rownames(csce@colData), ]$celltype;
# pdf(generate.filename('plottsne', paste0('corrected_group_', name), 'pdf'), width = 8)
# plotTSNE(csce, colour_by="group") + ggtitle("Corrected")
# dev.off()
