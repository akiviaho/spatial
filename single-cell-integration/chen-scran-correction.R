if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.8")

BiocManager::install("scran")
BiocManager::install("scater")
install.packages("remotes")

library(remotes)
install_version("Seurat","2.3.4")


library(scran);
library(Seurat);
library(scater);

setwd("~/Dropbox (Compbio)/prostate_spatial/results/single-cell-integration/")
mycount <- read.table(file = "~/Dropbox (Compbio)/prostate_spatial/data/single-cell-integration/GSM4203181_data.matrix.txt",header = TRUE)

mycount <- SingleCellExperiment(assays=list(logcounts = as.matrix(mycount)));


### block batch effect

count1 <- logcounts(mycount)[, grepl('1$',colnames(mycount))]
count2 <- logcounts(mycount)[, grepl('2',colnames(mycount))];
count3 <- logcounts(mycount)[, grepl('3',colnames(mycount))];
count4 <- logcounts(mycount)[, grepl('4',colnames(mycount))];
count5 <- logcounts(mycount)[, grepl('5',colnames(mycount))];
count6 <- logcounts(mycount)[, grepl('6',colnames(mycount))];
count7 <- logcounts(mycount)[, grepl('7',colnames(mycount))];
count8 <- logcounts(mycount)[, grepl('8',colnames(mycount))];
count9 <- logcounts(mycount)[, grepl('9',colnames(mycount))];
count10 <- logcounts(mycount)[, grepl('10',colnames(mycount))];
count11 <- logcounts(mycount)[, grepl('11',colnames(mycount))];
count12 <- logcounts(mycount)[, grepl('12',colnames(mycount))];
count13 <- logcounts(mycount)[, grepl('13',colnames(mycount))];

####
original <- list(count1, count2, count3, count4, count5, count6, count7, count8,
                 count9, count10, count11, count12, count13);
names(original) <- 1:13
####run fastMNN
out <- do.call(fastMNN, c(original, list(k = 20, d = 50, approximate = TRUE)));
####
combined <- do.call(cbind, original);
sce <- SingleCellExperiment(list(logcounts = combined));
reducedDim(sce, 'MNN') <- out$corrected;
#sce$Batch <- paste0('p', batch[order(batch)]);
sce$Batch <- epi.ne@meta.data[rownames(sce@colData), ]$orig.ident;

####run pca
set.seed(100);
osce <- runPCA(sce, ntop=Inf, method="irlba");
osce <- runTSNE(osce, use_dimred="PCA")
set.seed(100)
csce <- runTSNE(sce, use_dimred="MNN");

####
ot <- plotTSNE(osce, colour_by="Batch") + ggtitle("Original") 
ct <- plotTSNE(csce, colour_by="Batch") + ggtitle("Corrected")
pdf(generate.filename('plottsne', paste0('batch_compare_', name), 'pdf'), width = 15)
multiplot(ot, ct, cols=2);
dev.off();
###
csce$group <- epi.ne@meta.data[rownames(csce@colData), ]$celltype;
pdf(generate.filename('plottsne', paste0('corrected_group_', name), 'pdf'), width = 8)
plotTSNE(csce, colour_by="group") + ggtitle("Corrected")
dev.off()
