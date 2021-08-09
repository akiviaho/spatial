# Antti Kiviaho 4.8.2021
#
# The purpose of this script is to explore the cell-count and gene-expression
# relations in order to develop effective normalization methods for the data.
# 
# Correlation of cell-count, gene count and UMI count are calculated per sample
# Correlation of cell-count and per gene count are calculated per sample for all 
# genes.


library(Seurat)
library(hdf5r)


setwd("~/Dropbox (Compbio)/prostate_spatial/results/secondary-analysis-in-r/cell-count-analysis/")
sample_names <- list("BPH_651","BPH_688","CRPC_278","CRPC_489","CRPC_697","PC_03_6712","PC_15420OIK","PC_7875OIK","PC_4980")
sample_types <- list("BPH","BPH", "CRPC","CRPC","CRPC","PC","PC","PC","PC")
samples <- readRDS("~/Dropbox (Compbio)/prostate_spatial/data/samples-seurat-object.rds")


plot.correlation <- function(smpls)
{
  for(idx in 1:length(smpls))
  {
      data <- data.frame(cbind(smpls[[idx]]$nFeature_Spatial, smpls[[idx]]$nCount_Spatial, smpls[[idx]]$nCells))
      colnames(data) <- c("nFeature_Spatial","nCount_Spatial","nCells")
      cor_nCount <- cor.test(data$nCells, data$nCount_Spatial)
      cor_nFeature <- cor.test(data$nCells, data$nFeature_Spatial)
    p1 <- ggplot(data, mapping = aes(x = nCells, y = nCount_Spatial)) +
            geom_smooth(method = "gam") +
            geom_point() +
            ggtitle(names(smpls[idx])) +
            annotate(x=max(data$nCells)-4, y=max(data$nCount_Spatial), 
                     label=paste("R =", round(cor_nCount$estimate,2),"\n",
                                 "p =", signif(cor_nCount$p.value,2)), 
                     geom="text", size=5)
            
    p2 <- ggplot(data, mapping = aes(x = nCells, y = nFeature_Spatial)) +
            geom_smooth(method = "gam") +
            geom_point() +
            ggtitle(names(smpls[idx])) +
            annotate(x=max(data$nCells)-4, y=max(data$nFeature_Spatial),
                     label=paste("R =", round(cor_nFeature$estimate,2),"\n",
                                 "p =", signif(cor_nFeature$p.value,2)), 
                     geom="text", size=5)
          
    plot(p1)
    plot(p2)
  }
}

test.nCount.nFeature.correlation <- function(smpls)
{
  df <- data.frame(matrix(nrow = 4,ncol = 9))
  rownames(df) <- c("nCount_nCells_cor","nCount_nCells_pval","nFeature_nCells_cor","nFeature_nCells_pval")
  colnames(df) <- names(samples)
  for(idx in 1:length(smpls))
  {
    data <- data.frame(cbind(smpls[[idx]]$nFeature_Spatial, smpls[[idx]]$nCount_Spatial, smpls[[idx]]$nCells))
    colnames(data) <- c("nFeature_Spatial","nCount_Spatial","nCells")
    df[1,idx] <- cor.test(data$nCount_Spatial,data$nCells)$estimate
    df[2,idx] <- cor.test(data$nCount_Spatial,data$nCells)$p.value
    df[3,idx] <- cor.test(data$nFeature_Spatial,data$nCells)$estimate
    df[4,idx] <- cor.test(data$nFeature_Spatial,data$nCells)$p.value
  }
  return(df)
}

test.genewise.gene.count.correlation <- function(smpls)
{
  df <- data.frame(matrix(nrow = nrow(smpls[[1]]), ncol = 9))
  rownames(df) <- rownames(smpls[[1]])
  colnames(df) <- names(smpls)
  for(idx in 1:length(smpls))
  {
    data <- data.frame(cbind(smpls[[idx]]$nCells))
    genecounts <- t(as.matrix(GetAssayData(smpls[[idx]])))
    data <- data[order(match(colnames(data),colnames(genecounts)))]
    colnames(data) <- c("nCells")
    data <- cbind(data,genecounts)
    for(idx2 in 2:ncol(data))
    {
      df[idx2-1,idx] <- cor.test(data$nCells,data[,idx2])$p.value
    }
  }
  adjusted_df <- as.data.frame(matrix(p.adjust(as.matrix(df),"fdr"),ncol = ncol(df)))
  adjusted_df <- signif(adjusted_df)
  rownames(adjusted_df) <- rownames(df)
  colnames(adjusted_df) <- colnames(df)
  adjusted_df <- signif(adjusted_df,3)
  return(adjusted_df)
}

get.signif.genes <- function(pval_df, thr = 0.05)
{
  df <- data.frame(matrix(ncol = ncol(pval_df),nrow = 0))
  colnames(df) <- colnames(pval_df)
  for(idx in 1:nrow(pval_df))
  {
    if(!is.na(all(pval_df[idx,] < thr)))
    {
      if(all(pval_df[idx,] < thr))
      {
        df <- rbind(df,pval_df[idx,])
      }
    }
  }
  return(df)
}

pdf(file = "correlation-plots.pdf")
plot.correlation(samples)
dev.off()

res <- test.nCount.nFeature.correlation(samples)
res <- signif(res, digits=2)
res[c(1,3),] <- format(res[c(1,3),],scientific=F)
write.table(res,file = "correlation-test-samplewise.txt")


adj_genewise_corr <- test.genewise.gene.count.correlation(samples)
signif_genes_corr <- get.signif.genes(adj_genewise_corr[,-which(colnames(adj_genewise_corr) 
                                                                %in% c("BPH_688","CRPC_278","PC_03_6712"))])
write.table(signif_genes_corr,"correlation-test-genewise.txt")
# Cell count didn't correlate for CRPC_278 and PC_03_6712, segmentation was unsuccesful for BPH_688

