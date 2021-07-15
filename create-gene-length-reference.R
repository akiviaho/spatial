#This script takes a gtf file as input and creates a ENSEMBL-SYMBOL-LENGTH-table.
# AK 14.7.2021

library(GenomicFeatures)
library(biomaRt)
library(tibble)

#Variables
gtf_file_path <- "~/Dropbox (Compbio)/prostate_spatial/data/reference/genes.gtf"
export_file_path <- "~/Dropbox (Compbio)/prostate_spatial/data/reference/gene_lengths.txt"

#Functions
create.gene.length.reference <- function(file = gtf_file_path)
{
  txdb <- makeTxDbFromGFF(file,format="gtf")
  exons.list.per.gene <- exonsBy(txdb,by="gene")
  df <- as.data.frame(sum(width(reduce(exons.list.per.gene))))
  df <- rownames_to_column(df)
  colnames(df) <- c("ENSEMBL","LENGTH")
  return(df)
}

get.symbols <- function(df = gene_lengths)
{
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <- df$ENSEMBL
  G_list <- getBM(filters = "ensembl_gene_id", 
                  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
                  values = genes, mart = mart)
  G_list <- G_list[,1:2]
  colnames(G_list) <- c("ENSEMBL","SYMBOL")
  df_result <- merge(x = G_list ,y = df ,by="ENSEMBL", all =FALSE)
  
  return(df_result)
}

#Execute 
gene_lengths <- create.gene.length.reference()
ensembl_symbol_lengths <- get.symbols()

write.csv(ensembl_symbol_lengths,file = export_file_path)