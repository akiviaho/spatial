# Antti Kiviaho 2.9.2021
#
# Ad-hoc script to find CTCF-target genes from https://doi.org/10.1093/nar/gkz462
# These genes are extracted from Figure 4H with multiple evidence speaking in favor.
# Can be used in visualizing CRPC stuff


counts = read.table('~/Downloads/gkz462_supplemental_files (1)/supp_table2.csv',sep=',',header = T)
counts <- counts[!duplicated(counts$Gene),]
rownames(counts) <- NULL
counts = tibble::column_to_rownames(counts,var='Gene')

CTCF_knockout_DE_genes <- c('MYC','ZBTB39','AAAS','ANKRD46','LRP10','CREG1','MFGE8','ZNF580',
                            'BOLA3','FAM214B','SPSB1','RBPMS','ATP8B3','SAAL1','ZMAT3','LPIN1',
                            'RGS16','SLFN5','ATP5G1','ATP2B4','HINT3','FAM102B','RBM45','C19orf60',
                            'IFI6','TRIM14','MRPS2','C1QBP','HEXIM1','NPM3','C5orf56','CIZ1',
                            'DTWD2','BSCL6','DNPH1','HDDC3','PIGO','ALKBH2')
de_genes <- counts[CTCF_knockout_DE_genes,]
de_genes_up <- de_genes[which(de_genes$c3_with.IAA > de_genes$c3_no.IAA & de_genes$c16_with.IAA > de_genes$c16_no.IAA),]
de_genes_down <- de_genes[which(de_genes$c3_with.IAA < de_genes$c3_no.IAA & de_genes$c16_with.IAA < de_genes$c16_no.IAA),]

write.csv(de_genes_up,"ctcf_depletion_upreg_genes.csv")
write.csv(de_genes_down,"ctcf_depletion_downreg_genes.csv")