# Author: Antti Kiviaho
# Date: 8.3.2023
#
# A script for running infercnv on integrated single cell datasets.
# Make sure to identify reference cells via e.g. clustering.
# Also make sure to 'conda activate infercnvpy'
#

#import os
#os.chdir('/lustre/scratch/kiviaho/prostate_spatial')

import scanpy as sc
import infercnvpy as cnv
import anndata as ad
import matplotlib.pyplot as plt

sc.settings.set_figure_params(figsize=(5, 5))

if __name__ == '__main__':

    adata = sc.read_h5ad('./Epithelial.h5ad')

    # This to alleviate possible normalization discreprencies (more forceful normalization)
#    adata.X = adata.layers['counts'].copy()
#    sc.pp.normalize_total(adata)
#    sc.pp.log1p(adata)

    ## Merge gene features (var) with chromosomal coordinates before this!
    cnv.io.genomic_position_from_gtf('./gencode.v43.annotation.gtf', adata)

    cnv.tl.infercnv(adata,
                    reference_key='phenotype',
                    reference_cat='normal')

    adata.write_h5ad('Epithelial_cells_with_infercnv.h5ad')