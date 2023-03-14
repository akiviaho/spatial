# Author: Antti Kiviaho
# Date: 8.3.2023
#
# A script for running infercnv on integrated single cell datasets.
# Make sure to identify reference cells via e.g. clustering.
# Also make sure to 'conda activate infercnvpy'
#

import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
from scripts.utils import load_from_pickle, save_to_pickle

sc.settings.set_figure_params(figsize=(5, 5))

if __name__ == '__main__':

    adata = load_from_pickle('./epithelial-scvi-integrated-6-sc-datasets.pickle')

    ## Merge gene features (var) with chromosomal coordinates before this!
    cnv.io.genomic_position_from_gtf('./gencode.v43.annotation.gtf', adata)

    # Providing scvi-clusters with >20% normal sample content as 'normal'
    cnv.tl.infercnv(
    adata,
    # reference_key='broad_celltypes',
    # reference_cat=['B cell', 'Endothelial', 'Fibroblast', 'Mast', 'Monocytic', 'Myofibroblast', 'T cell'],
    reference_key="VI_clusters",
    reference_cat=['0','1','2','6','19','20','31','45'],
    window_size=250,
)

    save_to_pickle(adata,'./epithelial-scvi-integrated-6-sc-datasets-with-infercnv.pickle')