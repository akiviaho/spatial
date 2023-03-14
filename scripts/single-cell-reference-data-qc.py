# Author: Antti Kiviaho
# Date: 22.2.2023
# A script for running normalization and sample integration clustering.
# Uses the scbi integration environment and pipeline:
#
# 1. Cell and gene filtering
# 2. scran normalization through R interface using
# 3. batch-aware scaling with scib
# 4. batch-aware HVGs with scib
# 5. scanorama integration into PCA, clustering, UMAP

import numpy as np
import anndata as ad
import scanpy as sc
import pandas as pd
import seaborn as sns
import scib
import matplotlib.pyplot as plt
from scipy import sparse
from pathlib import Path
from utils import load_from_pickle, save_to_pickle
import warnings
warnings.filterwarnings('ignore')
import os

# Make sure we're operating in the correct directory
os.chdir('/lustre/scratch/kiviaho/prostate_spatial')


def qc_filters(adata, remove_doublets=True):
    # requires scib-pipline-R4.0 conda environment !
    # import scib
    # Filter out cells by using a hybrid of the original publications thresholds
    sc.pp.filter_cells(adata, min_counts=600)
    sc.pp.filter_cells(adata, min_genes = 300)
    sc.pp.filter_genes(adata, min_counts= 10)
    # Leave out cells with > 20% mitochondrial reads
    adata = adata[adata.obs.pct_counts_mt < 20, :]
    if remove_doublets:
        sc.external.pp.scrublet(adata)
        adata = adata[adata.obs['predicted_doublet']==False]
    
    return adata

if __name__ == '__main__':

    datasets = ['dong_2020','chen_2021','cheng_2022','chen_2022','song_2022','hirz_2023']
    adata_dict = {}
    for dataset_id in datasets:
        adata = sc.read_h5ad('./sc-reference/'+dataset_id+'/adata_obj.h5ad')
        adata_dict[dataset_id] = adata


    # Produce QC plots as done at https://scanpy-tutorials.readthedocs.io/en/latest/spatial/integration-scanorama.html
    # Save the QC plots to a path
    dir_path = './plots/qc-plots-for-sc-reference'
    Path(dir_path).mkdir(parents=True, exist_ok=True)
    for name in datasets:
        adata = adata_dict[name]
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        fig, axs = plt.subplots(1, 4, figsize=(24, 6))
        fig.suptitle(f"Covariates for filtering: {name}")

        sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
        sns.histplot(
            adata.obs["total_counts"][adata.obs["total_counts"] < 20000],
            kde=False,
            bins=40,
            ax=axs[1],
        )
        sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
        sns.histplot(
            adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
            kde=False,
            bins=60,
            ax=axs[3],
        )
        fig.tight_layout()
        fig.set_dpi(200)
        plt.savefig(dir_path+'/'+name+'_qc_metrics'+'.png')
        plt.clf()
        

    for dset in datasets:
        print('Processing dataset: '+dset)
        adata = adata_dict[dset].copy()
        if not sparse.issparse(adata.X):
            adata.X = sparse.csr_matrix(adata.X)
        
        adata = qc_filters(adata)
        
        print('Normalizing dataset: '+dset)
        scib.preprocessing.normalize(adata,precluster=False, sparsify=False)
        # add ids to the data for use after data concatenation
        adata.obs['dataset'] = dset
        adata_dict[dset] = adata
        del adata
        print(dset+' done!')

    print('saving the processed data to pickle...')
    save_to_pickle(adata_dict,'normalized_sc_6_datasets.pickle')
    
    # Produce QC plots as done at https://scanpy-tutorials.readthedocs.io/en/latest/spatial/integration-scanorama.html
    # Save the QC plots to a path
    dir_path = './plots/qc-plots-for-sc-reference-filtered'
    Path(dir_path).mkdir(parents=True, exist_ok=True)
    for name in datasets:
        adata = adata_dict[name]
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        fig, axs = plt.subplots(1, 4, figsize=(24, 6))
        fig.suptitle(f"Covariates for filtering: {name}")

        sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
        sns.histplot(
            adata.obs["total_counts"][adata.obs["total_counts"] < 20000],
            kde=False,
            bins=40,
            ax=axs[1],
        )
        sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
        sns.histplot(
            adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
            kde=False,
            bins=60,
            ax=axs[3],
        )
        fig.tight_layout()
        fig.set_dpi(200)
        plt.savefig(dir_path+'/'+name+'_qc_metrics'+'.png')
        plt.clf()
