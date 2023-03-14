from utils import load_from_pickle, save_to_pickle
import scanpy as sc
import anndata as ad
import scib
import scvi
import seaborn as sns
import torch
import os
from matplotlib import pyplot as plt

# Author: Antti Kiviaho
# Date: 27.2.2023
#
# A script for running scvi integration on single cell datasets
# Dong 2020, Chen 2021, Cheng 2022, Chen 2022, Song 2022, Hirz 2023

if __name__ == "__main__":

    ##############    
    # Some params
    prefix = ''
    filter_celltypes = False
    which_cells = ['']

    ###############

    # Load the data
    adata = load_from_pickle('./sc-reference/normalized_sc_6_datasets_with_annot.pickle')
    adata = ad.concat(adata)
    adata.obs_names_make_unique() # Some duplicate index persists

    if filter_celltypes:
        adata = adata[adata.obs['broad_celltypes'].isin(which_cells)]


    # Preprocess and scale
    adata.obs.dataset = adata.obs.dataset.astype('category')
    scib.preprocessing.scale_batch(adata,batch='dataset')
    print('Scaling done...')

    hvg_list = scib.preprocessing.hvg_batch(adata,batch_key='dataset',target_genes=2000,flavor='seurat',adataOut=False)
    print('HVGs calculated...')

    print('CUDA is available: ' + str(torch.cuda.is_available()))
    print('GPUs available: ' + str(torch.cuda.device_count()))

    print('Initiating training on GPU ...')
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="dataset")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")

    vae.train(use_gpu=True)

    vae.save(prefix+'-scvi-model')
    adata.obsm["X_scVI"] = vae.get_latent_representation()

    sc.pp.neighbors(adata, use_rep="X_scVI",random_state=745634)
    sc.tl.umap(adata,random_state=745634)
    sc.tl.leiden(adata, key_added="VI_clusters")
    print('NN graph, UMAP & Leiden ready...')

    save_to_pickle(adata,prefix+'-scvi-integrated-6-sc-datasets.pickle')
    print('SCVI integration saved...')

    #### EXTEND TO SCANVI (use harmonized annotations) ####

    lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="broad_celltypes",
    unlabeled_category="NA")

    lvae.train(max_epochs=20, n_samples_per_label=100, use_gpu=True)
    lvae.save(prefix+'-scanvi-model')

    adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)

    sc.pp.neighbors(adata, use_rep="X_scANVI",random_state=745634)
    sc.tl.umap(adata,random_state=745634)
    sc.tl.leiden(adata, key_added="ANVI_clusters")
    print('NN graph, UMAP & Leiden ready...')

    save_to_pickle(adata,prefix+'-scANVI-integrated-6-sc-datasets.pickle')
    print('scANVI integration saved...')

