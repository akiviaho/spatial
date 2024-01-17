from utils import load_from_pickle, save_to_pickle
import scanpy as sc
import anndata as ad
import scib
import scvi
import seaborn as sns
import torch
import os
from matplotlib import pyplot as plt
from datetime import datetime

# Author: Antti Kiviaho
# Date: 27.2.2023
#
# A script for running scvi integration on single cell datasets
# Henry 2018, Dong 2020, Chen 2021, Cheng 2022, Chen 2022, Song 2022, Wong 2022, Hirz 2023


if __name__ == "__main__":

    ##############    
    # Some params

    current_date = datetime.today().strftime('%Y%m%d')
    # Add a technology column:
    chemistry_version = {'henry_2018':'chromium_v2','dong_2020':'chromium_v2',
                         'chen_2021':'chromium_v2','cheng_2022':'chromium_v2',
                         'chen_2022':'chromium_v2','song_2022':'seq_well3',
                         'wong_2022':'chromium_v1','hirz_2023':'chromium_v2'}
    

    ###############

    # Load the data
    adata = load_from_pickle('./sc-reference/normalized_sc_8_datasets.pickle')

    # Updated for 8 datasets. Run this before concatenation (error fix)
    subs = adata['henry_2018'].copy()
    subs.var_names_make_unique()
    del subs.raw
    adata['henry_2018'] = subs.copy()
    
    # Downsample Hirz to 50k stop it from dominating the integration
    adata['hirz_2023'] = adata['hirz_2023'][adata['hirz_2023'].obs.sample(int(5e4),
    random_state=934537).index]

    adata = ad.concat(adata)
    adata.obs_names_make_unique() # Some duplicate index persists

    # Preprocess and scale
    adata.obs.dataset = adata.obs.dataset.astype('category')
    scib.preprocessing.scale_batch(adata,batch='dataset') # https://doi.org/10.1038/s41592-021-01336-8

    print('Scaling done...')

    # Map chemistry version according to the dataset
    adata.obs['chemistry_version'] = adata.obs['dataset'].map(chemistry_version).astype('category')

    adata.raw = adata

    adata = scib.preprocessing.hvg_batch(adata,batch_key='dataset',target_genes=2000,flavor='seurat',adataOut=True)
    print('HVGs calculated...')

    print('CUDA is available: ' + str(torch.cuda.is_available()))
    print('GPUs available: ' + str(torch.cuda.device_count()))

    print('Initiating training on GPU ...')
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="dataset",categorical_covariate_keys=['chemistry_version'])
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")

    vae.train(use_gpu=True)

    vae.save('scvi_model_'+current_date)
    adata.obsm["X_scVI"] = vae.get_latent_representation()

    sc.pp.neighbors(adata, use_rep="X_scVI",random_state=745634)
    sc.tl.umap(adata,random_state=745634)
    sc.tl.leiden(adata, key_added="VI_clusters")
    print('NN graph, UMAP & Leiden ready...')

    save_to_pickle(adata,'scvi_integrated_8_sc_datasets_'+current_date+'.pickle')
    print('SCVI integration saved...')

    #### EXTEND TO SCANVI (use harmonized annotations) ####
""" 
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

    save_to_pickle(adata,prefix+'-scANVI-integrated-7-sc-datasets.pickle')
    print('scANVI integration saved...')

 """