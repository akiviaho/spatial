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
# Date: 14.3.2023
#
# A script for running scvi integration on spatial data
# There are 17 samples in total

if __name__ == "__main__":


    # Load the data
    adata = load_from_pickle('./data/normalized_visium_data.pickle')
    adata = ad.concat(adata)

    # Preprocess and scale
    adata.obs.sample_id = adata.obs.sample_id.astype('category')
    scib.preprocessing.scale_batch(adata,batch='sample_id')
    print('Scaling done...')

    hvg_list = scib.preprocessing.hvg_batch(adata,batch_key='sample_id',target_genes=3000,flavor='seurat',adataOut=False)
    print('HVGs calculated...')

    print('CUDA is available: ' + str(torch.cuda.is_available()))
    print('GPUs available: ' + str(torch.cuda.device_count()))

    print('Initiating training on GPU ...')
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample_id")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")

    vae.train(use_gpu=True)

    vae.save('visium-scvi-model')
    adata.obsm["X_scVI"] = vae.get_latent_representation()

    sc.pp.neighbors(adata, use_rep="X_scVI",random_state=3463453)
    sc.tl.umap(adata,random_state=3463453)
    sc.tl.leiden(adata, key_added="VI_clusters")
    print('NN graph, UMAP & Leiden ready...')

    save_to_pickle(adata,'visium_after_scvi_integration.pickle')
    print('SCVI integration saved...')
