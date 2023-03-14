
import scanpy as sc
import scib
import anndata as ad
import matplotlib.pyplot as plt
from scipy import sparse
from scripts.utils import load_from_pickle, save_to_pickle
from pathlib import Path


if __name__ == "__main__":

    adata_dict = load_from_pickle('./normalized_sc_datasets.pickle')
    adata_dict.pop('song_2022')
    adata_concat = ad.concat(adata_dict)

    adata_concat.obs.dataset = adata_concat.obs.dataset.astype('category')
    scib.preprocessing.scale_batch(adata_concat,batch='dataset')
    print('Scaling done...')

    hvg_list = scib.preprocessing.hvg_batch(adata_concat,batch_key='dataset',target_genes=2000,flavor='seurat',adataOut=False)
    print('HVGs calculated...')

    adata_scanorama = scib.integration.scanorama(adata_concat,batch='dataset',hvg=hvg_list)
    print('Scanorama ready...')

    sc.pp.neighbors(adata_scanorama, use_rep="X_scanorama",random_state=745634)
    sc.tl.umap(adata_scanorama,random_state=745634)
    sc.tl.leiden(adata_scanorama, key_added="clusters")
    print('NN graph, UMAP & Leiden ready...')

    save_to_pickle(adata_scanorama,'scanorama_integrated_normalized_hirz_chen_datasets.pickle')
    print('Scanorama integration saved...')

    dir_path = './plots/sc-reference-integration-umaps-song-excluded'
    Path(dir_path).mkdir(parents=True, exist_ok=True)

    sc.set_figure_params(dpi=200,frameon=True,figsize=(12,12))
    fig,ax1 = plt.subplots(1,1)
    sc.pl.umap(adata_scanorama, color=["dataset"],ax=ax1,size=20,show=False)
    fig.tight_layout()
    fig.savefig(dir_path+'/'+'umap-with-datasets.png')
    plt.clf()

    sc.set_figure_params(dpi=200,frameon=True,figsize=(14,12))
    fig,ax1 = plt.subplots(1,1)
    sc.pl.umap(adata_scanorama, color=["clusters"],ax=ax1,size=20,show=False)
    fig.tight_layout()
    fig.savefig(dir_path+'/'+'umap-with-clusters.png')
    plt.clf()

    sc.set_figure_params(dpi=200,frameon=True,figsize=(16,12))
    fig,ax1 = plt.subplots(1,1)
    sc.pl.umap(adata_scanorama, color=["cells"],ax=ax1,size=20,show=False)
    fig.tight_layout()
    fig.savefig(dir_path+'/'+'umap-with-cell-types-hirz.png')
    plt.clf()


