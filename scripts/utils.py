import numpy as np
import anndata as ad
import scanpy as sc
import pickle


def get_sample_ids(sample_types=['BPH','CRPC','PC']):
    samples = []
    if 'BPH' in sample_types:
        samples = samples + ['BPH_651','BPH_665','BPH_688']

    if 'PC' in sample_types:
        samples = samples + ['PC-03-6712','PC_00_16338_II','PC_01_06342_VAS',
        'PC_01_14451_OIK','PC_02_05601_OIK','PC_02_10136_VAS',
        'PC_03_01669_TUTKV','PC_15420OIK','PC_4980','PC_7875OIK']

    if 'CRPC' in sample_types:
        samples = samples + ['CRPC-278','CRPC-489','CRPC-697','CRPC_531']

    return samples

def qc_and_normalize(adata):
    # requires scib-pipline-R4.0 conda environment !
    import scib
    # normalize and calculate leiden clustering
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.filter_cells(adata, min_counts=500) 
    scib.preprocessing.normalize(adata,precluster=False)
    return adata

def spatially_aware_clustering(adata,proximity_weight=0.3,res=1.0):
    import squidpy as sq
    # Define the joint adjacency weighting
    sc.pp.scale(adata)
    sc.pp.highly_variable_genes(adata)
    sc.pp.pca(adata, n_comps=15)
    sc.pp.neighbors(adata, random_state=36345)
    sq.gr.spatial_neighbors(adata, n_rings=2, coord_type="grid", n_neighs=6,transform='cosine')
    joint_adj = adata.obsp['spatial_connectivities']*proximity_weight + adata.obsp['connectivities']
    sc.tl.leiden(adata,adjacency=joint_adj,key_added='joint_leiden_clusters',resolution=res,random_state=36345)
    return adata

def save_to_pickle(obj,filename):
    import pickle
    with open(filename, 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

def load_from_pickle(filename):
    import pickle
    with open(filename, 'rb') as handle:
        obj = pickle.load(handle)
    return obj