import os 
os.chdir('/lustre/scratch/kiviaho/prostate_spatial/')

import scanpy as sc
import numpy as np
import squidpy as sq
import pandas as pd
import anndata as ad

import matplotlib.pyplot as plt
from utils import load_from_pickle, get_sample_ids
import matplotlib as mpl

import seaborn as sns
sns.set_theme()

sc.set_figure_params(figsize=(6,6))

import warnings
warnings.filterwarnings("ignore")


adata_slides = load_from_pickle('./data/clustered_visium_data.pickle')
nmf_derived_genes = pd.read_csv('nmf_derived_genes_epi_endo_fibro.csv')['0']
samples = get_sample_ids()

# Concatenate the visium data and subset it by using a set of NMF derived genes
adata_concat = ad.concat(adata_slides)
mask = ([g in pd.Index(nmf_derived_genes) for g in adata_concat.var_names])
adata_concat = adata_concat[:,mask]

# Calculate the spatial neighborhood weights for each individual sample
for s in samples:
    sq.gr.spatial_neighbors(adata_slides[s],n_rings=2, coord_type="grid", n_neighs=6,transform='cosine')
spatial_connectivity_matrix = ad.concat(adata_slides,pairwise=True).obsp['spatial_connectivities'].copy()
adata_concat.obsp['spatial_connections'] = spatial_connectivity_matrix

# Set the clustering parameters
proximity_weight = 0.3
res = 0.7

# Define the joint adjacency weighting
sc.pp.pca(adata_concat, n_comps=10)
sc.pp.neighbors(adata_concat, use_rep='X_pca', random_state=2692056)#
joint_adj = spatial_connectivity_matrix*proximity_weight + adata_concat.obsp['connectivities']
adata_concat.obsp['joint_adjacencies'] = joint_adj
sc.tl.leiden(adata_concat,adjacency=joint_adj,key_added='common_spatial_clusters',resolution=res,random_state=2692056)

# Save the data for following use
adata_concat.write('jointly_clustered_visium_as_anndata_20230803.h5ad')

