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


# Download data and format cell2location mapping results into obs columns in both aggregated adata and individual slides
#c2l-results/cell2location_map_20230721/
adata_vis = sc.read_h5ad('./c2l-results/visium_adata_with_c2l_mapping_20230721.h5ad')

adata_vis.obs.joint_leiden_clusters = adata_vis.obs.sample_id.astype(str) + '_' + adata_vis.obs.joint_leiden_clusters.astype(str)
adata_vis.obs['joint_leiden_clusters'] = pd.Categorical(adata_vis.obs['joint_leiden_clusters'])


adata_slides = load_from_pickle('./data/clustered_visium_data.pickle')
samples = get_sample_ids()

# Copy obsm (cell2location results) to adata object obs
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# Copy obsm (cell2location results) to individual slides
for sample in samples:
    vis_subset = adata_vis[adata_vis.obs['sample_id']==sample]
    
    if (vis_subset.obs_names == adata_slides[sample].obs_names).all():
        adata_slides[sample].obsm = vis_subset.obsm.copy()
        
        # add 5% quantile, representing confident cell abundance, 'at least this amount is present',
        # to adata.obs with nice names for plotting
        adata_slides[sample].uns['mod'] = vis_subset.uns['mod'].copy()
        adata_slides[sample].obs[adata_slides[sample].uns['mod']['factor_names']] = adata_slides[sample].obsm['q05_cell_abundance_w_sf']

# Creating an anndata-structure for inferred cell numbers

# Filters if necessary
ctypes = adata_vis.uns['mod']['factor_names']
#ctypes = [ctype for ctype in ctypes if 'epithelial' in ctype or 'fibroblast' in ctype]

cell_mapping_dat = sc.AnnData(adata_vis.obs[ctypes])
if (cell_mapping_dat.obs.index == adata_vis.obs.index).all():
    cell_mapping_dat.obs = adata_vis.obs[['sample_id','joint_leiden_clusters']]


if (adata_vis.obs.index == cell_mapping_dat.obs.index).all():
    adata_vis.obs['inferred_cell_counts'] = np.round(np.sum(cell_mapping_dat.X,axis=1),0).tolist()
    print('Inferred cell counts copied!')


# extract first string before '-' or '_'
cell_mapping_dat.obs['phenotype'] = cell_mapping_dat.obs['sample_id'].str.extract('([^\-_]+)', expand=False)

cell_mapping_dat.layers['counts'] = cell_mapping_dat.X.copy()

for s in samples:
    sq.gr.spatial_neighbors(adata_slides[s],n_rings=2, coord_type="grid", n_neighs=6,transform='cosine')
spatial_connectivity_matrix = ad.concat(adata_slides,pairwise=True).obsp['spatial_connectivities'].copy()
cell_mapping_dat.obsp['spatial_connections'] = spatial_connectivity_matrix

proximity_weight = 0.2
res = 0.8

# Define the joint adjacency weighting
sc.pp.pca(cell_mapping_dat, n_comps=10)
sc.pp.neighbors(cell_mapping_dat, use_rep='X_pca', random_state=2692056)#
joint_adj = spatial_connectivity_matrix*proximity_weight + cell_mapping_dat.obsp['connectivities']
cell_mapping_dat.obsp['joint_adjacencies'] = joint_adj
sc.tl.leiden(cell_mapping_dat,adjacency=joint_adj,key_added='mapping_based_spatial_clusters',resolution=res,random_state=2692056)

# Save the data for following use
cell_mapping_dat.write('c2l_mapping_as_anndata_20230721.h5ad')

