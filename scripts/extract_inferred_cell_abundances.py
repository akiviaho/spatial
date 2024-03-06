# Date: 28.2.2024
# Author: Antti Kiviaho

# A script for extracting the inferred cell type abundances from
# a concatenated anndata object resulting from cell2location

import os 
os.chdir('/lustre/scratch/kiviaho/prostate_spatial/')

import scanpy as sc
import pandas as pd
import anndata as ad
from utils import get_sample_ids_reorder

cell_type_renaming_dict = {
    'mesenchymal epithelium':'sensescent epithelium',
    'interferon signaling epithelium':'club_interferon response epithelium',
    'cancer epithelium':'tumor',
    'intermediate CRPC epithelium':'cycling epithelium 1',
    'cycling epithelium':'cycling epithelium 2',
    'FOSL1 tumor epithelium':'FOSL1 related epithelium',
    'fibroblasts':'fibroblast'
    }



if __name__ == '__main__':

    samples = get_sample_ids_reorder()

    # Change the run_name variable to select the appropriate iteration
    run_name = '20240125' # 20230908 run for the Tampere PC cohort

    # Download data and format cell2location mapping results into obs columns in both aggregated adata and individual slides
    adata_vis = sc.read_h5ad('./c2l-results/visium_adata_with_c2l_mapping_'+run_name+'.h5ad')

    # Copy obsm (cell2location results) to adata object obs
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

    ctype_abundances_all = adata_vis.obsm['q05_cell_abundance_w_sf'].copy()

    # Extract redundant information from the mapping names
    ctype_abundances_all.columns = [c.replace('q05cell_abundance_w_sf_','') for c in ctype_abundances_all.columns]
    
    # Rename some of the older column names
    ctype_abundances_all = ctype_abundances_all.rename(columns=cell_type_renaming_dict)

    ctype_abundances_all.to_csv('./data/c2l_run_'+run_name+'_inferred_celltype_abundances_unfiltered.csv')

    # Copy obsm (cell2location results) to individual slides
    for sample in samples:
        adata_sample = sc.read_h5ad('./data/normalized_visium/'+sample+'_normalized.h5ad')
        vis_subset = adata_vis[adata_sample.obs_names].copy()

        if (vis_subset.obs_names == adata_sample.obs_names).all():

            adata_sample.obsm = vis_subset.obsm.copy()

            # add 5% quantile, representing confident cell abundance, 'at least this amount is present',
            # to adata.obs with nice names for plotting
            adata_sample.uns['mod'] = vis_subset.uns['mod'].copy()
            adata_sample.obs[adata_sample.uns['mod']['factor_names']] = adata_sample.obsm['q05_cell_abundance_w_sf']

            abundance_df = adata_sample.obsm['q05_cell_abundance_w_sf'].copy()

            # Extract redundant information from the mapping names
            abundance_df.columns = [c.replace('q05cell_abundance_w_sf_','') for c in abundance_df.columns]
            
            # Rename some of the older column names
            abundance_df = abundance_df.rename(columns=cell_type_renaming_dict)

            abundance_df.to_csv('./data/inferred_celltype_abundances/'+sample+'_abundances.csv')

            print(sample + ' done')

