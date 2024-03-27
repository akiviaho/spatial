
import os 
os.chdir('/lustre/scratch/kiviaho/prostate_spatial/')

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad

from utils import get_sample_ids_reorder
from datetime import datetime   


#from statsmodels.stats.multitest import fdrcorrection
from tqdm import tqdm

import warnings
warnings.filterwarnings("ignore")

samples = get_sample_ids_reorder()

color_dict = {
    'Tumor': '#fc8d62',
 'Luminal epithelium': '#8da0cb',
 'Basal epithelium': '#66c2a5',
 'Club epithelium': '#ffd92f',
 'Immune': '#a6d854',
 'Endothelium': '#e78ac3',
 'Fibroblast': '#e5c494',
 'Muscle': '#b3b3b3'
 }

regions = list(color_dict.keys())
region_colors = list(color_dict.values())

all_obs = pd.read_csv('./data/pathology_filtered_obs_with_regions.csv',index_col=0)

sample_list = samples
region_degs_dict = {}

for region in regions:

    print('Processing region: '+ region)
    concat_adata = None

    for sample in tqdm(sample_list, desc="Processing samples", unit="sample"):

        sample_class = all_obs.loc[(all_obs['sample_id']==sample),'sample_class'].iloc[0]
        
        if (sample_class == 'TRNA') | (sample_class == 'NEADT'):
            slide = sc.read_h5ad('./data/visium_with_regions/'+sample+'_with_regions.h5ad')
            slide_subs = slide[slide.obs['predicted_region']==region].copy()
            if slide_subs.shape[0] >= 10:

                if concat_adata is None:
                    concat_adata = slide_subs.copy()
                else:
                    concat_adata = ad.concat([concat_adata,slide_subs],join='outer')

    # Run differential gene expression analysis
    sc.tl.rank_genes_groups(concat_adata,groupby='sample_class',use_raw=False,reference='TRNA',method='wilcoxon')

    # Write the results to an excel file
    res_df = sc.get.rank_genes_groups_df(concat_adata,group='NEADT')

    res_df = res_df[~res_df['names'].str.contains('RPL|RPS')] # Take out ribosomal genes
    res_df = res_df[~res_df['names'].str.startswith('MT')] # Take out mitochondrial genes

    # Save the deg res into a dictionary
    region_degs_dict[region] = res_df

    print('')


# Save the results into an excel-file

# Create a Pandas Excel writer using the file name
writer = pd.ExcelWriter('./supplementary_tables/NEADT_vs_TRNA_degs_by_region.xlsx', engine='xlsxwriter')

# Iterate through each key-value pair in the dictionary
for key, value in region_degs_dict.items():
    # Write each dataframe to a separate sheet in the Excel file
    value.to_excel(writer, sheet_name=key)

# Save and close the Excel writer
writer.save()
