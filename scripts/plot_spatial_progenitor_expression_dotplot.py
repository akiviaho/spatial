import os 
os.chdir('/lustre/scratch/kiviaho/prostate_spatial/')

import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from utils import get_sample_ids_reorder, get_sample_crop_coords, get_sample_id_mask
from datetime import datetime   

import seaborn as sns
sns.set_theme(style='white')

#from statsmodels.stats.multitest import fdrcorrection
from tqdm import tqdm

from scipy.stats import zscore

import warnings
warnings.filterwarnings("ignore")

if __name__ == '__main__':

    samples = get_sample_ids_reorder()
    sample_crop_coord = get_sample_crop_coords()
    sample_id_masks = get_sample_id_mask()


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


    # Fetch normalized expression and percentage of spots that express any given gene

    def fetch_region_normalized_expression_and_percentage(gene_markers, sample_list, group_id, regions_list = regions, use_unnormalized = True):
        '''
        Returns the per-sample mean expression of genes of interest across all samples that are defined in sample_list
        inside a region defined by region parameter.
        '''

        final_expression_df = pd.DataFrame()
        final_percentage_df = pd.DataFrame()
        for region in regions_list:

            print('Region: ' + region)
            region_expression_df = pd.DataFrame()
            #region_pct_df = pd.DataFrame()

            for sample in tqdm(sample_list, desc="Processing samples", unit="sample"):
                
                slide = sc.read_h5ad('./data/visium_with_regions/'+sample+'_with_regions.h5ad')
                slide_subs = slide[slide.obs['predicted_region']==region].copy()

                # Control for the number of data points belonging to a class
                if slide_subs.shape[0] >= 10:
                    
                    present_genes = [g for g in gene_markers if g in slide_subs.var_names]
                    missing_genes = [g for g in gene_markers if g not in slide_subs.var_names]

                    genes_all_arr_order_match = present_genes + missing_genes

                    #if use_unnormalized:
                    expr_without_missing_genes = slide_subs[:,present_genes].layers['counts'].copy().todense()
                    #else:
                        #expr_without_missing_genes = slide_subs[:,present_genes].X.copy()


                    expr_all_spots = np.concatenate((expr_without_missing_genes,np.full((expr_without_missing_genes.shape[0],len(missing_genes)), np.nan)),axis=1)

                    expr_as_df = pd.DataFrame(data=expr_all_spots.T,index=genes_all_arr_order_match,columns=slide_subs.obs_names)
                    expr_as_df = expr_as_df.loc[gene_markers]
                    
                    # Concatenate the counts from a single sample into a dataframe with all the spots
                    region_expression_df = pd.concat([region_expression_df,expr_as_df],axis=1)

                del slide, slide_subs

            # Added on 4.3.2024 prior to plotting the whole thing
            # Previously NaN's inflated the percentage, as they were non-zero
            region_expression_df = region_expression_df.fillna(0)
            # Put the "percentage of spots expressed in" information into a dataframe
            region_pct_df = pd.DataFrame((region_expression_df != 0).sum(axis=1)/region_expression_df.shape[1],columns=[region])

            
            # Here you concatenate the mean of all valid spots into a dataframe
            final_expression_df = pd.concat([final_expression_df,region_expression_df.mean(axis=1)],axis=1)
            final_percentage_df = pd.concat([final_percentage_df,region_pct_df],axis=1)
        

        final_expression_df.columns = [r + ' ' + group_id for r in regions_list]
        final_percentage_df.columns = [r + ' ' + group_id for r in regions_list]
        return(final_expression_df, final_percentage_df)




    # Get lists of samples in corresponding groupings
    #normal_samples = get_sample_ids_reorder(['BPH'])
    unt_samples = get_sample_ids_reorder(['untreated'])
    trt_samples = get_sample_ids_reorder(['bicalutamide','goserelin','degarelix','degarelix_apalutamide'])
    #crpc_samples = get_sample_ids_reorder(['CRPC'])

    ar_signaling_genes = ['AR','ABCC4','FKBP5','KLK3','MAF','NKX3-1','PMEPA1','div1', # AR regulated genes
    'KRT5','KRT15','TP63','div2', # Canonical basal markers
    'MMP7','PIGR','LTF','SCGB1A1','SCGB3A1','div3', # Club-like markers
    #'KRT4','PSCA','WFDC2','CYP2F1','TSPAN8','CLU','KRT19','KLF5','ANXA3','PPP1R1B','S100A11','KRT8','ATP1B1', # All progenitor markers (Baures et al. Cancers 2022)
    'S100A11','WFDC2','CLU','KRT19','KLF5','ATP1B1','KRT4', # Select representative progenitor markers for main plot
    'EGFR','MET' # Putative stemness-type receptors #,'IGF1R','ERBB4','ERBB2'
    ] 

    expr_unt, pct_unt = fetch_region_normalized_expression_and_percentage(ar_signaling_genes, unt_samples,group_id='(TRNA)',regions_list=regions[:4])
    expr_trt, pct_trt = fetch_region_normalized_expression_and_percentage(ar_signaling_genes, trt_samples,group_id='(NEADT)',regions_list=regions[:4])

    

    #plot_df = pd.concat([normal_ar,untreated_ar,treated_ar,crpc_ar],axis=1).T
    expr_df = pd.concat([expr_unt,expr_trt],axis=1).T
    pct_df = pd.concat([pct_unt,pct_trt],axis=1).T
    expr_df = expr_df.apply(lambda x: zscore(x, nan_policy='omit'))
    regions_mod = expr_df.index.tolist()

    # Format to long
    plot_df = expr_df.melt(ignore_index=False)
    plot_df.columns = ['gene','expression']

    # Format to long
    pct_df = pct_df.melt(ignore_index=False)
    pct_df.columns = ['gene','percentage']

    plot_df['percentage'] = pct_df['percentage'].copy()

    plot_df = plot_df.reset_index(names='region')


    # Create the dotplot
    sns.set_theme(style='white')

    width = 12
    height = 4

    fig, ax = plt.subplots(figsize=(width, height))
    yticks_list = list(np.arange(2,(len(regions_mod)*2)+2,2)[::-1])


    # Get control over interactions order and gap
    plot_df['region_y'] = plot_df['region'].map(dict(zip(regions_mod,yticks_list)))
    sns.scatterplot(x='gene', y='region_y', size='percentage', hue='expression', 
                    hue_norm=(-2,2),data=plot_df, sizes=(30, 300), palette='bwr', ax=ax,legend=True,
                    )

    plt.ylim(0,yticks_list[0]+2)
    plt.yticks(yticks_list,regions_mod)
    plt.xticks(rotation=45)
    plt.legend(loc='center left',handlelength=1.5, handleheight=1.5, bbox_to_anchor=(1.05, 0.5))
    plt.tight_layout()

    plt.savefig('plots/normalized_gene_expression_heatmaps/ar_club_markers_expression_nan_fixed.pdf')

