import os 
os.chdir('/lustre/scratch/kiviaho/prostate_spatial/')

import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from utils import load_from_pickle, get_sample_ids, get_sample_ids, get_sample_crop_coords, save_to_pickle
import gseapy as gp

import seaborn as sns
sns.set_theme(style='white')

import warnings
warnings.filterwarnings("ignore")



samples = get_sample_ids()
sample_crop_coord=get_sample_crop_coords()


if __name__ == '__main__':

########### SSGSEA #############
    """
    adata_slides = load_from_pickle('./data/slides_with_cell_mapping_based_regions.pkl')
    gene_sets = pd.read_csv('./custom_gene_lists_encode_fixed.gmt',sep='\t',index_col=0,header=None).index.tolist()

    custom_gene_set_enrichment_scores = {}

    # Calculate the enrichment scores for custom gene lists
    for sample in samples:
        adata = adata_slides[sample].copy()
        expr_values = pd.DataFrame(data=adata.layers['counts'].todense().copy(),
                                index=adata.obs_names.copy(),
                                columns=adata.var_names.copy(),dtype=int).T

        # txt, gct file input
        ssgsea_res = gp.ssgsea(data=expr_values,
            gene_sets='./custom_gene_lists_encode_fixed.gmt',
            sample_norm_method='rank',
            seed=246347)

        # Save the enrichment score object into dict
        custom_gene_set_enrichment_scores[sample] = ssgsea_res.res2d.copy()

        nes = ssgsea_res.res2d.pivot(index='Term', columns='Name', values='NES')
        nes_T = nes.T

        # Save the normalized enrichment scores into the obs column
        adata_slides[sample].obs = pd.concat([adata.obs.copy(),nes_T],axis=1)

        # Save the dict object
        save_to_pickle(custom_gene_set_enrichment_scores,'./data/spatial_ssgsea_results.pkl')
    """

################################## Reading results if they're ready-made ################################################## 

    gene_sets = pd.read_csv('./custom_gene_lists_encode_fixed.gmt',sep='\t',index_col=0,header=None).index.tolist()
    adata_slides = load_from_pickle('./data/slides_with_cell_mapping_based_regions.pkl')
    sample_crop_coord=get_sample_crop_coords()
    gsea_res = load_from_pickle('./data/spatial_ssgsea_results.pkl')

    # Add the ssgsea enrichment scores to slide obs columns
    for sample in samples:
        nes_T = gsea_res[sample].pivot(index='Term', columns='Name', values='NES').T

        # Check if all elements in gene_sets are columns in df
        missing_cols = set(gene_sets) - set(nes_T.columns)

        # Add missing columns to df with np.nan values
        for col in missing_cols:
            nes_T[col] = np.nan
        
        nes_T = nes_T[gene_sets].astype(float)

        if (adata_slides[sample].obs_names == nes_T.index).all():
            adata_slides[sample].obs = pd.concat([adata_slides[sample].obs.copy(),nes_T],axis=1)


################################## Plotting the results on spatial sections ################################################## 

    # Plot the scores
    for var_to_plot in gene_sets:

        it=0
        concat_obs = pd.DataFrame(columns=adata_slides[samples[it]].obs.columns)
        for slide in adata_slides:
            concat_obs = pd.concat([concat_obs,adata_slides[slide].obs.copy()],axis=0)

        fig, axs = plt.subplots(5, 8, figsize=(24, 15),dpi=120)

        min = concat_obs[var_to_plot].quantile(0.05)
        max = concat_obs[var_to_plot].quantile(0.95)
        median = concat_obs[var_to_plot].quantile(0.50)

        if np.abs(0 - median)<=0.1:
            colormap= 'bwr'
            if (min < 0) & (max > 0):
                center = 0
            else:
                center=median
        elif np.isnan(median):
            continue # No scores for this one, don't plot
        else:
            center = median
            colormap = 'viridis'

        for i in range(5):
            for j in range(8):
                
                if it < len(samples):
                    sc.pl.spatial(adata_slides[samples[it]],color=var_to_plot,title=samples[it],
                                vmin=min, vcenter=center,vmax=max, crop_coord=sample_crop_coord[samples[it]],
                                colorbar_loc=None, cmap=colormap, size=1.3, alpha_img=0.8, legend_loc=None,
                                ax=axs[i,j],show=False
                                )

                    axs[i,j].set_xlabel(None)
                    axs[i,j].set_ylabel(None)
                    it+=1
                else:
                    axs[i,j].set_visible(False)

        plt.tight_layout()

        # create a custom axes for the colorbar
        cax = fig.add_axes([0.66, 0.1, 0.3, 0.04])

        # draw the custom colorbar: if there are negative scores center to 0 and plot using bwr. Otherwise use viridis
        if colormap == 'bwr':
            # create a custom diverging colormap
            cmap = colors.LinearSegmentedColormap.from_list('my_cmap', ['blue', 'white', 'red'])
            # create a symmetric normalization around the center value
            norm_colors = colors.TwoSlopeNorm(vmin=min, vcenter=center, vmax=max)
            fig.colorbar(plt.cm.ScalarMappable(norm=norm_colors, cmap=colormap), cax=cax, orientation='horizontal')

        elif colormap == 'viridis':
            fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(vmin=min, vmax=max), cmap=colormap), cax=cax, orientation='horizontal')

        cax.set_title(var_to_plot)

        plt.savefig('./plots/all_samples_grid/custom_gene_lists/ssgsea/'+var_to_plot+'_on_all_spatial_sections.pdf')
        plt.clf()

