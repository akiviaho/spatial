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
sample_crop_coord = get_sample_crop_coords()


if __name__ == '__main__':

######### Scanpy scoring method ###########

    # Re-download slides to get rid of the orig scores
    adata_slides = load_from_pickle('./data/slides_with_cell_mapping_based_regions.pkl')
    gene_set_df = pd.read_csv('./custom_gene_lists_encode_fixed.gmt',sep='\t',index_col=0,header=None).T

    custom_gene_set_scanpy_scores = {}

    # Calculate the scanpy scores for custom gene lists
    for sample in samples:
        slide = adata_slides[sample]
        for col in gene_set_df.columns:

            sc.tl.score_genes(slide,gene_set_df[col],score_name=col,random_state=2531035)

        custom_gene_set_scanpy_scores[sample] = slide.obs[gene_set_df.columns].copy()

        # Save the dict object
        save_to_pickle(custom_gene_set_scanpy_scores,'./data/spatial_scanpy_score_results.pkl')


    # Plot the scores
    for var_to_plot in gene_set_df.columns:

        it=0
        concat_obs = pd.DataFrame(columns=adata_slides[samples[it]].obs.columns)
        for slide in adata_slides:
            concat_obs = pd.concat([concat_obs,adata_slides[slide].obs.copy()],axis=0)

        fig, axs = plt.subplots(5, 8, figsize=(24, 15),dpi=120)

        min = concat_obs[var_to_plot].quantile(0.05)
        max = concat_obs[var_to_plot].quantile(0.95)
        median = concat_obs[var_to_plot].quantile(0.50)

        if np.abs(0 - median)<=0.1: # If median is near 0, plot with bwr colormap and 0 center
            colormap= 'bwr'
            center = 0
        elif np.isnan(median):
            continue # No scores for this one, don't plot
        else: # If median is far from zero, center at median and plot with viridis
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

        plt.savefig('./plots/all_samples_grid/custom_gene_lists/scanpy_score/'+var_to_plot+'_on_all_spatial_sections.pdf')
        plt.clf()
