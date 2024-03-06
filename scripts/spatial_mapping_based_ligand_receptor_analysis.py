import os 
os.chdir('/lustre/scratch/kiviaho/prostate_spatial/')

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import squidpy as sq


from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
from utils import load_from_pickle, get_sample_ids_reorder, get_sample_crop_coords, get_sample_id_mask, save_to_pickle
from statsmodels.stats.multitest import multipletests
from itertools import combinations
from sklearn.cluster import KMeans

import seaborn as sns
sns.set_theme(style='white')

import warnings
warnings.filterwarnings("ignore")

samples = get_sample_ids_reorder()
sample_crop_coord = get_sample_crop_coords()
sample_id_masks = get_sample_id_mask()


## Define the target and the source
source = 'Club epithelium'
target = 'Muscle'

# Define functions

def get_spot_interfaces(dat, cluster_of_interest, interaction_cluster, annotation_key='predicted_region', added_key='proximity_analysis'):

    # Create an observation column for spatial segmentation
    dat.obs[added_key] = np.nan
    distance_mat = dat.obsp['spatial_distances'].todense()

    for idx, obs_name in enumerate(dat.obs_names):
        cl = dat.obs[annotation_key][idx]

        if cl in [cluster_of_interest, interaction_cluster]:
            first_nhbor_idxs = np.where(distance_mat[:, idx] == 1.0)[0]  # Get first-term neighbor indices

            try:
                n_cl_neighbors = dat[first_nhbor_idxs].obs[annotation_key].value_counts()[cl]  # find first-term neighbor cluster annotations
                all_nhbor_indices = np.where(distance_mat[:, idx] != 0)[0]

                if cl == cluster_of_interest:
                    if (n_cl_neighbors >= 0) & (sum(dat.obs[annotation_key][all_nhbor_indices] == interaction_cluster) >= 2):
                        dat.obs.at[obs_name, added_key] = cl

                elif cl == interaction_cluster:
                    if (n_cl_neighbors >= 0) & (sum(dat.obs[annotation_key][all_nhbor_indices] == cluster_of_interest) >= 2):
                        dat.obs.at[obs_name, added_key] = cl

            except:
                continue

    # Modify the colors to maintain the original cluster color
    dat.obs[added_key] = dat.obs[added_key].astype('category')

    return(dat)

if __name__ == '__main__':

    # This is the publication iteration of regions
    adata_slides = load_from_pickle('./data/slides_with_cell_mapping_based_regions_20240125_tampere_arneo.pkl')

    # Get unique region names
    regions = adata_slides[samples[5]].obs['predicted_region'].cat.categories.tolist()

    # Get unique regions colors
    region_colors = adata_slides[samples[5]].uns['predicted_region_colors']

    # Create color legend
    region_colors = adata_slides[samples[5]].uns['predicted_region_colors']
    color_dict = dict(zip(regions,region_colors))


    valid_samples = []
    it=0
    fig, axs = plt.subplots(8, 6, figsize=(18, 24),dpi=120)

    for i in range(8):
        for j in range(6):

            if it < len(samples) :
                
                slide = adata_slides[samples[it]]
                # Using three rings you get 6 immediate neigbors, 12 second neigbors and 18 third neighbors
                #sq.gr.spatial_neighbors(slide,n_neighs=6,n_rings=3)
        
                slide = get_spot_interfaces(slide, source, target)

                # Qualify sample only if there are 10 or more of both source and target spots
                if not (slide.obs['proximity_analysis'].isna().all()):
                    if (slide.obs['proximity_analysis'].str.contains(source).sum() >= 10) & (slide.obs['proximity_analysis'].str.contains(target).sum() >= 10):
                    
                        ## Plotting ##
                        slide.uns['proximity_analysis_colors'] = [color_dict[cat] for cat in slide.obs['proximity_analysis'].cat.categories]

                        # create spatial plot
                        if 'P320' not in samples[it]:
                            sc.pl.spatial(slide,color='proximity_analysis',title=samples[it],
                                                    crop_coord=sample_crop_coord[samples[it]],
                                                    size=1.5, alpha_img=0, legend_loc=None,na_color='whitesmoke',
                                                    ax=axs[i,j],show=False
                                                    )

                        else:
                            sc.pl.spatial(slide,color='proximity_analysis',title=samples[it],
                                                    size=1.5, alpha_img=0, legend_loc=None,na_color='whitesmoke',
                                                    ax=axs[i,j],show=False
                                                    )


                        axs[i,j].set_xlabel(None)
                        axs[i,j].set_ylabel(None)

                        # Append this sample to the list
                        valid_samples.append(samples[it])
                    else:
                        axs[i,j].set_visible(False)
                else:
                    axs[i,j].set_visible(False)
            else:
                axs[i,j].set_visible(False)
            
            it+=1


    plt.tight_layout()
    plt.savefig('./plots/receptor_ligand_interaction_analysis/'+source+'_'+target+'_proximity_regions.pdf')
    #plt.clf()
    plt.clf()


    #### Second part, using valid_samples to do ligrec #####
    ligrec_dict = {}
    for sample in valid_samples:

        slide = adata_slides[sample].copy()
        slide.obs['proximity_analysis'] = slide.obs['proximity_analysis'].cat.add_categories(['NA'])
        slide.obs['proximity_analysis'] = slide.obs['proximity_analysis'].fillna('NA')

        if (len(slide.obs['proximity_analysis'].cat.categories.tolist()) == 3):
        
            ligrec_res = sq.gr.ligrec(
                slide,
                cluster_key='proximity_analysis',
                clusters = [source,target],
                complex_policy='all',
                show_progress_bar = False,
                n_perms=1000,
                seed=4359345,
                copy=True,
                use_raw=False
            )

            ligrec_dict[sample] = ligrec_res

    # Save the proximal spot annotations
    proximity_spot_annots = pd.DataFrame()
    for s in valid_samples:
        proximity_spot_annots = pd.concat([
            proximity_spot_annots,
            adata_slides[s].obs.copy()],
            axis=0)

    proximity_spot_annots = proximity_spot_annots[['sample_id','predicted_region','proximity_analysis']]
    proximity_spot_annots = proximity_spot_annots.rename(columns={'proximity_analysis':'{}_{}_proximity'.format(source,target)})


    proximity_spot_annots.to_csv('./data/proximity_spot_ids/{}_to_{}_spot_annotation.csv'.format(source,target))
    save_to_pickle(ligrec_dict,'./data/region_ligrec_analysis/'+source+'_'+target+'_slides_with_ligrec.pkl')
    print('{} to {} ligand-receptor interaction analysis has been saved!'.format(source,target))