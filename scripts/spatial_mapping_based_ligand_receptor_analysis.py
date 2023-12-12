import os 
os.chdir('/lustre/scratch/kiviaho/prostate_spatial/')

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import squidpy as sq


from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
from utils import load_from_pickle, save_to_pickle, get_sample_ids
from matplotlib.backends.backend_pdf import PdfPages
from itertools import permutations

import warnings
warnings.filterwarnings("ignore")

# Define functions
def get_spot_interfaces(dat, cluster_of_interest, interaction_cluster, annotation_key = 'tissue_region',added_key='proximity_analysis'):

    # Create an observation column for spatial segmentation
    dat.obs[added_key] = 'Background'
    distance_mat = dat.obsp['spatial_distances'].todense()

    for idx, obs_name in enumerate(dat.obs_names):
        cl = dat.obs[annotation_key][idx]

        if cl == cluster_of_interest:

            first_nhbor_idxs = np.where(distance_mat[:,idx]==1.0)[0] # Get first-term neighbor indices

            try:
                # If try fails, there are no matching clusters as keys in value_counts
                n_cl_neighbors = dat[first_nhbor_idxs].obs[annotation_key].value_counts()[cl] # find first-term neighbor cluster annotations POSSIBLE ERROR IF CL NOT IN DICT

                # Added this clause to control that only those with 'close' interactions with the interaction cluster are included. 
                all_nhbor_indices = np.where(distance_mat[:, idx] != 0)[0]

                # Downgraded the number of required first neighbors to two
                if (n_cl_neighbors >= 1) & (sum(dat.obs[annotation_key][all_nhbor_indices] == interaction_cluster) >= 3):
                    dat.obs.at[obs_name,added_key] = cl

            except:
                continue

    # Make a second loop to make sure the final cluster-of-interest annotations
    # are what's used to define proximal spots
    for idx, obs_name in enumerate(dat.obs_names):
        cl = dat.obs[added_key][idx]

        if cl == cluster_of_interest:
            
            all_nhbor_indices = np.where(distance_mat[:, idx] != 0)[0] 

            # Get the indices where neighboring spots are not the interest cluster 
            indices = np.where((dat.obs[added_key][all_nhbor_indices] != cl) & (dat.obs[annotation_key][all_nhbor_indices] == interaction_cluster))[0]

            # Update the 'proximity_analysis' column for the specific indices
            dat.obs.loc[dat.obs_names[all_nhbor_indices[indices]], added_key] = 'Proximal ' + interaction_cluster

    # Modify the colors to maintain the original cluster color
    dat.obs[added_key] = dat.obs[added_key].astype('category')

    return(dat)

if __name__ == '__main__':


    adata_slides = load_from_pickle('./data/clustered_visium_data.pickle')
    samples = get_sample_ids()

    # Change the run_name variable to select the appropriate iteration
    run_name = '20230908'

    cell_mapping_dat = sc.read_h5ad('c2l_mapping_as_anndata_'+run_name+'.h5ad')
    cell_types = list(cell_mapping_dat.var_names)

    # Copy the annotation column into each member of the adata_slides object
    for sample in samples:
        slide = adata_slides[sample]   
        slide.obs['tissue_region'] = cell_mapping_dat.obs.loc[slide.obs_names]['tissue_region'].astype('category')

    colors_dict = cell_mapping_dat.uns['tissue_region_mapping_colors']
    tissue_region_names = list(colors_dict.keys())

    region_combinations = list(permutations(tissue_region_names,2))
    i = 1

    for sample in samples:
        # Using three rings you get 6 immediate neigbors, 12 second neigbors and 18 third neighbors
        sq.gr.spatial_neighbors(slide,n_neighs=6,n_rings=3)

    for region_pair in region_combinations:
        ## Subset proximal spots, calculate receptor-ligand analysis and plot
        source = region_pair[0]
        target = region_pair[1]

        # Color order dependes on which category comes first in alphabetical order
        if source < 'Proximal ':
            cols = ['#919191',colors_dict[source],colors_dict[target]]
        elif 'Proximal ' < source:
            cols = ['#919191',colors_dict[target],colors_dict[source]]

        interaction_dict = {}

        for sample in samples:

            slide = adata_slides[sample]

            if (slide.obs['tissue_region'].str.contains(source).any()) & (slide.obs['tissue_region'].str.contains(target).any()):


                slide = get_spot_interfaces(slide, source, target)
            
                # Make sure there is an interface in the sample
                if len(slide.obs['proximity_analysis'].cat.categories) == 3:
                    sq.gr.ligrec(
                        slide,
                        n_perms=100,
                        cluster_key="proximity_analysis",
                        show_progress_bar = False
                    )

                    proximal_spots = 'Proximal ' + target

                    pvals = slide.uns['proximity_analysis_ligrec']['pvalues'][source][proximal_spots]
                    means = slide.uns['proximity_analysis_ligrec']['means'][source][proximal_spots]
                    tuple_array = pd.DataFrame(means[pvals<0.01][means>1]).index.values

                    interaction_dict[sample] = tuple_array

                    ## Plotting ##
                    slide.uns['proximity_analysis_colors'] = cols

                    # set figure axis size and dpi
                    fig, ax = plt.subplots(figsize=(8, 8), dpi=120)

                    # create spatial plot
                    sc.pl.spatial(slide,color='proximity_analysis',size=1.3,alpha=0.8, ax = ax, show= False, title= sample)
                    plt.tight_layout()

                    # create filename with sample name
                    filename = 'plots/receptor_ligand_interaction_analysis/' + sample + '_'+ source +'_to_'+ target +'_clusters_communication.pdf'

                    # create output folder if it doesn't exist
                    if not os.path.exists(os.path.dirname(filename)):
                        os.makedirs(os.path.dirname(filename))

                    # save plot to pdf with filename
                    with PdfPages(filename) as pdf:
                        pdf.savefig(fig)
                        plt.clf()

        if len(interaction_dict) > 0:
            save_to_pickle(interaction_dict,'./data/'+source+'_to_'+target+'_ligand_receptor_proximity_interactions.pickle')
        
        if i == 1 | i%5==0:
            print(str(i)/str(len(region_combinations)))
            i+=1