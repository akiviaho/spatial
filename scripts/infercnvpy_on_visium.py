import os
os.chdir('/lustre/scratch/kiviaho/prostate_spatial')

import numpy as np
import anndata as ad
import scanpy as sc
import pandas as pd
import infercnvpy as cnv
import matplotlib.pyplot as plt
from pathlib import Path
from utils import load_from_pickle, save_to_pickle, get_sample_ids_reorder
import warnings
warnings.filterwarnings('ignore')

################################################################


chr_order = ['chr'+str(i) for i in np.arange(1,23)]

def sort_chromosomes(vars,chromsomes = chr_order):
    # Vars is a pandas dataframe with chromosomal annotation
    # This function drops redundant genes and sorts the df by chromoal position

    # Drop redundant from variables
    vars = vars.dropna()
    vars = vars[~((vars['chromosome'] == 'chrM')|(vars['chromosome'] == 'chrX') | (vars['chromosome'] == 'chrY'))]

    # Sort variable to genomic pos
    vars['chromosome'] = vars['chromosome'].astype('category').cat.set_categories(chromsomes)
    vars = vars.sort_values(['chromosome','start'])

    return(vars)


################################################################

if __name__ == '__main__':

    adata_slides = load_from_pickle('./data/slides_with_cell_mapping_based_regions.pkl')

    samples = get_sample_ids_reorder()

    # Create the baseline reference from BPHs
    ref_regions = ['Luminal epithelium','Basal epithelium','Intermediate epithelium']
    ref_samples =  []

    # Take the specific regions in these samples as a baseline
    for sample in samples[:4]:
        slide = adata_slides[sample]
        ref_samples.append(slide[slide.obs['predicted_region'].isin(ref_regions)].copy())

    # join the reference into a single anndata object
    ref_adata_concat = ad.concat(ref_samples,join='outer')

    ref_adata_concat.obs['cnv_ref'] = 'baseline'

    # No need to run this each time
    #cnv.io.genomic_position_from_gtf('./gencode.v43.annotation.gtf',ref_adata_concat)

    # Download the pre-formatted genomic positions
    genomic_pos = pd.read_csv('genomic_positions_infercnvpy_formatted.csv',index_col=0)
    genomic_pos = sort_chromosomes(genomic_pos)
    ################################################################
    # Run infercnv

    infercnv_dict = {}
    n_samples = len(samples)
    count = 1

    for sample in samples[4:]:

        slide = adata_slides[sample].copy()
        slide.obs['cnv_ref'] = 'sample'

        slide = ad.concat([slide,ref_adata_concat],join='outer')

        # Get the unnormalized layer and remove redundant layers to save space
        slide.X = slide.layers['counts'].copy()
        del slide.layers['counts']
        del slide.raw

        # Put together the layers
        slide = slide[:,genomic_pos.index]
        if (slide.var_names == genomic_pos.index).all():
            slide.var = genomic_pos.copy()

        # Normalize and log transform the data jointly using scanpy
        sc.pp.normalize_total(slide)
        sc.pp.log1p(slide)

        cnv.tl.infercnv(
        slide,
        reference_key='cnv_ref',
        reference_cat=['baseline'],
        window_size=100,
        step=10
        )

        infercnv_dict[sample] = slide
        print('Sample ' +sample+ ' processed: '+ str(count)+'/'+str(n_samples))
        print('')
        count+=1
        save_to_pickle(infercnv_dict,'dict_with_visium_and_infercnv.pkl')