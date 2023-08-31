import os 
os.chdir('/lustre/scratch/kiviaho/prostate_spatial')
import numpy as np
import pandas as pd
import scanpy as sc
from itertools import combinations
import nimfa
import argparse


def save_to_pickle(obj,filename):
    import pickle
    with open(filename, 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

def nmf_preprocessing(dat,n_var_genes=2000):

    dat.X = dat.layers['counts'].copy()
    sc.pp.filter_genes(dat,min_counts=10)
    sc.pp.normalize_total(dat)
    sc.pp.scale(dat)
    sc.pp.highly_variable_genes(dat,n_top_genes=n_var_genes,flavor='seurat_v3',
                                subset=True, layer='counts')

    # Replace negative entries with zero
    dat.X[dat.X < 0] = 0

    return(dat)



################# Params ########################

parser = argparse.ArgumentParser()
parser.add_argument('--filename', type=str, help='Name of the file to be used')
parser.add_argument('--n_var_genes', type=int,default=2000, help='The target number of HVGs')
parser.add_argument('--n_top_genes_per_nmf_component', type=int,default=50, help='The number of top contributing genes to save for each component')

args = parser.parse_args()

filename = args.filename
n_hvgs = args.n_var_genes
n_top_genes = args.n_top_genes_per_nmf_component

# Threshold for retaining cells
min_cells = 100

# As defined in Gavish et al. Nature 2023
nmf_comp_range = np.arange(4,10)

# Define saving name
save_name = filename.replace('.h5ad','_top_'+str(n_top_genes)+'_genes_for_nmf_components_4_to_9.pkl')

if __name__ == '__main__':

    # Download data
    ctype_dat = sc.read_h5ad(filename)

    # Initialize dict object
    nmf_sample_dict = {}

    # Get samples that pass the min_cells threshold
    val_counts = ctype_dat.obs['sample'].value_counts()
    qualified_samples = list(val_counts[val_counts >= min_cells].index)
    dq_samples = list(val_counts[val_counts < min_cells].index)


    for s in qualified_samples:
    
        # Subset the data to only include a single sample
        dat = ctype_dat[ctype_dat.obs['sample'] == s]

        # Do a sample-specific preprocessing
        dat = nmf_preprocessing(dat,n_hvgs)

        # Create an empty dictionary to store the top genes for each column
        top_genes_dict = {}

        for n_comps in nmf_comp_range:

            # Perform non-negative matrix factorization using nsNMF
            f = nimfa.Nsnmf(dat.X, rank=n_comps)
            f_fit = f()
                
            # Extract the resulting matrices into variables W and H
            W = np.array(f_fit.basis())
            H = np.array(f_fit.coef())

            res = {'samples':dat.obs,'genes':dat.var,'sample_weights':W,'gene_weights':H}

            gene_w = res['gene_weights'].T
            genes = list(res['genes'].index)

            # Get the number of columns in gene_w
            num_columns = gene_w.shape[1]

            # Iterate over each column of gene_w
            for col in range(num_columns):
                # Get the top 50 gene indexes for the current column
                top_50_indexes = np.argsort(gene_w[:, col])[-n_top_genes:][::-1]
                
                # Get the top 50 genes for the current column using the indexes
                top_genes = [genes[idx] for idx in top_50_indexes]
                
                # Add the list of top genes to the dictionary
                key = str(n_comps)+'_'+ str(col)
                top_genes_dict[key] = top_genes

        print('Processed sample ' + s + ' (' + str(val_counts[s]) + ' cells)')

        nmf_sample_dict[s] = pd.DataFrame(top_genes_dict)

        # Saving the final dict object with 39-column dataframes for each sample   
        save_to_pickle(nmf_sample_dict,save_name)

    print('')
    print('Disqualified samples: ')
    for s in dq_samples:
        print(s + ' (' + str(val_counts[s]) + ' cells)')


