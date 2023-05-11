import os 
os.chdir('/lustre/scratch/kiviaho/prostate_spatial')
import numpy as np
import scanpy as sc
import scib
#from sklearn.decomposition import NMF
import nimfa
import upsetplot as ups
from matplotlib import pyplot as plt
import argparse

def save_to_pickle(obj,filename):
    import pickle
    with open(filename, 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)


def calculate_matrix_orders(arr):
    # get the indices that would sort each row of the array
    sort_indices = np.argsort(-arr, axis=1)

    # create an array to mark the sorted order
    sorted_order = np.empty_like(sort_indices)
    rows, cols = np.indices(arr.shape)
    sorted_order[rows, sort_indices] = cols

    # replace each entry in the original array with its index in the sorted order
    result = sorted_order.astype(int)

    return result


### Params ###
parser = argparse.ArgumentParser()
parser.add_argument('--filename', type=str, default='Fibroblast.h5ad', help='Name of the file to be used')
parser.add_argument('--n_var_genes', type=int,default=1000, help='The target number of HVGs')
args = parser.parse_args()

filename = args.filename
n_var_genes = args.n_var_genes
n_comps_range =  np.arange(2,16)

path_for_plots = 'sc_modules/plots/'
path_for_data = 'sc_modules/data/'

### Read and process the cell type subsetted data ###

dat = sc.read_h5ad(filename)
dat.X = dat.layers['counts'].copy()


sc.pp.filter_genes(dat,min_counts=1)
sc.pp.normalize_total(dat)
sc.pp.scale(dat)
sc.pp.highly_variable_genes(dat,n_top_genes=n_var_genes,flavor='seurat_v3',
                            subset=True, layer='counts')

# Replace negative entries with zero
dat.X[dat.X < 0] = 0

### Calulating the NMF ###

for n_comps in n_comps_range:

#    model = NMF(n_components=n_comps, init='nndsvda', solver='mu', 
#                beta_loss='kullback-leibler',random_state=35735)
#    W = model.fit_transform(dat.X)
#    H = model.components_


    # Perform non-negative matrix factorization using nsNMF
    f = nimfa.Nsnmf(dat.X, rank=n_comps)
    f_fit = f()
        
    # Extract the resulting matrices into variables W and H
    W = np.array(f_fit.basis())
    H = np.array(f_fit.coef())
    
    res = {'samples':dat.obs,'genes':dat.var,'sample_weights':W,'gene_weights':H}

    save_to_pickle(res,path_for_data+filename.replace('.h5ad','_'+str(n_comps)+'_modules_nsnmf_results.pickle'))

    gene_w = res['gene_weights'].T
    genes = list(res['genes'].index)

    mat_1 = calculate_matrix_orders(gene_w.T).T
    mat_2 = calculate_matrix_orders(gene_w)

    # Iterate through the factors (columns)
    genes_by_factors = {}
    for factor in range(gene_w.shape[1]):
        valid_genes = list()
        # Iterate through the genes, starting from the highest weighted gene of this factor
        for i in np.argsort(mat_1[:,factor]): 
            if mat_2[i,factor] ==0: # Is this the factor the gene effects the most?
                
                # If yes, add it to the list of genes
                valid_genes.append(genes[i])
            else:
                # if not, stop adding genes into the list, move on to the next factor
                break 
        # Append the list of valid genes into a dictionary under the appropriate key
        genes_by_factors['factor'+str(factor)] = valid_genes


    ups_data = ups.from_contents(genes_by_factors)
    ax_dict = ups.UpSet(ups_data, subset_size='count',show_counts=True).plot()

    plt.savefig(path_for_plots+'upsetplot_'+filename.replace('.h5ad','')+'_'+str(n_comps)+'_modules.png')

