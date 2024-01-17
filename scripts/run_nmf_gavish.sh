#!/bin/bash
#SBATCH -t 1-23:59:00
#SBATCH -J Endothelial-nmf
#SBATCH -o Endothelial-nmf.out.%j
#SBATCH -e Endothelial-nmf.err.%j
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G

# load conda environment
module load anaconda
source activate scib-pipeline-R4.0
echo "conda activated"

python -u ./scripts/nmf_on_celltype_gavish.py \
--filename ./nmf_annotation/Endothelial.h5ad \
--n_var_genes 2000 
