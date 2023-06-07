#!/bin/bash
#SBATCH -t 1-23:59:00
#SBATCH -J celltype-nmf
#SBATCH -o celltype-nmf.out.%j
#SBATCH -e celltype-nmf.err.%j
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --mail-type=END
#SBATCH --mail-user=antti.kiviaho@tuni.fi

# load conda environment
module load anaconda
source activate scib-pipeline-R4.0
echo "conda activated"

python -u ./scripts/nmf-on-celltype.py \
--filename Myeloid.h5ad \
--n_var_genes 2000