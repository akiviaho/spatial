#!/bin/bash
#SBATCH -t 03:59:00
#SBATCH -J rec_lig_analysis
#SBATCH -o rec_lig_analysis.out.%j
#SBATCH -e rec_lig_analysis.err.%j
#SBATCH --partition=test
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=80G
#SBATCH --mail-type=END
#SBATCH --mail-user=antti.kiviaho@tuni.fi

# load conda environment
module load anaconda
source activate squidpy
echo "conda activated"

python -u ./scripts/spatial_mapping_based_ligand_receptor_analysis.py
