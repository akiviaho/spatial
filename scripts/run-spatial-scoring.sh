#!/bin/bash
#SBATCH -t 1-23:59:00
#SBATCH -J scanpy_score
#SBATCH -o scanpy_score.out.%j
#SBATCH -e scanpy_score.err.%j
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --mail-type=END
#SBATCH --mail-user=antti.kiviaho@tuni.fi

# load conda environment
module load anaconda
source activate squidpy

python -u ./scripts/spatial_scanpy_enrichment_and_plots.py