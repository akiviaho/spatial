#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --partition=normal
#SBATCH --job-name=epith-sc-infercnv
#SBATCH --output=epith-sc-infercnv-%j.out
#SBATCH --error=epith-sc-infercnv-%j.err
#SBATCH --mail-user=antti.kiviaho@tuni.fi
#SBATCH --mail-type=END

# AK 8.3.2023
# SBATCH script to run infercnv on single cell datasets

module load anaconda
source activate infercnvpy

python -u ./single-cell-infercnv.py