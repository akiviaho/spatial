#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --partition=normal
#SBATCH --job-name=SCS-qc
#SBATCH --output=SCS-qc-%j.out
#SBATCH --error=SCS-qc-%j.err
#SBATCH --mail-user=antti.kiviaho@tuni.fi
# #SBATCH --mail-type=END

# AK 23.2.2023
# SBATCH script to run single cell refernce qc

module load anaconda
source activate scib-pipeline-R4.0

python -u /lustre/scratch/kiviaho/prostate_spatial/scripts/single-cell-reference-data-qc.py