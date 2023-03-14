#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --partition=normal
#SBATCH --job-name=spotclean
#SBATCH --output=spotclean.out
#SBATCH --error=spotclean.err
#SBATCH --mail-user=antti.kiviaho@tuni.fi
#SBATCH --mail-type=ALL


# AK 2.12.2022
# SBATCH script to run spotclean on all samples

module load anaconda
source activate r_spatial

Rscript spotclean.R > spotclean_out.txt

