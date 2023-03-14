#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --partition=normal
#SBATCH --job-name=sc-integrate
#SBATCH --output=sc-integrate-2.out
#SBATCH --error=sc-integrate-2.err
#SBATCH --mail-user=antti.kiviaho@tuni.fi
# #SBATCH --mail-type=ALL

# AK 23.2.2023
# SBATCH script to run integrate single cell reference datasets

module load anaconda
source activate scib-pipeline-R4.0

python -u ./scanmorama-integrate-sc-datasets.py