#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=40000
#SBATCH --partition=normal
#SBATCH --job-name=all-sr-aggregate


#AK 8.7.2021
# Script to aggregate multiple runs of spaceranger

module purge
module load compbio/spaceranger

cd /lustre/scratch/kiviaho/prostate_spatial/results/

spaceranger aggr --id=AGGREGATE_ALL \
--csv=/lustre/scratch/kiviaho/prostate_spatial/data/aggregate_all.csv \
 > "./aggregate_all.out"

module purge
