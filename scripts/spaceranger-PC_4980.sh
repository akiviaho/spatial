#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --partition=normal
#SBATCH --job-name=sr-PC-4980
#SBATCH --output=PC_4980_slurm.out
#SBATCH --error=PC_4980_slurm.err
#SBATCH --mail-user=antti.kiviaho@tuni.fi
#SBATCH --mail-type=ALL


# AK 24.5.2021
# Batch script to run spaceranger on a node
# Refer to sample-info.xlsx for variables (in Dropbox)

SHORT=PC_4980
SAMPLE=PC_4980
SLIDE=V10F03-107
AREA=D1



module purge
module load compbio/spaceranger

spaceranger count --id=$SAMPLE \
--transcriptome=/lustre/scratch/kiviaho/prostate_spatial/data/refdata-cellranger-GRCh38-3.0.0/ \
--fastqs=/lustre/scratch/kiviaho/prostate_spatial/data/reads/ \
--sample="${AREA}_${SAMPLE}" \
--image="/lustre/scratch/kiviaho/prostate_spatial/data/images/${SAMPLE}.tiff" \
--slide=$SLIDE \
--area=$AREA \
--loupe-alignment="/lustre/scratch/kiviaho/prostate_spatial/data/loupe-alignments/${SAMPLE}.json" > "./${SAMPLE}.out"

module purge
