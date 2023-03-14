#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --partition=normal
#SBATCH --job-name=CRPC_531
#SBATCH --mail-user=antti.kiviaho@tuni.fi
#SBATCH --mail-type=ALL


# AK 29.11.2022
# Batch script to run spaceranger on a node
# Refer to sample-info.xlsx for variables (in Dropbox)

SAMPLE=CRPC_531
SLIDE=V12U07-009
AREA=D1



module purge
module load compbio/spaceranger

spaceranger count --id="${SAMPLE}" \
--transcriptome=/lustre/scratch/kiviaho/prostate_spatial/data/refdata-cellranger-GRCh38-3.0.0/ \
--fastqs=/lustre/scratch/kiviaho/prostate_spatial/data/reads/2022-reads/${SAMPLE} \
--sample="${SAMPLE}" \
--image="/lustre/scratch/kiviaho/prostate_spatial/data/images/${SAMPLE}_20x_zoom.tiff" \
--slide=$SLIDE \
--area=$AREA \
--loupe-alignment="/lustre/scratch/kiviaho/prostate_spatial/data/loupe-alignments/${SLIDE}-${AREA}.json" > "./${SAMPLE}.out"

module purge
