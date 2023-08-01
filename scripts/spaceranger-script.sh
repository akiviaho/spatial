#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --partition=normal
#SBATCH --job-name=V12N14-331-slide
#SBATCH --output=V12N14-331-spaceranger.out
#SBATCH --mail-user=antti.kiviaho@tuni.fi
#SBATCH --mail-type=END


# AK 29.11.2022
# Batch script to run spaceranger on a node
# Refer to sample-info.xlsx for variables (in Dropbox)

SLIDE=V12N14-331
SAMPLES=(PC_06_04077_OIK_ANT_2 PC_06_04581_OIK_POST_0 PC_06_16086_VAS_POST_2)
AREAS=(A1 B1 D1)

module purge
module load compbio/spaceranger

for i in "${!SAMPLES[@]}"; do
  SAMPLE="${SAMPLES[$i]}"
  AREA="${AREAS[$i]}"
  spaceranger count --id="${SAMPLE}" \
  --transcriptome=/lustre/scratch/kiviaho/prostate_spatial/data/refdata-cellranger-GRCh38-3.0.0/ \
  --fastqs=/lustre/scratch/kiviaho/prostate_spatial/data/reads/2023-reads/${SAMPLE} \
  --sample="${SAMPLE}" \
  --image="/lustre/scratch/kiviaho/prostate_spatial/data/images/${SAMPLE}_20x_zoom.tiff" \
  --slide=$SLIDE \
  --area=$AREA \
  --loupe-alignment="/lustre/scratch/kiviaho/prostate_spatial/data/loupe-alignments/${SLIDE}-${AREA}.json"
done

module purge
