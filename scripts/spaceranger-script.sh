#!/bin/bash
#SBATCH --time=5-23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --partition=normal
#SBATCH --job-name=V13Y08-090-slide
#SBATCH --output=V13Y08-090-spaceranger.out
#SBATCH --mail-user=antti.kiviaho@tuni.fi
#SBATCH --mail-type=END


# AK 29.11.2022
# Batch script to run spaceranger on a node
# Refer to sample-info.xlsx for variables (in Dropbox)

SLIDE=V13Y08-090
SAMPLES=(MET_A3 MET_GP12 MET_A14 MET_A16)
AREAS=(A1 B1 C1 D1)

module purge
module load compbio/spaceranger

for i in "${!SAMPLES[@]}"; do
  SAMPLE="${SAMPLES[$i]}"
  AREA="${AREAS[$i]}"
  spaceranger count --id="${SAMPLE}" \
  --transcriptome=/lustre/scratch/kiviaho/prostate_spatial/data/refdata-cellranger-GRCh38-3.0.0/ \
  --fastqs=/lustre/scratch/kiviaho/prostate_spatial/data/reads/2024-met-reads/X208SC24033399-Z01-F001/01.RawData/${AREA}_${SAMPLE} \
  --sample="${SAMPLE}" \
  --image="/lustre/scratch/kiviaho/prostate_spatial/data/images/${SAMPLE}_20x_zoom.tiff" \
  --slide=$SLIDE \
  --area=$AREA \
  --loupe-alignment="/lustre/scratch/kiviaho/prostate_spatial/data/loupe-alignments/${SLIDE}-${AREA}.json"
done

module purge
