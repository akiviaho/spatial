#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH --mem=60G
#SBATCH --job-name=seg-temp
#SBATCH --time=02:00:00
#SBATCH --output=seg-temp.out

RUNNING_DIR=/lustre/scratch/kiviaho/prostate_spatial/tools/segmentation
cd /lustre/scratch/kiviaho/prostate_spatial/tools/segmentation

module load miniconda3

source activate seg
python ${RUNNING_DIR}/seg-sample.py

 
