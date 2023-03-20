#!/bin/bash
#SBATCH -t 1-23:59:00
#SBATCH -J spatial-scvi
#SBATCH -o spatial-scvi.out.%j
#SBATCH -e spatial-scvi.err.%j
#SBATCH --partition=gpu
#SBATCH --gres=gpu:teslav100:1
# #SBATCH --nodelist=nag16
#SBATCH --exclude=meg[10-12],nag[01-09]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=200G
#SBATCH --mail-type=END
#SBATCH --mail-user=antti.kiviaho@tuni.fi

module load CUDA/11.2
echo "CUDA loaded"

# load conda environment
module load anaconda
source activate scib-pipeline-R4.0
echo "conda activated"

python -u ./scripts/spatial-scvi-integrate.py