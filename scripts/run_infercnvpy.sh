#!/bin/bash
#SBATCH -t 1-23:59:00
#SBATCH -J infercnvpy
#SBATCH -o infercnvpy.out.%j
#SBATCH -e infercnvpy.err.%j
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --mail-type=END
#SBATCH --mail-user=antti.kiviaho@tuni.fi

# load conda environment
module load anaconda
source activate infercnvpy
echo "conda activated"

python -u ./scripts/infercnvpy_on_visium.py
