#!/bin/bash
#SBATCH -t 1-23:59:00
#SBATCH -J spatial-cluster
#SBATCH -o spatial-cluster.out.%j
#SBATCH -e spatial-cluster.err.%j
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --mail-type=END
#SBATCH --mail-user=antti.kiviaho@tuni.fi

# load conda environment
module load anaconda
source activate squidpy
echo "conda activated"

python -u ./scripts/gex-based-joint-clustering.py

# ./scripts/spatial-mapping-based-clustering.py 
