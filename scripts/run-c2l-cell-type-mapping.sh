#!/bin/bash
#SBATCH -t 23:59:00
#SBATCH -J c2l-map
#SBATCH -o c2l-map.out.%j
#SBATCH -e c2l-map.err.%j
#SBATCH --partition=gpu
#SBATCH --gres=gpu:teslav100:1
#SBATCH --constraint=gpumem_32
# #SBATCH --nodelist=nag13
# #SBATCH --exclude=meg[10-12],nag[01-09]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --mail-type=END
#SBATCH --mail-user=antti.kiviaho@tuni.fi

module load CUDA/11.2
echo "CUDA loaded"

export PYTHONNOUSERSITE="literallyanyletters"

module load anaconda
source activate cell2loc_env

python -u ./scripts/c2l-cell-type-mapping.py