#!/bin/bash
#SBATCH -t 1-23:59:00
#SBATCH -J c2l-infer
#SBATCH -o c2l-infer.out.%j
#SBATCH -e c2l-infer.err.%j
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

export PYTHONNOUSERSITE="literallyanyletters"

module load anaconda
source activate cell2loc_env

python -u ./scripts/c2l-cell-type-mapping.py