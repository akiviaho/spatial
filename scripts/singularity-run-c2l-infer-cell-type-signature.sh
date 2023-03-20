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
#SBATCH --mem=100G
#SBATCH --mail-type=END
#SBATCH --mail-user=antti.kiviaho@tuni.fi

module load CUDA/11.2
echo "CUDA loaded"

singularity exec \
  --no-home  \
  --nv \
  -f \
  -e \
  -B /lustre/scratch/kiviaho/prostate_spatial:/working_directory \
  /home/ak431480/apps/singularity/cell2location-v0.06-alpha.sif \
  /bin/bash -c "cd /working_directory && chmod 777 /tmp/numba_cache && NUMBA_CACHE_DIR=/tmp/numba_cache && python -u scripts/c2l-infer-cell-type-reference.py"