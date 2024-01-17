#!/bin/bash/

read -p "Enter python script name: " filename
srun \
    --job-name squidpy_env \
    --partition normal \
    --mem 120G \
    --ntasks 1 \
    --cpus-per-task 8 \
    --time 12:00:00 \
    bash -c "module load anaconda && source activate squidpy && python -u $filename"

