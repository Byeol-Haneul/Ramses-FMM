#!/bin/bash -l
#SBATCH --job-name=sedov3d
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=100G
#SBATCH --gres=gpu:1
#SBATCH --constraint=gpu80

module purge
module load nvhpc/25.5

nsys profile \
    --delay=60 \
    --duration=20 \
    --trace=cuda,nvtx,osrt \
    --output=sedov3d_profile \
    ./ramses3d sedov3d_amr.nml > run.log