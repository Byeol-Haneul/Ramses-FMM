#!/bin/bash -l
#SBATCH --job-name=fmm-smoke
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu=20G
#SBATCH --gres=gpu:1
#SBATCH --constraint=gpu80

set -euo pipefail

module purge || true
module load nvhpc/25.5

run_id=${SLURM_JOB_ID:-manual}
mkdir -p "runs/${run_id}"

(
  cd "runs/${run_id}"
  ../../../unigrid_coeur_test/fmm_coeur ../../fmm_two_particle.nml > ../../fmm_run.log
)
