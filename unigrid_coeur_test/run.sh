#!/bin/bash -l
#SBATCH --job-name=coeur-unigrid
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=100G
#SBATCH --gres=gpu:1
#SBATCH --constraint=gpu80

set -euo pipefail

module purge || true
module load nvhpc/25.5

run_id=${SLURM_JOB_ID:-manual}
mkdir -p "runs/${run_id}/mg" "runs/${run_id}/fmm"

(
  cd "runs/${run_id}/mg"
  ../../../mg_coeur ../../../coeur_unigrid.nml > ../../../mg_run.log
)

(
  cd "runs/${run_id}/fmm"
  ../../../fmm_coeur ../../../coeur_unigrid.nml > ../../../fmm_run.log
)
