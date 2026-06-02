#!/bin/bash -l
#SBATCH --job-name=mg-coeur-rb
#SBATCH --partition=hackathon
#SBATCH --reservation=hackathon
#SBATCH --time=00:10:00
#SBATCH --mem=96G
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

set -uo pipefail

module purge || true
module load nvhpc/25.5

run_id=${SLURM_JOB_ID:-manual}
run_root="runs/${run_id}_mg_coeur_rollback"
namelist="../../../coeur_unigrid.nml"

mkdir -p "${run_root}/mg"

(
    cd "${run_root}/mg"
    ../../../mg_coeur "${namelist}" \
        > ../../../mg_coeur_rollback_run.log 2>&1
)
