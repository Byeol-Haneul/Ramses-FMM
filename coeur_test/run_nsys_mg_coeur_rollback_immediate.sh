#!/bin/bash -l
#SBATCH --job-name=nsys-mg-coeur-rb
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
run_root="runs/${run_id}_nsys_mg_coeur_rollback_immediate"
namelist="../../../coeur_unigrid.nml"

mkdir -p "${run_root}/mg"

(
    cd "${run_root}/mg"
    nsys profile --force-overwrite=true --trace=cuda,nvtx,osrt --sample=none --cpuctxsw=none \
        --output=mg_coeur_rollback_unigrid ../../../mg_coeur "${namelist}" \
        > ../../../mg_nsys_coeur_rollback_immediate_run.log 2>&1
)
mg_status=$?

if [ -f "${run_root}/mg/mg_coeur_rollback_unigrid.nsys-rep" ]; then
    nsys stats --force-export=true --report cuda_gpu_kern_sum \
        "${run_root}/mg/mg_coeur_rollback_unigrid.nsys-rep" \
        > "${run_root}/mg/mg_cuda_gpu_kern_sum.txt" 2>&1 || true
    nsys stats --force-export=true --report nvtx_sum \
        "${run_root}/mg/mg_coeur_rollback_unigrid.nsys-rep" \
        > "${run_root}/mg/mg_nvtx_sum.txt" 2>&1 || true
fi

exit "${mg_status}"
