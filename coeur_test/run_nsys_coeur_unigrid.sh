#!/bin/bash -l
#SBATCH --job-name=nsys-coeur
#SBATCH --partition=hackathon
#SBATCH --reservation=hackathon
#SBATCH --time=00:20:00
#SBATCH --mem=96G
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

set -uo pipefail

module purge || true
module load nvhpc/25.5

run_id=${SLURM_JOB_ID:-manual}
run_root="runs/${run_id}_nsys_coeur_unigrid_delay90_dur30"
namelist="../../../coeur_unigrid.nml"

mkdir -p "${run_root}/mg" "${run_root}/fmm"

nsys_opts=(
    profile
    --force-overwrite=true
    --trace=cuda,nvtx,osrt
    --sample=none
    --cpuctxsw=none
    --delay=90
    --duration=30
)

(
    cd "${run_root}/mg"
    nsys "${nsys_opts[@]}" --output=mg_coeur_unigrid ../../../mg_coeur "${namelist}" \
        > ../../../mg_nsys_coeur_unigrid_run.log 2>&1
)
mg_status=$?
if [ "${mg_status}" -ne 0 ] && [ -f "${run_root}/mg/mg_coeur_unigrid.nsys-rep" ]; then
    mg_status=0
fi

(
    cd "${run_root}/fmm"
    nsys "${nsys_opts[@]}" --output=fmm_coeur_unigrid ../../../fmm_coeur "${namelist}" \
        > ../../../fmm_nsys_coeur_unigrid_run.log 2>&1
)
fmm_status=$?
if [ "${fmm_status}" -ne 0 ] && [ -f "${run_root}/fmm/fmm_coeur_unigrid.nsys-rep" ]; then
    fmm_status=0
fi

if [ -f "${run_root}/mg/mg_coeur_unigrid.nsys-rep" ]; then
    nsys stats --force-overwrite=true --report cuda_gpu_kern_sum \
        "${run_root}/mg/mg_coeur_unigrid.nsys-rep" \
        > "${run_root}/mg/mg_cuda_gpu_kern_sum.txt" 2>&1 || true
    nsys stats --force-overwrite=true --report nvtx_sum \
        "${run_root}/mg/mg_coeur_unigrid.nsys-rep" \
        > "${run_root}/mg/mg_nvtx_sum.txt" 2>&1 || true
fi

if [ -f "${run_root}/fmm/fmm_coeur_unigrid.nsys-rep" ]; then
    nsys stats --force-overwrite=true --report cuda_gpu_kern_sum \
        "${run_root}/fmm/fmm_coeur_unigrid.nsys-rep" \
        > "${run_root}/fmm/fmm_cuda_gpu_kern_sum.txt" 2>&1 || true
    nsys stats --force-overwrite=true --report nvtx_sum \
        "${run_root}/fmm/fmm_coeur_unigrid.nsys-rep" \
        > "${run_root}/fmm/fmm_nvtx_sum.txt" 2>&1 || true
fi

exit $((mg_status != 0 || fmm_status != 0))
