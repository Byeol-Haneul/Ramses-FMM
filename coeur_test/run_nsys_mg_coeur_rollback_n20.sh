#!/bin/bash -l
#SBATCH --job-name=nsys-mg-coeur-n20
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
run_root="runs/${run_id}_nsys_mg_coeur_rollback_n20"
namelist="../coeur_unigrid_n20.nml"

mkdir -p "${run_root}/mg"
cp coeur_unigrid.nml "${run_root}/coeur_unigrid_n20.nml"
sed -i '/^ncontrol=/a nstepmax=20' "${run_root}/coeur_unigrid_n20.nml"

(
    cd "${run_root}/mg"
    nsys profile --force-overwrite=true --trace=cuda,nvtx,osrt --sample=none --cpuctxsw=none \
        --output=mg_coeur_rollback_n20 ../../../mg_coeur "${namelist}" \
        > ../../../mg_nsys_coeur_rollback_n20_run.log 2>&1
)
mg_status=$?

if [ -f "${run_root}/mg/mg_coeur_rollback_n20.nsys-rep" ]; then
    nsys stats --force-export=true --report cuda_gpu_kern_sum \
        "${run_root}/mg/mg_coeur_rollback_n20.nsys-rep" \
        > "${run_root}/mg/mg_cuda_gpu_kern_sum.txt" 2>&1 || true
    nsys stats --force-export=true --report nvtx_sum \
        "${run_root}/mg/mg_coeur_rollback_n20.nsys-rep" \
        > "${run_root}/mg/mg_nvtx_sum.txt" 2>&1 || true
fi

exit "${mg_status}"
