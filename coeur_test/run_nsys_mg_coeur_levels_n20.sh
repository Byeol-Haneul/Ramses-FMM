#!/bin/bash -l
#SBATCH --job-name=nsys-mg-coeur-lvl
#SBATCH --partition=hackathon
#SBATCH --reservation=hackathon
#SBATCH --time=00:10:00
#SBATCH --mem=96G
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=7-10%1

set -uo pipefail

module purge || true
module load nvhpc/25.5

level=${SLURM_ARRAY_TASK_ID:-${1:-8}}
if [ "${level}" -lt 7 ] || [ "${level}" -gt 10 ]; then
    echo "Expected level in [7, 10], got ${level}" >&2
    exit 2
fi

run_id=${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}
task_id=${SLURM_ARRAY_TASK_ID:-${level}}
run_root="runs/${run_id}_nsys_mg_coeur_levels_n20/L${level}"
namelist="../coeur_unigrid_L${level}_n20.nml"

mkdir -p "${run_root}/mg"
cp coeur_unigrid.nml "${run_root}/coeur_unigrid_L${level}_n20.nml"

# One AMR oct stores 8 cells 8 cells. Keep headroom for the MG coarser tree,
# matching the successful level-8 setup where ngridmax is about 4x the base octs.
base_ngrid=$((1 << (3 * level - 3)))
ngridmax=$((4 * base_ngrid))
ncachemax=$((ngridmax / 4))
if [ "${ncachemax}" -lt 2000000 ]; then
    ncachemax=2000000
fi

sed -i \
    -e "/^ncontrol=/a nstepmax=20" \
    -e "s/^levelmin=.*/levelmin=${level}/" \
    -e "s/^levelmax=.*/levelmax=${level}/" \
    -e "s/^ngridmax=.*/ngridmax=${ngridmax}/" \
    -e "s/^ncachemax=.*/ncachemax=${ncachemax}/" \
    "${run_root}/coeur_unigrid_L${level}_n20.nml"

(
    cd "${run_root}/mg"
    nsys profile --force-overwrite=true --trace=cuda,nvtx,osrt --sample=none --cpuctxsw=none \
        --output="mg_coeur_L${level}_n20" ../../../../mg_coeur "${namelist}" \
        > "../../../mg_nsys_coeur_L${level}_n20_run.log" 2>&1
)
mg_status=$?

if [ -f "${run_root}/mg/mg_coeur_L${level}_n20.nsys-rep" ]; then
    nsys stats --force-export=true --report cuda_gpu_kern_sum \
        "${run_root}/mg/mg_coeur_L${level}_n20.nsys-rep" \
        > "${run_root}/mg/mg_cuda_gpu_kern_sum.txt" 2>&1 || true
    nsys stats --force-export=true --report nvtx_sum \
        "${run_root}/mg/mg_coeur_L${level}_n20.nsys-rep" \
        > "${run_root}/mg/mg_nvtx_sum.txt" 2>&1 || true
fi

exit "${mg_status}"
