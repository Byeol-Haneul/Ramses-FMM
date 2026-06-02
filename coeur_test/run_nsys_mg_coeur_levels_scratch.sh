#!/bin/bash -l
#SBATCH --job-name=nsys-mg-coeur
#SBATCH --partition=hackathon
#SBATCH --reservation=hackathon
#SBATCH --time=00:30:00
#SBATCH --mem=96G
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=5-10%1

set -euo pipefail

module purge || true
module load nvhpc/25.5

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    repo_root="${SLURM_SUBMIT_DIR}"
else
    repo_root=$(cd "$(dirname "$0")/.." && pwd)
fi
scratch_base=${SCRATCH_BASE:-/scratch/gpfs/TEYSSIER/jl4415}
level=${SLURM_ARRAY_TASK_ID:-${1:-8}}
nstepmax=${NSTEPMAX:-20}

if [ "${level}" -lt 5 ] || [ "${level}" -gt 10 ]; then
    echo "Expected level in [5, 10], got ${level}" >&2
    exit 2
fi

run_id=${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}
run_root="${scratch_base}/ct/${run_id}_nsys_mg_n${nstepmax}/L${level}"
namelist="${run_root}/coeur.nml"
profile_name="mg_coeur_L${level}_n${nstepmax}"

mkdir -p "${run_root}"
cp "${repo_root}/coeur_test/coeur_unigrid.nml" "${namelist}"

# One AMR oct stores 8 cells. Keep MG headroom similar to the successful
# scratch sweep and avoid regular dumps so profiler output stays manageable.
base_ngrid=$((1 << (3 * level - 3)))
ngridmax=$((4 * base_ngrid))
if [ "${level}" -eq 10 ]; then
    ngridmax=1181116007
fi
ncachemax=$((ngridmax / 4))
if [ "${ncachemax}" -lt 2000000 ]; then
    ncachemax=2000000
fi

sed -i \
    -e "/^ncontrol=/a nstepmax=${nstepmax}" \
    -e "s/^levelmin=.*/levelmin=${level}/" \
    -e "s/^levelmax=.*/levelmax=${level}/" \
    -e "s/^ngridmax=.*/ngridmax=${ngridmax}/" \
    -e "s/^ncachemax=.*/ncachemax=${ncachemax}/" \
    -e "s/^delta_tout *=.*/delta_tout = 100.0/" \
    -e "/^&OUTPUT_PARAMS/a foutput=0\noutput_part=.false.\noutput_grav=.false.\noutput_hydro=.false.\noutput_amr=.false." \
    "${namelist}"

cd "${run_root}"
nsys profile --force-overwrite=true --trace=cuda,nvtx,osrt --sample=none --cpuctxsw=none \
    --output="${profile_name}" "${repo_root}/coeur_test/mg_coeur" "coeur.nml" \
    > mg_nsys_run.log 2>&1
mg_status=$?

if [ -f "${profile_name}.nsys-rep" ]; then
    nsys stats --force-export=true --report cuda_gpu_kern_sum \
        "${profile_name}.nsys-rep" > mg_cuda_gpu_kern_sum.txt 2>&1 || true
    nsys stats --force-export=true --report nvtx_sum \
        "${profile_name}.nsys-rep" > mg_nvtx_sum.txt 2>&1 || true
fi

exit "${mg_status}"
