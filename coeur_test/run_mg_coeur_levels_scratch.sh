#!/bin/bash -l
#SBATCH --job-name=mg-coeur-lvl
#SBATCH --partition=hackathon
#SBATCH --reservation=hackathon
#SBATCH --time=00:20:00
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
run_root="${scratch_base}/ct/${run_id}_mg_n${nstepmax}/L${level}"
namelist="${run_root}/coeur.nml"

mkdir -p "${run_root}"
cp "${repo_root}/coeur_test/coeur_unigrid.nml" "${namelist}"

# One AMR oct stores 8 cells. Keep MG headroom similar to the existing level
# sweep, but cap the cache floor so smaller levels do not under-allocate.
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
/usr/bin/time -f "%e" -o time_mg.txt \
    "${repo_root}/coeur_test/mg_coeur" "coeur.nml" \
    > mg_run.log 2>&1
