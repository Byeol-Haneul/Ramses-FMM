#!/bin/bash
#SBATCH --job-name=ramses-collapse
#SBATCH --output=/home/jl4415/scratch/collapse/mg_orfmmm/collapse_%j.out
#SBATCH --error=/home/jl4415/scratch/collapse/mg_orfmmm/collapse_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --mail-user=jl4415@princeton.edu
#SBATCH --mail-type=END,FAIL

# ============================
# Environment Setup
# ============================
set -euo pipefail

module purge
module load openmpi/gcc/4.1.2

cd /home/jl4415/mini-ramses || exit 1

export SLURM_CPU_BIND=cores
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp

NP=96
RUN_ROOT="/home/jl4415/scratch/collapse/mg_orfmmm"
FMM_DIR="${RUN_ROOT}/fmm"
MG_DIR="${RUN_ROOT}/mg"
EXE_DIR="${RUN_ROOT}/bin"
NML_FILE="/home/jl4415/mini-ramses/namelist/boss_bodenheimer1979.nml"
FMM_EXE_SRC="/home/jl4415/mini-ramses/bin/collapse_fmm3d"
MG_EXE_SRC="/home/jl4415/mini-ramses/bin/collapse_mg3d"
FMM_EXE="${EXE_DIR}/collapse_fmm3d"
MG_EXE="${EXE_DIR}/collapse_mg3d"

mkdir -p "$FMM_DIR" "$MG_DIR" "$EXE_DIR"

cp -p "$FMM_EXE_SRC" "$FMM_EXE"
cp -p "$MG_EXE_SRC" "$MG_EXE"
chmod u+x "$FMM_EXE" "$MG_EXE"
cp -p "$NML_FILE" "${RUN_ROOT}/boss_bodenheimer1979.nml"

# ============================
# FMM run
# ============================
cd "$FMM_DIR" || exit 1
srun -n "$NP" "$FMM_EXE" "$NML_FILE" \
    > fmm_${SLURM_JOB_ID}.out \
    2> fmm_${SLURM_JOB_ID}.err

# ============================
# MG run
# ============================
cd "$MG_DIR" || exit 1
srun -n "$NP" "$MG_EXE" "$NML_FILE" \
    > mg_${SLURM_JOB_ID}.out \
    2> mg_${SLURM_JOB_ID}.err
