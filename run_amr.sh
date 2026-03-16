#!/bin/bash
#SBATCH --job-name=ramses-amr-bench
#SBATCH --output=/home/jl4415/mini-ramses/bench_%j.out
#SBATCH --error=/home/jl4415/mini-ramses/bench_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --exclusive
#SBATCH --mail-user=jl4415@princeton.edu
#SBATCH --mail-type=END,FAIL

# ============================
# Environment Setup
# ============================
module purge
module load openmpi/gcc/4.1.2

cd /home/jl4415/mini-ramses || exit 1

export SLURM_CPU_BIND=cores
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp

timestamp() { date +"%Y%m%d_%H%M%S"; }

# ============================
# Benchmark settings
# ============================
NP_LIST=(4 8 16 32 64)
NML_FILE="namelist/benchmark/lvl9_13.nml"
while true; do
    sleep 1
done

# ============================
# OUTER LOOP
# ============================
for rep in {1..5}; do
    RUNSTAMP=$(timestamp)
    OUTDIR="/home/jl4415/mini-ramses/benchmark_amr/run_${RUNSTAMP}"
    mkdir -p "$OUTDIR"

    echo "=============================================="
    echo "   AMR BENCHMARK RUN $rep    ($RUNSTAMP)"
    echo "=============================================="

    for NP in "${NP_LIST[@]}"; do
        echo "--- Strong scaling: MPI ranks $NP | FMM-AMR ---"
        srun -n "$NP" ./ramses_fmm_amr "$NML_FILE"

        if [ -f time_fmm.txt ]; then
            mv time_fmm.txt "${OUTDIR}/strong_fmm_amr_lvl9_14_P${NP}_${RUNSTAMP}.txt"
        else
            echo "Warning: time_fmm.txt not found for FMM-AMR, NP=$NP"
        fi

        #echo "--- Strong scaling: MPI ranks $NP | MG ---"
        #srun -n "$NP" ./ramses_mg "$NML_FILE"

        #if [ -f time_mg.txt ]; then
        #    mv time_mg.txt "${OUTDIR}/strong_mg_lvl9_14_P${NP}_${RUNSTAMP}.txt"
        #else
        #    echo "Warning: time_mg.txt not found for MG, NP=$NP"
        #fi
    done

    echo "===== RUN $rep COMPLETE ====="
done

echo "===== ALL BENCHMARK RUNS COMPLETE ====="