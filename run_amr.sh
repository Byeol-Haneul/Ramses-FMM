#!/bin/bash
#SBATCH --job-name=ramses-amr-bench
#SBATCH --output=/home/jl4415/mini-ramses/benchmark_amr_stable/bench_%j.out
#SBATCH --error=/home/jl4415/mini-ramses/benchmark_amr_stable/bench_%j.err
#SBATCH --time=00:30:00
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
cd /home/jl4415/mini-ramses/benchmark_amr_stable

# ============================
# Benchmark settings
# ============================
NP_LIST=(64 32 16 8 4 2 1)
NML_FILE_LIST=("../namelist/benchmark/noast/lvl8_13.nml") #"../namelist/benchmark/noast/lvl8_14.nml")

# ============================
# OUTER LOOP
# ============================
for rep in {1}; do
    RUNSTAMP=$(timestamp)

    echo "=============================================="
    echo "   AMR BENCHMARK RUN $rep    ($RUNSTAMP)"
    echo "=============================================="

    # Loop over NML files
    for NML_FILE in "${NML_FILE_LIST[@]}"; do

        # Extract short name (e.g. lvl9_13)
        NML_TAG=$(basename "$NML_FILE" .nml)

        OUTDIR="/home/jl4415/mini-ramses/benchmark_amr_stable/${NML_TAG}/run_${RUNSTAMP}"
        mkdir -p "$OUTDIR"

        echo "##############################################"
        echo "   Running NML: $NML_FILE"
        echo "##############################################"

        for NP in "${NP_LIST[@]}"; do
            echo "--- Strong scaling: MPI ranks $NP | MERGE-FMM-AMR ---"
            srun -n "$NP" ./fmm_amr "$NML_FILE"

            if [ -f time_fmm.txt ]; then
                mv time_fmm.txt "${OUTDIR}/strong_fmm_${NML_TAG}_P${NP}_${RUNSTAMP}.txt"
            else
                echo "Warning: time_fmm.txt not found for FMM-AMR, NP=$NP"
            fi

            echo "--- Strong scaling: MPI ranks $NP | MG ---"
            srun -n "$NP" ./mg_amr "$NML_FILE"

            if [ -f time_mg.txt ]; then
                mv time_mg.txt "${OUTDIR}/strong_mg_${NML_TAG}_P${NP}_${RUNSTAMP}.txt"
            else
                echo "Warning: time_mg.txt not found for MG-AMR, NP=$NP"
            fi
        done

    done

    echo "===== RUN $rep COMPLETE ====="
done

echo "===== ALL BENCHMARK RUNS COMPLETE ====="