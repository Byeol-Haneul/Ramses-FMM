#!/bin/bash
#SBATCH --job-name=ramses-bench
#SBATCH --output=/home/jl4415/mini-ramses/bench_%j.out
#SBATCH --error=/home/jl4415/mini-ramses/bench_%j.err
#SBATCH --qos=debug
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
cd /home/jl4415/mini-ramses

# Deterministic CPU binding + MPI backend
export SLURM_CPU_BIND=cores
export OMP_PLACES=cores
export OMP_PROC_BIND=close

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp

CORES_PER_NODE=96
NODES_ALLOCATED=1
MAX_RANKS=$(( CORES_PER_NODE * NODES_ALLOCATED ))
OUTDIR="/home/jl4415/mini-ramses/benchmark_out"
mkdir -p $OUTDIR

timestamp() { date +"%Y%m%d_%H%M%S"; }

# ============================
# Generate powers-of-two ranks
# ============================
generate_ranks() {
    local max=$1
    local ranks=()
    local val=16
    while :; do
        if (( val <= max )); then
            ranks+=("$val")
        else
            ranks+=("$max")
            break
        fi
        # compute next value; if it exceeds max, stop next iteration
        local next=$((val * 2))
        if (( next > max )); then
            ranks+=("$max")
            break
        fi
        val=$next
    done
    echo "${ranks[@]}"
}


WEAK_RANKS=($(generate_ranks $MAX_RANKS))
STRONG_RANKS=("${WEAK_RANKS[@]}")

# ============================
# Strong Scaling (fixed level 9)
# ============================
STRONG_LVL=10
echo "===== STRONG SCALING | Level $STRONG_LVL ====="
for NP in "${STRONG_RANKS[@]}"; do
    if (( NP > MAX_RANKS )); then
        echo "Skipping strong scaling ranks $NP: exceeds allocation $MAX_RANKS"
        continue
    fi

    echo "--- Strong scaling: MPI ranks $NP ---"

    # Timed Runs
    srun --mpi=pmix --cpu-bind=cores -n ${NP} /home/jl4415/mini-ramses/ramses_mg \
         /home/jl4415/mini-ramses/namelist/benchmark/lvl${STRONG_LVL}_fmm1.nml
    mv time_mg.txt ${OUTDIR}/strong_mg_lvl${STRONG_LVL}_P${NP}_$(timestamp).txt

    srun --mpi=pmix --cpu-bind=cores -n ${NP} /home/jl4415/mini-ramses/ramses_fmm \
         /home/jl4415/mini-ramses/namelist/benchmark/lvl${STRONG_LVL}_fmm1.nml
    mv time_fmm.txt ${OUTDIR}/strong_fmm1_lvl${STRONG_LVL}_P${NP}_$(timestamp).txt

    srun --mpi=pmix --cpu-bind=cores -n ${NP} /home/jl4415/mini-ramses/ramses_fmm \
         /home/jl4415/mini-ramses/namelist/benchmark/lvl${STRONG_LVL}_fmm2.nml
    mv time_fmm.txt ${OUTDIR}/strong_fmm2_lvl${STRONG_LVL}_P${NP}_$(timestamp).txt
done

echo "===== BENCHMARK COMPLETE ====="
