#!/bin/bash
#SBATCH --job-name=ramses-bench
#SBATCH --output=/home/jl4415/mini-ramses/bench_%j.out
#SBATCH --error=/home/jl4415/mini-ramses/bench_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=2
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

export SLURM_CPU_BIND=cores
export OMP_PLACES=cores
export OMP_PROC_BIND=close

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp

CORES_PER_NODE=96
NODES_ALLOCATED=2
MAX_RANKS=$(( CORES_PER_NODE * NODES_ALLOCATED ))

timestamp() { date +"%Y%m%d_%H%M%S"; }

# ============================
# Rank generators
# ============================

generate_ranks() {
    local max=$1
    local ranks=()
    local val=8
    while :; do
        if (( val <= max )); then
            ranks+=("$val")
        else
            ranks+=("$max")
            break
        fi
        local next=$((val * 2))
        if (( next > max )); then
            ranks+=("$max")
            break
        fi
        val=$next
    done
    echo "${ranks[@]}"
}

generate_weak_ranks() {
    local max=$1
    local ranks=(1)
    local val=1
    while :; do
        val=$(( val * 8 ))
        if (( val >= max )); then
            ranks+=("$max")
            break
        fi
        ranks+=("$val")
    done
    echo "${ranks[@]}"
}

STRONG_RANKS=($(generate_ranks $MAX_RANKS))
WEAK_RANKS=($(generate_weak_ranks $MAX_RANKS))
WEAK_LEVELS=(8 9 10)

# ============================
# OUTER LOOP — RUN 5 TIMES
# ============================
for rep in {1}; do
    RUNSTAMP=$(timestamp)
    OUTDIR="/home/jl4415/mini-ramses/benchmark_stable/run_${RUNSTAMP}"
    #mkdir -p $OUTDIR

    echo "=============================================="
    echo "   BENCHMARK RUN $rep    ($RUNSTAMP)"
    echo "=============================================="
    # -------------------------------------------
    # STRONG SCALING: LEVEL 10
    # -------------------------------------------
    STRONG_LVL=10
    echo "===== STRONG SCALING | Level $STRONG_LVL ====="

    for NP in "${STRONG_RANKS[@]}"; do
        [[ $NP -gt $MAX_RANKS ]] && continue

        echo "--- Strong scaling: MPI ranks $NP ---"

        srun -n ${NP} ./ramses_fmm namelist/benchmark/lvl${STRONG_LVL}_fmm2.nml
        mv time_fmm.txt ${OUTDIR}/strong_fmm2_lvl${STRONG_LVL}_P${NP}_${RUNSTAMP}.txt

        srun -n ${NP} ./ramses_fmm namelist/benchmark/lvl${STRONG_LVL}_fmm1.nml
        mv time_fmm.txt ${OUTDIR}/strong_fmm1_lvl${STRONG_LVL}_P${NP}_${RUNSTAMP}.txt

        srun -n ${NP} ./ramses_mg namelist/benchmark/lvl${STRONG_LVL}_fmm1.nml
        mv time_mg.txt  ${OUTDIR}/strong_mg_lvl${STRONG_LVL}_P${NP}_${RUNSTAMP}.txt
    done

    echo "===== RUN $rep COMPLETE ====="
done

echo "===== ALL 5 BENCHMARK RUNS COMPLETE ====="
