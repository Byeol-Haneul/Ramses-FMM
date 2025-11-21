#!/bin/bash
#SBATCH --job-name=ramses-bench
#SBATCH --output=/home/jl4415/mini-ramses/bench_%j.out
#SBATCH --error=/home/jl4415/mini-ramses/bench_%j.err
#SBATCH --qos=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=3
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
MAX_RANKS=$(( CORES_PER_NODE * NODES_ALLOCATED ))   # = 192
OUTDIR="/home/jl4415/mini-ramses/complete_benchmark"
mkdir -p $OUTDIR

timestamp() { date +"%Y%m%d_%H%M%S"; }

# ============================
# Rank Generators
# ============================

# Existing strong-scaling generator (powers of two, starting from 16)
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
        local next=$((val * 2))
        if (( next > max )); then
            ranks+=("$max")
            break
        fi
        val=$next
    done
    echo "${ranks[@]}"
}

# Weak-scaling generator (8 → 32 → 128 → 512, capped at user max)
generate_weak_ranks() {
    local max=$1
    local ranks=(8)
    local val=8
    while :; do
        val=$(( val * 4 ))
        if (( val >= max )); then
            ranks+=("$max")
            break
        fi
        ranks+=("$val")
    done
    echo "${ranks[@]}"
}

STRONG_RANKS=($(generate_ranks $MAX_RANKS))
WEAK_RANKS=($(generate_weak_ranks 512))   # up to 512 ranks available

# ============================
# WEAK SCALING (Levels 8 → 11)
# ============================
'''
echo "===== WEAK SCALING | Levels 8 → 11 ====="

WEAK_LEVELS=(8 9 10 11)

for i in "${!WEAK_RANKS[@]}"; do
    NP=${WEAK_RANKS[$i]}
    LVL=${WEAK_LEVELS[$i]}

    echo "--- Weak scaling: MPI ranks $NP | Level $LVL ---"

    srun  -n ${NP} ./ramses_mg \
         namelist/benchmark/lvl${LVL}_fmm1.nml
    mv time_mg.txt  ${OUTDIR}/weak_mg_lvl${LVL}_P${NP}_$(timestamp).txt

    srun  -n ${NP} ./ramses_fmm \
         namelist/benchmark/lvl${LVL}_fmm1.nml
    mv time_fmm.txt ${OUTDIR}/weak_fmm1_lvl${LVL}_P${NP}_$(timestamp).txt

    srun  -n ${NP} ./ramses_fmm \
         namelist/benchmark/lvl${LVL}_fmm2.nml
    mv time_fmm.txt ${OUTDIR}/weak_fmm2_lvl${LVL}_P${NP}_$(timestamp).txt
done
'''
# ============================
# STRONG SCALING (fixed level 9)
# ============================

STRONG_LVL=10
echo "===== STRONG SCALING | Level $STRONG_LVL ====="

'''
for NP in "${STRONG_RANKS[@]}"; do
    if (( NP > MAX_RANKS )); then
        echo "Skipping $NP (exceeds $MAX_RANKS)"
        continue
    fi

    echo "--- Strong scaling: MPI ranks $NP ---"

    srun  -n ${NP} ./ramses_mg \
         namelist/benchmark/lvl${STRONG_LVL}_fmm1.nml
    mv time_mg.txt  ${OUTDIR}/strong_mg_lvl${STRONG_LVL}_P${NP}_$(timestamp).txt

    srun  -n ${NP} ./ramses_fmm \
         namelist/benchmark/lvl${STRONG_LVL}_fmm1.nml
    mv time_fmm.txt ${OUTDIR}/strong_fmm1_lvl${STRONG_LVL}_P${NP}_$(timestamp).txt

    srun  -n ${NP} ./ramses_fmm \
         namelist/benchmark/lvl${STRONG_LVL}_fmm2.nml
    mv time_fmm.txt ${OUTDIR}/strong_fmm2_lvl${STRONG_LVL}_P${NP}_$(timestamp).txt
done
'''

NP=128
STRONG_LVL=10
srun  -n ${NP} ./ramses_mg \
      namelist/benchmark/lvl${STRONG_LVL}_fmm1.nml
mv time_mg.txt  ${OUTDIR}/strong_mg_lvl${STRONG_LVL}_P${NP}_$(timestamp).txt

srun  -n ${NP} ./ramses_fmm \
      namelist/benchmark/lvl${STRONG_LVL}_fmm1.nml
mv time_fmm.txt ${OUTDIR}/strong_fmm1_lvl${STRONG_LVL}_P${NP}_$(timestamp).txt

srun  -n ${NP} ./ramses_fmm \
      namelist/benchmark/lvl${STRONG_LVL}_fmm2.nml
mv time_fmm.txt ${OUTDIR}/strong_fmm2_lvl${STRONG_LVL}_P${NP}_$(timestamp).txt
echo "===== BENCHMARK COMPLETE ====="
