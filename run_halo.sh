#!/bin/bash
#SBATCH --job-name=ramses-amr-halo
#SBATCH --output=/home/jl4415/mini-ramses/halo_%j.out
#SBATCH --error=/home/jl4415/mini-ramses/halo_%j.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --mail-user=jl4415@princeton.edu
#SBATCH --mail-type=END,FAIL

module purge
module load openmpi/gcc/4.1.2

cd /home/jl4415/mini-ramses || exit 1

export SLURM_CPU_BIND=cores
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp

NP=96
NML_FILE="/home/jl4415/mini-ramses/namelist/amr_tests/isolated_halo.nml"

# ============================
# MG run
# ============================
cd /home/jl4415/mini-ramses/halo_test/halo_mg || exit 1
srun -n "$NP" ./mg_halo "$NML_FILE" \
    > mg_${SLURM_JOB_ID}.out \
    2> mg_${SLURM_JOB_ID}.err

# ============================
# FMM run
# ============================
cd /home/jl4415/mini-ramses/halo_test/halo_fmm || exit 1
srun -n "$NP" ./fmm_halo "$NML_FILE" \
    > fmm_${SLURM_JOB_ID}.out \
    2> fmm_${SLURM_JOB_ID}.err