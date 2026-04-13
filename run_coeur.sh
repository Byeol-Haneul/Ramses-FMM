#!/bin/bash
#SBATCH --job-name=ramses-amr-coeur
#SBATCH --output=/home/jl4415/mini-ramses/coeur_%j.out
#SBATCH --error=/home/jl4415/mini-ramses/coeur_%j.err
#SBATCH --time=24:00:00
#SBATCH --qos=pu-short-stellar
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
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

NP=96
NML_FILE="/home/jl4415/mini-ramses/namelist/amr_tests/coeur.nml"

'''
cd /home/jl4415/mini-ramses/coeur_fmm
srun -n "$NP" ./fmm_coeur "$NML_FILE"
'''

cd /home/jl4415/mini-ramses/coeur_mg
srun -n "$NP" ./mg_coeur "$NML_FILE" \
    > mg_${SLURM_JOB_ID}.out \
    2> mg_${SLURM_JOB_ID}.err
