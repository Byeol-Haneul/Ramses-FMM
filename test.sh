#!/bin/bash
#SBATCH --job-name=ramses-test
#SBATCH --output=/home/jl4415/mini-ramses/run_%j.out
#SBATCH --error=/home/jl4415/mini-ramses/run_%j.err
#SBATCH --qos=debug
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --mail-user=jl4415@princeton.edu
#SBATCH --mail-type=END,FAIL


# Optional: specify account if required
# #SBATCH --account=stellar
#!/bin/bash

module purge

# Run your program
module load openmpi/gcc/4.1.2
module load anaconda3/2025.6


srun -n 96 ramses3d namelist/fmm_tests/isolated_halo.nml
