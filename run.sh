#!/bin/bash
#SBATCH --job-name=ramses-test
#SBATCH --output=/home/jl4415/mini-ramses/run_%j.out
#SBATCH --error=/home/jl4415/mini-ramses/run_%j.err
#SBATCH --qos=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mail-user=jl4415@princeton.edu
#SBATCH --mail-type=END,FAIL


# Optional: specify account if required
# #SBATCH --account=stellar
#!/bin/bash

module purge

# Run your program
module load openmpi/gcc/4.1.2

srun -n 64 bin/ramses3d namelist/benchmark/lvl8_fmm1.nml
#srun -n 16 ramses_mg_dump  namelist/benchmark/lvl8_fmm1.nml

