#!/bin/bash
#SBATCH --job-name=direct-phi
#SBATCH --output=direct_phi_%j.out
#SBATCH --error=direct_phi_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mem=512G
#SBATCH --ntasks-per-node=96
#SBATCH --exclusive


module load anaconda3/2024.2
module load openmpi/gcc/4.1.2
conda activate directphi

mpirun -np 16 python3 direct.py \
  --ic-file ../nfw/ic_part \
  --box-cut 300.0 \
  --lvl 9 \
  --z-cell 128 \
  --G 1.0 \
  --out-dir ./txt \
  --prefix direct_phi_slice
