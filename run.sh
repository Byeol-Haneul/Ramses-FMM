#!/bin/bash
#SBATCH --job-name=stellar-debug
#SBATCH --output=logs/stellar-debug_%j.out
#SBATCH --error=logs/stellar-debug_%j.err
#SBATCH --time=00:30:00              # Max allowed for stellar-debug
#SBATCH --qos=stellar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@princeton.edu

# Optional: specify account if required
# #SBATCH --account=stellar
#!/bin/bash

while true; do
    echo "This will run forever..."
    sleep 1  # optional, to avoid spamming the terminal
done

module purge

# Run your program
module load openmpi/gcc/4.1.2

#cd bin ;  make clean; make NVAR=6 FMM=0; cd .. ; bin/ramses3d namelist/fmm_tests/isolated_offset.nml
cd bin  ;  make NVAR=6 FMM=1 MPI=1 HYDRO=0; cd .. ; mpirun -np 2 bin/ramses3d namelist/isolated_halo.nml
#cd bin ;  make NVAR=6 FMM=1 MPI=0; cd .. ; bin/ramses3d namelist/halo.nml
#cd bin ;  make NVAR=6 FMM=1 MPI=0; cd .. ; bin/ramses3d namelist/single.nml


