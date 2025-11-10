#!/bin/bash
#cd bin ;  make NVAR=6 FMM=1; cd .. ; bin/ramses3d namelist/isolated.nml
#cd bin ;  make clean; make NVAR=6 FMM=0; cd .. ; bin/ramses3d namelist/isolated.nml
#cd bin  ; make NVAR=6 FMM=1 MPI=1; cd .. ; mpirun -np 2 bin/ramses3d namelist/single.nml
#cd bin ; make clean; make NVAR=6 FMM=0; cd .. ; bin/ramses3d namelist/eight.nml
cd bin ;  make NVAR=6 FMM=1 MPI=0; cd .. ; bin/ramses3d namelist/halo.nml
#cd bin ; make clean; make NVAR=6 FMM=0; cd .. ; bin/ramses3d namelist/isolated_halo.nml

