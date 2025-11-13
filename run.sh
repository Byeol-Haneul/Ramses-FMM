#!/bin/bash
#cd bin ; make clean; make NVAR=6 FMM=1 MPI=0; cd .. ; bin/ramses3d namelist/fmm_tests/single1.nml
#bin/ramses3d namelist/fmm_tests/single2.nml
#cd bin ; make NVAR=6 FMM=1 MPI=1; cd .. ; mpirun -np 2 bin/ramses3d namelist/fmm_tests/halo1.nml
#mpirun -np 2 bin/ramses3d namelist/fmm_tests/halo2.nml
#bin/ramses3d namelist/fmm_tests/single_offset1.nml
#bin/ramses3d namelist/fmm_tests/single_offset2.nml


cd bin ;  make clean; make NVAR=6 FMM=0 MPI=1; cd .. ; mpirun -np 2 bin/ramses3d namelist/fmm_tests/isolated_offset.nml
#cd bin  ;  make clean; make NVAR=6 FMM=0 MPI=1; cd .. ; mpirun -np 2 bin/ramses3d namelist/isolated_halo.nml
#cd bin ;  make NVAR=6 FMM=1 MPI=0; cd .. ; bin/ramses3d namelist/halo.nml
#cd bin ;  make NVAR=6 FMM=1 MPI=0; cd .. ; bin/ramses3d namelist/single.nml


