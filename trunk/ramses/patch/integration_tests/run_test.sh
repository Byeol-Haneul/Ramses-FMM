#! /bin/bash

export OMPI_FC=ifort

cp  ../../amr/update_time.f90 ./
patch -R update_time.f90 < update_time.f90_diff 

cd ../../bin/
make clean > log.txt
make -f ../patch/integration_tests/Makefile > log.txt
cd ../patch/integration_tests/

../../bin/ramses3d param_file.nml > log.txt
if diff particles.txt particles.txt > /dev/null; then
    echo "=============================="
    echo "NON MPI TEST OK"
    echo "=============================="
else
    echo "=============================="
    echo "NON MPI TEST FAILED!!!!!!!!!!"
    echo "=============================="
fi
mpirun -np 3 ../../bin/ramses3d param_file.nml > log.txt
if diff particles.txt particles.txt > /dev/null; then
    echo "=============================="
    echo "MPI TEST OK"
    echo "=============================="
else
    echo "=============================="
    echo "MPI TEST FAILED!!!!!!!!!!!!!!"
    echo "=============================="
fi
rm particles.txt
rm update_time.f90
