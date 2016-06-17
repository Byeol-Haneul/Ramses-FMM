#! /bin/bash

cp  ../../amr/update_time.f90 ./
cp  ../../bin/Makefile ./

ending=`cat Makefile | grep NDIM -m 1 | sed 's/[^0-9]//g'`d

echo $ending
cp halo_ics_small/ic_part$ending halo_ics_small/ic_part

patch -R update_time.f90 < update_time.f90_diff
patch -R Makefile < Makefile_diff 

cd ../../bin/
make clean 
rm ramses$ending
make -f ../patch/integration_tests/Makefile
cd ../patch/integration_tests/

../../bin/ramses$ending param_file.nml > run_serial.log
if diff particles.txt particles.txt_comp_$ending > /dev/null; then
    echo "=============================="
    echo "NON MPI TEST OK"
    echo "=============================="
    rm particles.txt
else
    echo "=============================="
    echo "NON MPI TEST FAILED!!!!!!!!!!"
    echo "=============================="
    mv particles.txt particles_serial.txt
fi
mpirun -np 3 ../../bin/ramses$ending param_file.nml > run_parallel.log
if diff particles.txt particles.txt_comp_$ending > /dev/null; then
    echo "=============================="
    echo "MPI TEST OK"
    echo "=============================="
    rm particles.txt
    rm update_time.f90
    rm Makefile
else
    echo "=============================="
    echo "MPI TEST FAILED!!!!!!!!!!!!!!"
    echo "=============================="
    mv particles.txt particles_parallel.txt
fi
