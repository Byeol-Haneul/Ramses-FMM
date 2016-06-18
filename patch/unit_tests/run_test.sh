#! /bin/bash


cp  ../../bin/Makefile ./
ending=`cat Makefile | grep NDIM -m 1 | sed 's/[^0-9]//g'`d

patch -R Makefile < Makefile_diff 

cd ../../bin/
make clean 
make -f ../patch/unit_tests/Makefile 
cd ../patch/unit_tests/

../../bin/ramses_unit_tests$ending param_file.nml

#rm Makefile
