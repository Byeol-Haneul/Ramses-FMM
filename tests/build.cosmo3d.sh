#!/bin/bash

cd bin
make NDIM=3 MPI=1 UNITS=cosmo CONDINIT=default HYDRO=1 GRAV=1 NVAR=7

