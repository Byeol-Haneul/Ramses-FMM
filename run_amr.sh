#!/bin/bash
rm -rf output_00002; cd bin; make; cd ..; bin/ramses3d namelist/spheres.nml; mv output_00002 output_00013
rm -rf output_00002; ./ramses_mg_dump namelist/spheres.nml; mv output_00002 output_00014