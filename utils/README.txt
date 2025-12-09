This directory contains useful tools related to mini-ramses, mostly for preparing ICs and postprocessing.
Add tools from the original RAMSES utils here once they have been adapted to the new RAMSES data structure.

part2cube (mini-ramses)
- Path: utils/f90/part2cube.f90
- Build (GNU):
  gfortran -O3 -ffree-line-length-none -o part2cube utils/f90/part2cube.f90
- Usage:
  ./part2cube -inp /abs/path/to/output_00001 -pre part \
              -dep CIC|TSC|PCS -nx 256 -ny 256 -nz 256 -per F \
              -xmi 0 -xma 1 -ymi 0 -yma 1 -zmi 0 -zma 1
- Output: writes prefix.cube (unformatted), with dims then real*4 cube data
