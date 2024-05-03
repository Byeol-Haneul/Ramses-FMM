We discuss here the various compilation time parameters. These parameters are used with the Makefile available in directory `bin/`.
You can compile the code by setting these parameters to your preffered value using:
`make PARAM=value`. You can see many examples in the directory `namelist/`.

| Variable name, syntax, default value | Type | Description |
|:-------- |:---- |:----------------- |
| `COMPILER = GNU` | `GNU` or `INTEL` | This sets the Fortran compiler you would like to use. Only 2 options are available, namely `gfortran` from the Gnu Compiler Collection and `ifort` from Intel, Inc. |
| `EXEC = ramses`  | `string` | This sets the name of the executable for the compiled code. The number of dimensions used at compilation via parameter `NDIM` is appended at the end. |
| `NHILBERT = 1`   | `1`, `2` or `3`  | This sets the number of long integers uswed to code the Hilbert key. This choice sets also the maximum level of refinement. In 3D, `NHILBERT=1` corresponds to `levelmax<21`, `NHILBERT=2` to `levelmax<42` and `NHILBERT=3` to `levelmax<63`.|
| `NVECTOR = 32`   | `integer`  | This sets the size of the vector sweeps used in many subroutines of the code. This is used for optimization purposes. This is highly problem and architecture dependant. |
| `NPRE = 8`       | `4` or `8`  | This sets the number of bytes used to code floating point numbers. `4` means single precision (not recommanded) `8` means double precision. |
| `NDIM = 3`       | `1`, `2` or `3`  | This sets the dimensionality of the simulation. |
| `NVAR = 5`       | `integer`  | Number of hydro variables used in the hydro solver. The default (and minimum) value is `5`. You can use more variables for passive scalars. |
| `NENER = 0`      | `integer`  | Number of non-thermal energies used in the hydro solver. The default value is `0`. These are different than passive scalars. |
| `HYDRO = 1`      | `0` or `1`  | This turns on or off the hydro solver. |
| `MHD = 0`        | `0` or `1`  | This turns on or off the MHD solver. |
| `GRAV = 0`       | `0` or `1`  | This turns on or off the self-gravity solver. |
| `MPI = 0`        | `0` or `1`  | This turns on or off parallel computing using the Message Passing Interface (MPI) library. MPI must be properly installed on your machine. In particular, you have to check that the MPI library was configured for your Fortran compiler. |
| `MPIF90 = mpif90` | `string`  | This sets the name of the MPI Fortran compiler on your system. |
| `INIT = `        | `string`  | This sets the adopted initial condition for your simulation. This parameter can be an empty string (default) wich uses the default initial condition generator in the code, as explained in the corresponding namelist block documentation. Possible values are `COEUR`, `INSTA`, `DOUBLEMACH`, `OT`, `PONO` and more. These initial conditions are set in the file `hydro/condinit.f90`. You can add your own here and share them with the community later via a pull request. |
| `UNITS = `       | `string`  | This sets the adopted unit system to convert variables from code units to cgs units. This parameter can be an empty string (default) wich uses the default unit system that can be set via the namelist, as explained in the corresponding namelist block documentation. You can also use predefined unit systems. Possible values are `COEUR`, `COSMO`, `MERGER` and more. These different unit systems are set in the file `amr/units.f90`. You can add your own here and share them with the community later via a pull request. |
| `DEBUG = 0`      | `0` or `1`  | This turns on or off the debug mode. |





