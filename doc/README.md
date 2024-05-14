
## Documentation ##

This directory contains the basic documentation for mini-ramses.

We describe the meaning and default values of the namelist variables in each of the following namelist blocks:

[Compilation parameters](./compilation_params.md): parameters used to compile the code using the `make` command.

[Global run parameters](./run_params.md): Main global parameters of the code controlling its execution. 

[Output parameters](./output_params.md): Parameters controlling the output strategy.

[AMR grid parameters](./amr_params.md): Parameters controlling mesh generation and memory.

[Poisson solver parameters](./poisson_params.md): Parameters controlling the Poisson solver for self-gravity.

[Initial conditions parameters](./init_params.md): Parameters used to setup the initial conditions.

[Hydro and MHD solver parameters](./hydro_params.md): Parameters controlling the second-order Godunov solver (aka MUSCL-Hancock) for solving the Euler and ideal MHD equations.

[Mesh refinement parameters](./refine_params.md): Parameters controlling the refinement strategy of the AMR grid. 

[Unit system parameters](./units_params.md): Parameters controlling the conversion between code units and cgs units.

[Boundary conditions parameters](./boundary_params.md): Parameters controlling the boundary conditions of the hydro solver.

[Cooling parameters](./cooling_params.md): Parameters controlling the cooling function used in the code.

[Star formation parameters](./star_params.md): Parameters controlling the subgrid star formation model used in the code.

[Clump finder parameters](./clump_params.md): Parameters controlling the clump finder.








