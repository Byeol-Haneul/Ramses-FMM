This sets of parameters, contained in the namelist block `&AMR_PARAMS`, controls the AMR grid global properties. Parameters specifying the refinement strategy are described in the namelist block `&REFINE_PARAMS` and are used only if `levelmax>levelmin`.
 
| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `levelmin=1`                 |  `integer`    | Minimum level of refinement. This parameter sets the size of the coarse (or base) grid to `nx=2**levelmin`.|
| `levelmax=1`                 |  `integer`    | Maximum level of refinement. If `levelmax=levelmin`, the simulation will be executed on a standad Cartesian grid of linear size `nx=2**levelmin`|
| `ngridmax=0`                 |  `integer`    | Maximum number of grids (or octs) that can be allocated during the run within each MPI process. |
| `npartmax=0`                 |  `integer`    | Maximum number of dark matter particles that can be allocated during the run within each MPI process. |
| `nexpand=1`                  |  `integer`    | Number of times the mesh expansion is applied to the refinement map (see mesh smoothing).|
| `boxlen=1.0`                 |  `real`       | Logical box size in code units. It corresponds to the square box in which Cartesian indices as well as Hilbert indices are defined.|
| `box_size=0.0`               |  `real`       | Physical length of the domain in the x-direction in code units. This is used only for non-periodic boundary conditions. |
| `box_xmin=0`                 |  `integer`    | Cartesian integer x-coordinate of the left boundary of the physical domain. This integer coordinate is defined for octs at levelmin. This is used only for non-periodic boundary conditions. |
| `box_xmax=0`                 |  `integer`    | Cartesian integer x-coordinate of the right boundary of the physical domain. This integer coordinate is defined for octs at levelmin. This is used only for non-periodic boundary conditions. |
| `box_ymin=0`                 |  `integer`    | Cartesian integer y-coordinate of the bottom boundary of the physical domain. This integer coordinate is defined for octs at levelmin. This is used only for non-periodic boundary conditions. |
| `box_ymax=0`                 |  `integer`    | Cartesian integer y-coordinate of the top boundary of the physical domain. This integer coordinate is defined for octs at levelmin. This is used only for non-periodic boundary conditions. |
| `box_zmin=0`                 |  `integer`    | Cartesian integer z-coordinate of the back boundary of the physical domain. This integer coordinate is defined for octs at levelmin. This is used only for non-periodic boundary conditions. |
| `box_zmax=0`                 |  `integer`    | Cartesian integer z-coordinate of the front boundary of the physical domain. This integer coordinate is defined for octs at levelmin. This is used only for non-periodic boundary conditions. |

