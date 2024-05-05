This sets of parameters, contained in the namelist block `&INIT_PARAMS`. This is used to set up the initial conditions.
 
| Variable name, syntax, default value | Fortran type | Description |
|:---------------------------- |:------------- |:------------------------- |
| `initfile=' ' ` | `80*char` | Directory where IC files are stored.
| `filetype='ascii'`  | `20*char` | Type of initial conditions file for particles and gas cells. Possible choices are `ascii`, `gadget`, `grafic` or `grafic_zoom`. Default value is `ascii`. In this case, only particles are read in the ascii file called `ic_part`, `ic_star`or `ic_sink`. For gas cells, initial conditions are set using file `condinit.f90` via compilation time parameters or via a patch. If nothing is done to set user defined initial conditions, then the default setup strategy is basd on uniform regions as explained now. |
| `nregion=1`  | `integer` | Number of independent regions in the computational box used to set up initial flow variables. |
| `region_type='square'`  | `10*char` | Geometry defining each region. `square` defines a generalized ellipsoidal shape, while `point` defines a delta function. |
| `x_center=0.0`  | `real arrays` | X coordinate of the center of each region. |
| `y_center=0.0`  | `real arrays` | Y coordinate of the center of each region. |
| `z_center=0.0`  | `real arrays` | Z coordinate of the center of each region. |
| `length_x=0.0`  | `real arrays` | Size in X direction of each region. Used only for `square` regions.|
| `length_y=0.0`  | `real arrays` | Size in Y direction of each region. Used only for `square` regions.|
| `length_z=0.0`  | `real arrays` | Size in Z direction of each region. Used only for `square` regions.|
| `exp_region=2.0`  | `real arrays` | Exponent defining the norm used to compute distances for the generalized ellipsoid. `exp_region=2` corresponds to a spheroid, `exp_region=1` to a diamond shape, `exp_region>=10` to a perfect square. |
| `d_region=0.0`  | `real arrays` | Density. For `point` regions this is used to define mass. |
| `u_region=0.0`  | `real arrays` | X velocity. For `point` regions this is used to define momentum. |
| `v_region=0.0`  | `real arrays` | Y velocity. For `point` regions this is used to define momentum. |
| `w_region=0.0`  | `real arrays` | Z velocity. For `point` regions this is used to define momentum. |
| `A_region=0.0`  | `real arrays` | X component of the magnetic field. |
| `B_region=0.0`  | `real arrays` | Y component of the magnetic field. |
| `C_region=0.0`  | `real arrays` | Z component of the magnetic field. |
| `p_region=0.0`  | `real arrays` | Pressure. For `point` regions this is used to define energy. |
| `prad_region=0.0` | `real arrays` | Non-thermal pressure. For `point` regions this is used to define energy. This is used only if compilation option `NENER>0`. |
| `var_region=0.0` | `real arrays` | Passive scalars. This is used only if compilation option `NVAR>5+NENER`. |
| `aexp_ini=10.0`  | `real` | This parameter sets the starting expansion factor for cosmology runs only. Default value is read in the IC file. |
| `omega_b=0.045`  | `real` | This parameter sets the baryonic density parameter for cosmology runs only. Default value is `omega_b=0.045`.  |


