This namelist is called &HYDRO_PARAMS, and is used to specify runtime parameters for the hydro and MHD solver. These parameters are based on rather standard concepts in computational fluid dynamics. We briefly describe them now.

| Variable name, syntax, default value | Fortran type  | Description               |
|:---------------------------- |:------------- |:------------------------- |
| `gamma=1.4`&nbsp;&nbsp;&nbsp;&nbsp;           |  `Real`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    | Adiabatic exponent for the perfect gas EOS |
| `gamma_rad=1.333`     |  `Real array` | Adiabatic exponent for each non-thermal pressure EOS (used only if `NENER>1`) |
| `courant_factor=0.5`  |  `Real`    | CFL number for time step control (less than 1) |
| `constant_gravity=0.0` |  `Real array` | 3-component of the constant gravitational acceleration in case `poisson=.false.` |
| `smallr=1d-10 `       |  `Real`    | Minimum density to prevent floating exceptions |
| `smallc=1d-10 `       |  `Real`    | Minimum sound speed to prevent floating exceptions |
| `riemann=’llf’`       |  `Character LEN=20`| Name of the Riemann solver. For the hydro solver (`MHD=0`), possible choices are `llf`, `hll` or `hllc`. For the MHD solver (`MHD=1`), possible choices are `llf`, `hll`, `roe`, `hlld` and `upwind`. |
| `riemann2d=’llf’`     |   `Characher LEN=20`| Name of the 2D Riemann solver for the induction equation (`MHD=1` only). Possible choices are `llf`, `hll`, `roe`,`hlla`, `hllf`, `hlld` and `upwind`. |
| `slope_type=1` |  `Integer`    | Type of slope limiter used for the piecewise linear reconstruction of volume-averaged quantities: `slope_type=0`: first order scheme, `slope_type=1`: MinMod limiter, `slope_type=2`: MonCen limiter, `slope_type=3`: Multi-dimensional MonCen limiter. In 1D runs only, it is also possible to choose: `slope_type=4`: Superbee limiter, `slope_type=5`: Ultrabee limiter. |
| `slope_mag_type=1` |  `Integer`    | Type of slope limiter used for the piecewise linear reconstruction of face-averaged quantities (`MHD=1` only). |
| `turb=.false.` |  `Logical`    | Solve the turbulent kinetic energy equation using a standard LES model. |
| `induction=.false.` |  `Logical`    | Limit the MHD solver to just the induction equation. The velocity field is set in the induction equation but none of the hydro variables are updated. |
| `entropy=.false.` |  `Logical`    | Solve for the conservation of entropy as an additional passive scalar. You must compile the code with `NVAR>5`. It can also be used in conjunction with the `dual_energy` parameter. |
| `dual_energy=-1` |  `Real`    | Activate dual energy scheme for high-Mach flows if `dual_energy>=0`. Useful to prevent negative temperatures. The chosen value is used to set the fraction of the energy truncation error. It must be used in conjunction with ``entropy=.true.`. |
| `difmag=0d0` |  `Real`    | Add explicit diffusion for all volume-averaged conservative variable. | 
| `etamag=0d0` |  `Real`    | Add explicit magnetic diffusivity (Ohm's law). |
