The namelist block `&POISSON_PARAMS` is used to specify runtime parameters for the gravitational acceleration. The `&RUN_PARAMS` parameter `poisson=.true.` will activate self-gravity. If `poisson=.false.` then gravity is considered an external force with constant acceleration specified in the `&HYDRO_PARAMS`namelist.

Two different Poisson solvers are available in RAMSES: conjugate gradient (CG) and multigrid (MG). Unlike the CG solver, MG has an initialization overhead cost (at every call of the solver), but is much more efficient on very big levels with few "holes".  MG is always used for `levelmin`. MG can also be used on refined levels in conjuction with CG. The parameter `cg_levelmin` selects the Poisson solver as follows:

* The coarse level at `levelmin` is solved with MG
* Refined levels with `level<cg_levelmin` are solved with MG
* Refined levels with `level>=cg_levelmin` are solved with CG

| Variable name | Fortran type | Default value  | Description      |
|:------------------- |:-------|:----- |:------------------------- |
| `epsilon`           | `real`  | 1e-4  | Stopping criterion for the iterative Poisson solver: residual 2-norm should be lower than `epsilon` times the right hand side 2-norm. |
| `cg_levelmin`       | `integer`  | 999 | Minimum level from which the Conjugate Gradient solver is used in place of the Multigrid solver. |



