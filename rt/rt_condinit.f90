#ifdef RT
!================================================================
!================================================================
!================================================================
!================================================================
subroutine rt_condinit(r,g,x,q,dx,nn)
  use amr_parameters, only: dp, ndim, nvector
  use rt_parameters, only: nrtvar
  use amr_commons, only: run_t, global_t
  use input_rt_condinit_module, only: rt_region_condinit
  implicit none
  type(run_t)::r
  type(global_t)::g
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nrtvar)::q ! RT variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user (aka code) units:
  ! x(i,1:ndim) are in [0,boxlen]**ndim.
  ! Q is the RT variable vector in code units.
  !================================================================
  ! Call built-in initial condition generator
  call rt_region_condinit(r,g,x,q,dx,nn)
  
  ! Add here, if you wish, some user-defined initial conditions
  ! ........


end subroutine rt_condinit
#endif