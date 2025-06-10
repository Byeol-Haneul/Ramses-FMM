!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(r,g,x,u,dx,ibound,ncell)
  use amr_parameters, only: ndim, nvector
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t, global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  integer ::ibound                            ! Index of boundary region
  integer ::ncell                             ! Number of active cells
  real(kind=8)::dx                            ! Cell size
  real(kind=8),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(kind=8),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user (aka code) units:
  ! x(i,1:ndim) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w and U(i,5): E.
  ! U is in user (aka code) units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================
  integer::ivar,i

  ! Add here, if you wish, some user-defined boudary conditions
  ! ........

end subroutine boundana
