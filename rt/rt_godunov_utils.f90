!###########################################################
!###########################################################
!###########################################################
subroutine get_rt_courant_dt(r,g,dt,ilevel)
!-------------------------------------------------------------------------
! Determine the coarse RT timestep length set by the Courant condition
!-------------------------------------------------------------------------
  use amr_parameters
  use rt_parameters
  use amr_commons, only: run_t, global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  integer::ilevel
  real(kind=8)::dt,dx
!-------------------------------------------------------------------------
  ! Mesh spacing at coarse level
  dx = 0.5D0**ilevel*r%boxlen
  dt = r%rt_courant_factor*dx/3d0/g%rt_c(ilevel)
end subroutine get_rt_courant_dt

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine rt_refine(r,ug,um,ud,ok)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t
  use rt_parameters, only: nrtvar, nrtgrp
  implicit none
  ! routine arguments
  type(run_t)::r
  real(kind=8)::ug(1:nrtvar)
  real(kind=8)::um(1:nrtvar)
  real(kind=8)::ud(1:nrtvar)
  logical ::ok
  ! local variables
  integer::igroup, ivar
  real(kind=8)::ng,nm,nd,error

  ! Compute errors
  do igroup=1, nrtgrp
    if(r%rt_err_grad_cn(igroup) >= 0.)then
      ivar = 1+(igroup-1)*(ndim+1)
      ng=ug(ivar); nm=um(ivar); nd=ud(ivar)
      error=2.0d0*MAX( &
          & ABS((nd-nm)/(nd+nm+r%rt_floor_cn(igroup))) , &
          & ABS((nm-ng)/(nm+ng+r%rt_floor_cn(igroup))) )
      ok = ok .or. error > r%rt_err_grad_cn(igroup)
    end if
  end do

end subroutine rt_refine
