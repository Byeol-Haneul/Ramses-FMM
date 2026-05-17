!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine rho_ana(x,d,dx,gravity_params)
  use amr_parameters, only: ndim
  implicit none
  real(kind=8),dimension(1:10)::gravity_params
  real(kind=8)::dx                  ! Cell size
  real(kind=8)::d                   ! Density
  real(kind=8),dimension(1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates analytical Poisson source term.
  ! Positions are in user units:
  ! x(1:ndim) are in [0,box_size]**ndim.
  ! d is the density field in user units.
  !================================================================
  real(kind=8)::dmass,emass,xmass,ymass,zmass,rr,rx,ry,rz,dd

  emass=gravity_params(1) ! Softening length
  xmass=gravity_params(2) ! Point mass coordinates
  ymass=gravity_params(3)
  zmass=gravity_params(4)
  dmass=1.0/(emass*(1.0+emass)**2)

  rx=0.0d0; ry=0.0d0; rz=0.0d0
  rx=x(1)-xmass
#if NDIM>1
  ry=x(2)-ymass
#endif
#if NDIM>2
  rz=x(3)-zmass
#endif
  rr=sqrt(rx**2+ry**2+rz**2)
  dd=1.0/(rr*(1.0+rr)**2)
  d=MIN(dd,dmass)

end subroutine rho_ana
