!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(r,g,x,f,dx,ncell)
  use amr_parameters, only: ndim, nvector
  use amr_commons, only: run_t, global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  integer::ncell                              ! Size of input arrays
  real(kind=8)::dx                            ! Cell size
  real(kind=8),dimension(1:nvector,1:ndim)::f ! Gravitational acceleration
  real(kind=8),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  integer::idim,i
  real(kind=8)::gmass,emass,xmass,ymass,zmass,rr,rx,ry,rz

  ! Point mass
  gmass=1.0
  emass=0.0
  xmass=r%boxlen/2d0
  ymass=r%boxlen/2d0
  zmass=r%boxlen/2d0
  do i=1,ncell
     rx=0.0d0; ry=0.0d0; rz=0.0d0
     rx=x(i,1)-xmass
     ry=x(i,2)-ymass
     rz=x(i,3)-zmass
     rr=sqrt(rx**2+ry**2+rz**2+emass**2)
     f(i,1)=-gmass*rx/rr**3
     f(i,2)=-gmass*ry/rr**3
     f(i,3)=-gmass*rz/rr**3
  end do

end subroutine gravana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine phiana(r,g,x,phi,dx,ncell)
  use amr_parameters, only: ndim, nvector
  use amr_commons, only: run_t, global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  integer::ncell                              ! Size of input arrays
  real(kind=8)::dx                            ! Cell size
  real(kind=8),dimension(1:nvector)::phi      ! Gravitational potential
  real(kind=8),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the potential using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! phi(i is the gravitational potential in user units.
  !================================================================
  integer :: i
  real(kind=8)::fourpi,rx,ry,rz,rr

  fourpi=4.D0*ACOS(-1.0D0)

  do i=1,ncell
     rx=0.0d0; ry=0.0d0; rz=0.0d0
     rx=x(i,1)-g%multipole%q(2)/g%multipole%q(1)
#if NDIM>1
     ry=x(i,2)-g%multipole%q(3)/g%multipole%q(1)
#endif
#if NDIM>2
     rz=x(i,3)-g%multipole%q(4)/g%multipole%q(1)
#endif
     rr=MAX(sqrt(rx**2+ry**2+rz**2),dx)
#if NDIM==1
     phi(i)=g%multipole%q(1)*fourpi/2d0*rr
#endif
#if NDIM==2
     phi(i)=g%multipole%q(1)*2d0*log(rr)
#endif
#if NDIM==3
     phi(i)=-g%multipole%q(1)/rr
#endif
  end do

end subroutine phiana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
