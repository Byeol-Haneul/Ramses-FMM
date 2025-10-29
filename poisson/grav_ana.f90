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

  ! Multipole expansion for isolated boundary conditions
  if(r%gravity_type==0)then
     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        rx=x(i,1)-g%multipole%q(2)/g%multipole%q(1)
#if NDIM>1
        ry=x(i,2)-g%multipole%q(3)/g%multipole%q(1)
#endif
#if NDIM>2
        rz=x(i,3)-g%multipole%q(4)/g%multipole%q(1)
#endif
        rr=sqrt(rx**2+ry**2+rz**2)
#if NDIM==1
        f(i,1)=-g%multipole%q(1)*2d0*ACOS(-1d0)*rx/rr
#endif
#if NDIM==2
        f(i,1)=-g%multipole%q(1)*2d0*rx/rr*2
        f(i,2)=-g%multipole%q(1)*2d0*ry/rr*2
#endif
#if NDIM==3
        f(i,1)=-g%multipole%q(1)*rx/rr*3
        f(i,2)=-g%multipole%q(1)*ry/rr*3
        f(i,3)=-g%multipole%q(1)*rz/rr*3
#endif
     end do
  end if

  ! Constant vector
  if(r%gravity_type==1)then
     do idim=1,ndim
        do i=1,ncell
           f(i,idim)=r%gravity_params(idim)
        end do
     end do
  end if

  ! Point mass
  if(r%gravity_type==2)then
     gmass=r%gravity_params(1) ! GM
     emass=r%gravity_params(2) ! Softening length
     xmass=r%gravity_params(3) ! Point mass coordinates
     ymass=r%gravity_params(4)
     zmass=r%gravity_params(5)
     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        rx=x(i,1)-xmass
#if NDIM>1
        ry=x(i,2)-ymass
#endif
#if NDIM>2
        rz=x(i,3)-zmass
#endif
        rr=sqrt(rx**2+ry**2+rz**2+emass**2)
#if NDIM==1
        f(i,1)=-gmass*rx/rr
#endif
#if NDIM==2
        f(i,1)=-gmass*ry/rr**2
        f(i,2)=-gmass*ry/rr**2
#endif
#if NDIM==3
        f(i,1)=-gmass*rx/rr**3
        f(i,2)=-gmass*ry/rr**3
        f(i,3)=-gmass*rz/rr**3
#endif
     end do
  end if

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
     rr=sqrt(rx**2+ry**2+rz**2)
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
