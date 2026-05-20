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
  ! x(i,1:ndim) are cell center position in [0,box_size] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  integer::idim,i,j
  real(kind=8)::gmass,emass,xmass,ymass,zmass,rr,rx,ry,rz
  real(kind=8)::mass, rr2, rr3, rr5, rr7, trace_c, scalar_q
  real(kind=8),dimension(3)::center, rr_vec, cr
  real(kind=8),dimension(3,3)::q_central
  
  ! Default to zero acceleration for unsupported gravity_type values.
  f(1:ncell,1:ndim)=0.0d0

  ! Multipole expansion for isolated boundary conditions
  if(r%gravity_type==0)then
     if (g%multipole%q(1) <= 0.0d0) return
     mass = g%multipole%q(1)
     center = 0.0d0
     center(1:ndim) = g%multipole%q(2:ndim+1)/mass
#ifdef FMM
#if NDIM==3
     q_central = 0.0d0
     q_central(1,1) = g%multipole%q(5)
     q_central(1,2) = g%multipole%q(6)
     q_central(2,1) = q_central(1,2)
     q_central(1,3) = g%multipole%q(7)
     q_central(3,1) = q_central(1,3)
     q_central(2,2) = g%multipole%q(8)
     q_central(2,3) = g%multipole%q(9)
     q_central(3,2) = q_central(2,3)
     q_central(3,3) = g%multipole%q(10)
     trace_c = q_central(1,1) + q_central(2,2) + q_central(3,3)
#endif
#endif
     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        rx=x(i,1)-center(1)
#if NDIM>1
        ry=x(i,2)-center(2)
#endif
#if NDIM>2
        rz=x(i,3)-center(3)
#endif
        rr=MAX(sqrt(rx**2+ry**2+rz**2),dx)
#if NDIM==1
        f(i,1)=-mass*2d0*ACOS(-1d0)*rx/rr
#endif
#if NDIM==2
        f(i,1)=-mass*2d0*rx/rr*2
        f(i,2)=-mass*2d0*ry/rr*2
#endif
#if NDIM==3
#ifdef FMM
        rr_vec = 0.0d0
        rr_vec(1) = rx
        rr_vec(2) = ry
        rr_vec(3) = rz
        rr2 = rr*rr
        rr3 = rr2*rr
        rr5 = rr3*rr2
        rr7 = rr5*rr2
        cr = 0.0d0
        do idim=1,3
           do j=1,3
              cr(idim) = cr(idim) + q_central(idim,j)*rr_vec(j)
           end do
        end do
        scalar_q = 0.0d0
        do idim=1,3
           do j=1,3
              scalar_q = scalar_q + rr_vec(idim)*q_central(idim,j)*rr_vec(j)
           end do
        end do
        f(i,1:3) = -mass*rr_vec(1:3)/rr3 + (3.0d0*cr(1:3) - trace_c*rr_vec(1:3))/rr5 - &
             & 2.5d0*(3.0d0*scalar_q - trace_c*rr2)*rr_vec(1:3)/rr7
#else
        f(i,1)=-mass*rx/rr**3
        f(i,2)=-mass*ry/rr**3
        f(i,3)=-mass*rz/rr**3
#endif
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
     xmass=r%gravity_params(3) ! Point mass x-coordinate
     ymass=r%gravity_params(4) ! Point mass y-coordinate
     zmass=r%gravity_params(5) ! Point mass z-coordinate
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
  ! x(i,1:ndim) are cell center position in [0,box_size] (user units).
  ! phi(i is the gravitational potential in user units.
  !================================================================
  integer :: i,idim,j
  real(kind=8)::fourpi,rx,ry,rz,rr
  real(kind=8)::mass, rr2, rr5, trace_c, scalar_q
  real(kind=8),dimension(3)::center, rr_vec
  real(kind=8),dimension(3,3)::q_central

  fourpi=4.D0*ACOS(-1.0D0)
  phi(1:ncell)=0.0d0

  if (g%multipole%q(1) <= 0.0d0) return

  mass = g%multipole%q(1)
  center = 0.0d0
  center(1:ndim) = g%multipole%q(2:ndim+1)/mass
#ifdef FMM
#if NDIM==3
  q_central = 0.0d0
  q_central(1,1) = g%multipole%q(5)
  q_central(1,2) = g%multipole%q(6)
  q_central(2,1) = q_central(1,2)
  q_central(1,3) = g%multipole%q(7)
  q_central(3,1) = q_central(1,3)
  q_central(2,2) = g%multipole%q(8)
  q_central(2,3) = g%multipole%q(9)
  q_central(3,2) = q_central(2,3)
  q_central(3,3) = g%multipole%q(10)
  trace_c = q_central(1,1) + q_central(2,2) + q_central(3,3)
#endif
#endif

  do i=1,ncell
     rx=0.0d0; ry=0.0d0; rz=0.0d0
     rx=x(i,1)-center(1)
#if NDIM>1
     ry=x(i,2)-center(2)
#endif
#if NDIM>2
     rz=x(i,3)-center(3)
#endif
     rr=MAX(sqrt(rx**2+ry**2+rz**2),dx)
#if NDIM==1
     phi(i)=mass*fourpi/2d0*rr
#endif
#if NDIM==2
     phi(i)=mass*2d0*log(rr)
#endif
#if NDIM==3
#ifdef FMM
     rr_vec = 0.0d0
     rr_vec(1) = rx
     rr_vec(2) = ry
     rr_vec(3) = rz
     rr2 = rr*rr
     rr5 = rr2*rr2*rr
     scalar_q = 0.0d0
     do idim=1,3
        do j=1,3
           scalar_q = scalar_q + rr_vec(idim)*q_central(idim,j)*rr_vec(j)
        end do
     end do
     phi(i)=-mass/rr - 0.5d0*(3.0d0*scalar_q - trace_c*rr2)/rr5
#else
     phi(i)=-mass/rr
#endif
#endif
  end do

end subroutine phiana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
