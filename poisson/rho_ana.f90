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

  ! =============================================================
  ! Added variables for two-sphere initial condition
  ! =============================================================
  real(kind=8)::x1mass,y1mass,z1mass
  real(kind=8)::x2mass,y2mass,z2mass
  real(kind=8)::r1,r2
  real(kind=8)::radius1,radius2
  real(kind=8)::rho_val

  ! =============================================================
  ! Branch for custom two-sphere IC when gravity_type(0) = -1
  ! =============================================================

  if (gravity_params(1) < -0.5) then
     ! Sphere 1 center from parameters
     x1mass = gravity_params(2)
     y1mass = gravity_params(3)
     z1mass = gravity_params(4)

     ! Sphere 2 fixed center
     x2mass = gravity_params(5)
     y2mass = gravity_params(6)
     z2mass = gravity_params(7)

     ! Two separate radii
     radius1 = gravity_params(8)
     radius2 = gravity_params(9)

     ! Fixed density
     rho_val = 1.0d0

     ! Distance to sphere 1
     r1 = sqrt( (x(1)-x1mass)**2 &
#if NDIM>1
              + (x(2)-y1mass)**2 &
#endif
#if NDIM>2
              + (x(3)-z1mass)**2 &
#endif
         )

     ! Distance to sphere 2
     r2 = sqrt( (x(1)-x2mass)**2 &
#if NDIM>1
          + (x(2)-y2mass)**2 &
#endif
#if NDIM>2
          + (x(3)-z2mass)**2 &
#endif
        )

     ! Density assignment
     if (r1 <= radius1 .or. r2 <= radius2) then
        d = rho_val
     else
        d = 0.0d0
     endif

     return
  endif

  ! =============================================================
  ! Original analytic density profile (unchanged)
  ! =============================================================

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
