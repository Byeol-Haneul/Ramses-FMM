!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(r,g,x,q,dx,nn)
  use amr_parameters, only: ndim, nvector
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t, global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  integer ::nn                            ! Number of cells
  real(kind=8)::dx                         ! Cell size
#ifdef MHD
  real(kind=8),dimension(1:nvector,1:nvar+3-ndim)::q ! Primitive variables
#else
  real(kind=8),dimension(1:nvector,1:nvar)::q ! Primitive variables
#endif
  real(kind=8),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! Kelvin-Helmholtz instability test: density step in z-direction,
  ! velocity shear in x-direction (high-density: +1, low-density: -1)
  ! Sinusoidal z-velocity perturbation at boundaries with parabolic profile
  !================================================================
  integer::i
  real(kind=8)::zz,xx,rho,p0,vx,vz
  real(kind=8)::z_bound1,z_bound2,width_half
  real(kind=8)::profile1,profile2,perturb_amplitude

  p0 = 1.0d0
  z_bound1 = 1.0d0/3.0d0
  z_bound2 = 2.0d0/3.0d0
  width_half = 0.025d0  ! Half-width of 0.05
  perturb_amplitude = 0.1d0  ! Small amplitude for perturbation

  do i=1,nn
     xx = x(i,1)/r%box_size(1)
     zz = x(i,ndim)/r%box_size(ndim)
     rho = 1.0d0
     vx = -1.0d0
     vz = 0.0d0
     if(zz >= z_bound1 .and. zz < z_bound2) then
        rho = 2.0d0
        vx = 1.0d0
     end if
     
     ! Add sinusoidal z-velocity perturbation at boundaries with parabolic profile
     ! Profile centered at z = 1/3
     if(abs(zz - z_bound1) < width_half) then
        profile1 = max(0.0d0, 1.0d0 - ((zz - z_bound1)/width_half)**2)
        vz = vz + perturb_amplitude * profile1 * sin(2.0d0*3.141592653589793d0*xx)
     endif
     ! Profile centered at z = 2/3
     if(abs(zz - z_bound2) < width_half) then
        profile2 = max(0.0d0, 1.0d0 - ((zz - z_bound2)/width_half)**2)
        vz = vz + perturb_amplitude * profile2 * sin(2.0d0*3.141592653589793d0*xx)
     endif
     
     q(i,1) = rho
     q(i,2) = vx
     q(i,3) = 0.0d0
     q(i,4) = 0.0d0
     if(ndim>=2)q(i,3) = 0.0d0  ! y-velocity
     if(ndim>=3)q(i,1+ndim) = vz  ! z-velocity
     q(i,5) = p0
#ifdef MHD
     q(i,nvar+1:nvar+3-ndim) = 0.0d0
#endif
  end do

  ! Compute entropy if needed
  if(r%entropy)then
     q(1:nn,r%ientropy)=q(1:nn,5)/q(1:nn,1)**r%gamma
  endif

  ! Compute metallicity if needed
  if(r%metal)then
     q(1:nn,r%imetal)=r%z_ave*0.02
  endif

end subroutine condinit
