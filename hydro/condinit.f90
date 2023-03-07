!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(r,g,x,u,dx,nn)
  use amr_parameters, only: dp, ndim, nvector
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t, global_t
  use input_hydro_condinit_module, only: region_condinit
  implicit none
  type(run_t)::r
  type(global_t)::g
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
#if NVAR>NDIM+2+NENER
  integer::ivar
#endif
#if NENER>0
  integer::irad
#endif
  real(dp),dimension(1:nvector,1:nvar)::q   ! Primitive variables
#ifdef CONDINITcoeur
  integer::i
  real(dp)::r2,rx,ry,rz,d,p,vx,vy,vz,r_trunc,r2_trunc,c2
  real(dp)::omega_code,AU,Msol,pi,M,sigma,r_min,r2_min,omega_const,r_vortex,invr2_vortex
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
#endif
#ifdef CONDINITinsta
  integer::i,id,iu,iv,iw,ip,ix,iy
  real(dp)::x0,lambday,ky,lambdaz,kz,rho1,rho2,p0,v0,v1,v2
#endif
#ifdef CONDINITdoublemach
  integer::i,id,iu,iv,ip
  real(dp)::pi,xp
#endif
  
#ifdef CONDINITdefault
  ! Call built-in initial condition generator
  call region_condinit(r,g,x,q,dx,nn)
#endif
  
  ! Add here, if you wish, some user-defined initial conditions
  ! ........

#ifdef CONDINITcoeur
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3
  ! constants
  AU=1.49598d13
  Msol= 1.98892d33
  pi=3.14159
  ! mass, radius, and ratio of rotational to gravitational energy
  r_trunc=25*4000.*AU/scale_l
  r_min=10.*AU/scale_l
  r_vortex=4000.*AU/scale_l
  M=100.*Msol/scale_m
  sigma=M/(4*pi*r_trunc)
  r2_trunc=r_trunc**2
  r2_min=r_min**2
  invr2_vortex=1./r_vortex**2
  omega_const=0.1*sqrt(1./r2_trunc+invr2_vortex)*sqrt(M/r_trunc)
  c2=(18939.2/(scale_l/scale_t))**2
  do i=1,nn
     rx=x(i,1)-r%boxlen/2.
     ry=x(i,2)-r%boxlen/2.
     rz=x(i,3)-r%boxlen/2.
     !density
     r2=rx**2+ry**2+rz**2
     d=sigma/(r2+r2_min)
     omega_code=omega_const/sqrt(1.+invr2_vortex*r2)
     if (r2>=r2_trunc)then
        d=d*1.e-4
        omega_code=omega_code/sqrt(r2)*exp(10.*(r2_trunc-r2))
     end if
     !pressure
     p=d*c2
     !velocity
     vx=-omega_code*ry
     vy=omega_code*rx
     vz=0.
     ! primitive variables
     q(i,1)=d
     q(i,2)=vx
     q(i,3)=vy
     q(i,4)=vz
     q(i,5)=p
  end do
#endif

#ifdef CONDINITinsta
  id=1; iu=2; iv=3; iw=4; ip=ndim+2
  x0=r%x_center(1)
  if(r%constant_gravity(2) .ne. 0)then
     ix=2
     iy=1
     iu=3
     iv=2
  else
     ix=1
     iy=2
     iu=2
     iv=3
  endif

  lambday=0.25
  ky=2.*acos(-1.0d0)/lambday
  lambdaz=0.25
  kz=2.*acos(-1.0d0)/lambdaz
  rho1=r%d_region(1)
  rho2=r%d_region(2)
  v1=r%v_region(1)
  v2=r%v_region(2)
  v0=0.1
  p0=10.
  do i=1,nn
     if(x(i,ix) < x0)then
        q(i,id)=rho1
        q(i,iu)=0.0
        if(ndim>1)q(i,iu)=v0*cos(ky*(x(i,iy)-lambday/2.))*exp(+ky*(x(i,ix)-x0))
        if(ndim>1)q(i,iv)=v1
        if(ndim>2)q(i,iw)=0.0D0
        q(i,ip)=p0+rho1*r%constant_gravity(ix)*x(i,ix)
     else
        q(i,id)=rho2
        q(i,iu)=0.0
        if(ndim>1)q(i,iu)=v0*cos(ky*(x(i,iy)-lambday/2.))*exp(-ky*(x(i,ix)-x0))
        if(ndim>1)q(i,iv)=v2
        if(ndim>2)q(i,iw)=0.0D0
        q(i,ip)=p0+rho1*r%constant_gravity(ix)*x0+rho2*r%constant_gravity(ix)*(x(i,ix)-x0)
     endif
  end do
#endif

#ifdef CONDINITdoublemach
  id=1; iu=2; iv=3; ip=ndim+2
  pi=acos(-1.0d0)
  do i=1,nn
     xp=x(i,1)-x(i,2)/tan(pi/3.0)-10./sin(pi/3.0)*g%t
     if(xp<1./6.)then
        q(i,id)=8.
        q(i,iu)=7.145
        q(i,iv)=-4.125
        q(i,ip)=116.5
     else
        q(i,id)=r%gamma
        q(i,iu)=0.0
        q(i,iv)=0.0
        q(i,ip)=1.0
     endif
  end do
#endif

  ! Compute entropy if needed
  if(r%entropy)then
     q(1:nn,r%ientropy)=q(1:nn,ndim+2)/q(1:nn,1)**r%gamma
  endif

  ! Compute metallicity if needed
  if(r%metal)then
     q(1:nn,r%imetal)=r%z_ave*0.02
  endif

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! thermal pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(r%gamma-1.0d0)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do irad=1,nener
     u(1:nn,ndim+2+irad)=q(1:nn,ndim+2+irad)/(r%gamma_rad(irad)-1.0d0)
     u(1:nn,ndim+2)=u(1:nn,ndim+2)+u(1:nn,ndim+2+irad)
  enddo
#endif
#if NVAR>NDIM+2+NENER
  ! passive scalars
  do ivar=ndim+3+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit
