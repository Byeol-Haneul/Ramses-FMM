!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(r,g,x,u,dx,nn)
  use amr_parameters, only: dp, ndim, nvector
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t, global_t
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

#if CONDINIT==coeur
  integer::i
  real(dp)::r2,rx,ry,rz,d,p,vx,vy,vz,r_trunc,r2_trunc,c2
  real(dp)::omega_code,AU,Msol,pi,M,sigma,r_min,r2_min,omega_const,r_vortex,invr2_vortex
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
#endif

  ! Call built-in initial condition generator
  call region_condinit(r,g,x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........

#if CONDINIT==coeur

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

  sigma=M/(4*pi*r_trunc)!*1.00015!*1.016

  r2_trunc=r_trunc**2
  r2_min=r_min**2
  invr2_vortex=1./r_vortex**2

  omega_const=0.1*sqrt(1./r2_trunc+invr2_vortex)*sqrt(M/r_trunc)

  c2=(18939.2/(scale_l/scale_t))**2

  do i=1,nn
     rx=x(i,1)-r%boxlen/2.
     ry=x(i,2)-r%boxlen/2.
     rz=x(i,3)-r%boxlen/2.

     !density/pressure
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
!     vx=-1./(3.**0.5)*omega_code*(ry-rz)
!     vy=-1./(3.**0.5)*omega_code*(rz-rx)
!     vz=-1./(3.**0.5)*omega_code*(rx-ry)

     !velocity
     vx=-1.*omega_code*ry
     vy=omega_code*rx
     vz=0.

     q(i,1)=d
     q(i,2)=vx
     q(i,3)=vy
     q(i,4)=vz
     q(i,5)=p
  end do

#endif

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
