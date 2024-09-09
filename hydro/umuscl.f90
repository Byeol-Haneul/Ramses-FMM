! ---------------------------------------------------------------
!  UNSPLIT     Unsplit second order Godunov integrator for
!              polytropic gas dynamics using either
!              MUSCL-HANCOCK scheme or Collela's PLMDE scheme
!              with various slope limiters.
!
!  inputs/outputs
!  uin         => (const)  input state
!  gravin      => (const)  input gravitational acceleration
!  iu1,iu2     => (const)  first and last index of input array,
!  ju1,ju2     => (const)  cell centered,    
!  ku1,ku2     => (const)  including buffer cells.
!  flux       <=  (modify) return fluxes in the 3 coord directions
!  if1,if2     => (const)  first and last index of output array,
!  jf1,jf2     => (const)  edge centered,
!  kf1,kf2     => (const)  for active cells only.
!  dx,dy,dz    => (const)  (dx,dy,dz)
!  dt          => (const)  time step
!  ndim        => (const)  number of dimensions
! ----------------------------------------------------------------
subroutine unsplit(uin,gravin,qin,cin,flux,tmp,dq,qm,qp,fx,tx,divu,&
#ifdef MHD
     & bin,emfx,emfy,emfz,bf,dbf,Ex,Ey,Ez,qRT,qRB,qLT,qLB,&
     & etamag,induction, &
#endif
     & dx,dy,dz,dt,iu1,iu2,ju1,ju2,ku1,ku2,if1,if2,jf1,jf2,kf1,kf2,&
     & gamma,gamma_rad,smallr,smallc,slope_type,slope_mag_type,riemann,riemann2d,difmag)
  use amr_parameters, only: dp, ndim
  use hydro_parameters, only: nvar, nprim, nener
  use const
  implicit none

  ! Input parameters
  real(dp)::dx,dy,dz,dt
  real(dp)::gamma,smallr,smallc,difmag,etamag
  logical::induction
  real(dp),dimension(1:nener+1)::gamma_rad
  integer::slope_type,slope_mag_type,riemann,riemann2d
  integer::iu1,iu2,ju1,ju2,ku1,ku2
  integer::if1,if2,jf1,jf2,kf1,kf2

  ! Input states
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin

  ! Output fluxes
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2,1:nprim,1:ndim)::flux
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2,1:2    ,1:ndim)::tmp
#ifdef MHD
  ! Input left and right face magnetic field
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:6)::bin

  ! Face-centered magnetic field
  real(dp),dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf

  ! Face-centered magnetic field slopes
  real(dp),dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3,1:2)::dbf

  ! Output electromotive force
  REAL(dp),DIMENSION(if1:if2,jf1:jf2,kf1:kf2)::emfx
  REAL(dp),DIMENSION(if1:if2,jf1:jf2,kf1:kf2)::emfy
  REAL(dp),DIMENSION(if1:if2,jf1:jf2,kf1:kf2)::emfz

  ! Intermediate electromotive force
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2)::Ex
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2)::Ey
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2)::Ez

  ! Edge-averaged left-right and top-bottom state arrays
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:3)::qRT
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:3)::qRB
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:3)::qLT
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:3)::qLB
#endif
  ! Primitive variables
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim)::qin
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2        )::cin

  ! Slopes
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:ndim)::dq

  ! Left and right state arrays
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:ndim)::qm
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:ndim)::qp

  ! Intermediate fluxes
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim)::fx
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:2    )::tx

  ! Velocity divergence
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2)::divu

  ! Local scalar variables
  integer::i,j,k,ivar
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  ! Translate to primative variables, compute sound speeds  
  call ctoprim(uin,qin,cin,gravin, &
#ifdef MHD
       & bin,bf, &
#endif
       & dt,iu1,iu2,ju1,ju2,ku1,ku2,gamma,gamma_rad,smallr,smallc)

  ! Compute TVD slopes
  call uslope(qin,dq,&
#ifdef MHD
       & bf,dbf, &
#endif
       & dx,dt,iu1,iu2,ju1,ju2,ku1,ku2,slope_type,slope_mag_type)

  ! Compute 3D traced-states in all three directions
#if NDIM==1
  call trace1d(qin,dq,qm,qp,dx,dt,iu1,iu2,ju1,ju2,ku1,ku2,gamma,gamma_rad,smallr,smallc)
#endif
#if NDIM==2
  call trace2d(qin,dq,qm,qp,&
#ifdef MHD
       & bf,dbf,Ez,qRT,qRB,qLT,qLB, &
#endif
       & dx,dy,dt,iu1,iu2,ju1,ju2,ku1,ku2,gamma,gamma_rad,smallr,smallc)
#endif
#if NDIM==3
  call trace3d(qin,dq,qm,qp,&
#ifdef MHD
       & bf,dbf,Ex,Ey,Ez,qRT,qRB,qLT,qLB, &
       & induction, &
#endif
       & dx,dy,dz,dt,iu1,iu2,ju1,ju2,ku1,ku2,gamma,gamma_rad,smallr,smallc)
#endif

  ! Solve for 1D flux in X direction
  call cmpflxm(qm,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          if1  ,if2  ,jlo  ,jhi  ,klo  ,khi  , &
#ifdef MHD
       &       6,7,8, &
#endif
       &       2,3,4,fx,tx,gamma,gamma_rad,smallr,smallc,riemann)
  ! Save flux in output array
  do ivar=1,nprim
     do k=klo,khi
        do j=jlo,jhi
           do i=if1,if2
              flux(i,j,k,ivar,1)=fx(i,j,k,ivar)*dt/dx
           end do
        end do
     end do
  end do
  do k=klo,khi
     do j=jlo,jhi
        do i=if1,if2
           do ivar=1,2
              tmp(i,j,k,ivar,1)=tx(i,j,k,ivar)*dt/dx
           end do
        end do
     end do
  end do

  ! Solve for 1D flux in Y direction
#if NDIM>1
  call cmpflxm(qm,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jf1  ,jf2  ,klo  ,khi  , &
#ifdef MHD
       &       7,6,8, &
#endif
       &       3,2,4,fx,tx,gamma,gamma_rad,smallr,smallc,riemann)
  ! Save flux in output array
  do ivar=1,nprim
     do k=klo,khi
        do j=jf1,jf2
           do i=ilo,ihi
              flux(i,j,k,ivar,2)=fx(i,j,k,ivar)*dt/dy
           end do
        end do
     end do
  end do
  do ivar=1,2
     do k=klo,khi
        do j=jf1,jf2
           do i=ilo,ihi
              tmp(i,j,k,ivar,2)=tx(i,j,k,ivar)*dt/dy
           end do
        end do
     end do
  end do
#endif

  ! Solve for 1D flux in Z direction
#if NDIM>2
  call cmpflxm(qm,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jlo  ,jhi  ,kf1  ,kf2  , &
#ifdef MHD
       &       8,6,7, &
#endif
       &       4,2,3,fx,tx,gamma,gamma_rad,smallr,smallc,riemann)
  ! Save flux in output array
  do ivar=1,nprim
     do k=kf1,kf2
        do j=jlo,jhi
           do i=ilo,ihi
              flux(i,j,k,ivar,3)=fx(i,j,k,ivar)*dt/dz
           end do
        end do
     end do
  end do
  do ivar=1,2
     do k=kf1,kf2
        do j=jlo,jhi
           do i=ilo,ihi
              tmp(i,j,k,ivar,3)=tx(i,j,k,ivar)*dt/dz
           end do
        end do
     end do
  end do
#endif
#ifdef MHD
#if NDIM>1
  ! Solve for EMF in Z direction
  CALL cmp_mag_flx(qRT,iu1+1,iu2+1,ju1+1,ju2+1,ku1  ,ku2  , &
       &           qRB,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &           qLT,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &           qLB,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &               if1  ,if2  ,jf1  ,jf2  ,klo  ,khi  , &
       &           2,3,4,6,7,8,Ez,gamma,gamma_rad,smallr,smallc,riemann2d)

  ! Save vector in output array
  do k=klo,khi
     do j=jf1,jf2
        do i=if1,if2
           emfz(i,j,k)=Ez(i,j,k)*dt/dx
        end do
     end do
  end do
#endif
#if NDIM>2
  ! Solve for EMF in Y direction
  CALL cmp_mag_flx(qRT,iu1+1,iu2+1,ju1,ju2,ku1+1,ku2+1, &
       &           qLT,iu1  ,iu2  ,ju1,ju2,ku1+1,ku2+1, &
       &           qRB,iu1+1,iu2+1,ju1,ju2,ku1  ,ku2  , &
       &           qLB,iu1  ,iu2  ,ju1,ju2,ku1  ,ku2  , &
       &               if1  ,if2  ,jlo,jhi,kf1  ,kf2  , &
       &           4,2,3,8,6,7,Ey,gamma,gamma_rad,smallr,smallc,riemann2d)
  ! Save vector in output array
  do k=kf1,kf2
     do j=jlo,jhi
        do i=if1,if2
           emfy(i,j,k)=Ey(i,j,k)*dt/dx
        end do
     end do
  end do
  ! Solve for EMF in X direction
  CALL cmp_mag_flx(qRT,iu1,iu2,ju1+1,ju2+1,ku1+1,ku2+1, &
       &           qRB,iu1,iu2,ju1+1,ju2+1,ku1  ,ku2  , &
       &           qLT,iu1,iu2,ju1  ,ju2  ,ku1+1,ku2+1, &
       &           qLB,iu1,iu2,ju1  ,ju2  ,ku1  ,ku2  , &
       &               ilo,ihi,jf1  ,jf2  ,kf1  ,kf2  , &
       &           3,4,2,7,8,6,Ex,gamma,gamma_rad,smallr,smallc,riemann2d)
  ! Save vector in output array
  do k=kf1,kf2
     do j=jf1,jf2
        do i=ilo,ihi
           emfx(i,j,k)=Ex(i,j,k)*dt/dx
        end do
     end do
  end do
#endif
#endif
  if(difmag>0.0)then
     call cmpdivu(qin,divu,dx,dy,dz,&
          & iu1,iu2,ju1,ju2,ku1,ku2,&
          & if1,if2,jf1,jf2,kf1,kf2)
     call cmpdiff(uin,flux,divu,dt,&
#ifdef MHD
          & bin, &
#endif
          & iu1,iu2,ju1,ju2,ku1,ku2,&
          & if1,if2,jf1,jf2,kf1,kf2,difmag)
  endif
#ifdef MHD
  if(etamag>0.0)then
     call cmpcurrent(bf,emfx,emfy,emfz,dx,dy,dz,dt,&
          & iu1,iu2,ju1,ju2,ku1,ku2,&
          & if1,if2,jf1,jf2,kf1,kf2,etamag)
  endif
#endif

end subroutine unsplit
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine trace1d(q,dq,qm,qp,dx,dt,iu1,iu2,ju1,ju2,ku1,ku2,gamma,gamma_rad,smallr,smallc)
  use amr_parameters, only: dp, ndim
  use hydro_parameters, only: nprim, nener, ie
  use const
  implicit none

  real(dp)::dx, dt
  integer::iu1, iu2, ju1, ju2, ku1, ku2
  real(dp)::gamma, smallr, smallc
  real(dp),dimension(1:nener+1)::gamma_rad
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim)::q
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:ndim)::dq
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:ndim)::qm
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:ndim)::qp

  ! Local variables
  integer::i, j, k, n
  integer::ilo, ihi, jlo, jhi, klo, khi
  integer::ir, iu, iv, iw, ip
  real(dp)::dtdx
  real(dp)::r, u, v, w, p
  real(dp)::drx, dux, dvx, dwx, dpx
  real(dp)::sr0, su0, sv0, sw0, sp0
#ifdef MHD
  integer::iA, iB, iC
  real(dp)::A, B, C
  real(dp)::dAx, dBx, dCx
  real(dp)::sA0, sB0, sC0
#endif
#if NENER>0
  integer::irad
  real(dp),dimension(1:nener)::e, dex, se0
#endif

  dtdx = dt/dx

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  ir=1; iu=2; iv=3; iw=4; ip=5
#ifdef MHD
  iA=6; iB=7; iC=8
#endif

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi

           ! Cell centered values
           r = q(i,j,k,ir)
           u = q(i,j,k,iu)
           v = q(i,j,k,iv)
           w = q(i,j,k,iw)
           p = q(i,j,k,ip)
#ifdef MHD
           A = q(i,j,k,iA)
           B = q(i,j,k,iB)
           C = q(i,j,k,iC)
#endif
#if NENER>0
           do irad = 1,nener
              e(irad) = q(i,j,k,ie+irad)
           end do
#endif
           ! TVD slopes in X direction
           drx = half*dq(i,j,k,ir,1)
           dux = half*dq(i,j,k,iu,1)
           dvx = half*dq(i,j,k,iv,1)
           dwx = half*dq(i,j,k,iw,1)
           dpx = half*dq(i,j,k,ip,1)
#ifdef MHD
           dBx = half*dq(i,j,k,iB,1)
           dCx = half*dq(i,j,k,iC,1)
#endif
#if NENER>0
           do irad = 1,nener
              dex(irad) = half*dq(i,j,k,ie+irad,1)
           end do
#endif
           ! Source terms (including transverse derivatives)
           sr0 = -u*drx - (dux)*r
           sp0 = -u*dpx - (dux)*gamma*p
           su0 = -u*dux - (dpx)/r
           sv0 = -u*dvx
           sw0 = -u*dwx
#ifdef MHD
           su0 = su0 - (B*dBx+C*dCx)/r
           sv0 = sv0 + (A*dBx)/r
           sw0 = sw0 + (A*dCx)/r
           sB0 = -u*dBx + A*dvx - B*dux
           sC0 = -u*dCx + A*dwx - C*dux
#endif
#if NENER>0
           do irad = 1,nener
              su0 = su0 - (dex(irad))/r
              se0(irad) = -u*dex(irad) - (dux)*gamma_rad(irad)*e(irad)
           end do
#endif
           ! Cell-centered predicted states
           r = r + sr0*dtdx
           u = u + su0*dtdx
           v = v + sv0*dtdx
           w = w + sw0*dtdx
           p = p + sp0*dtdx
#ifdef MHD
           B = B + sB0*dtdx
           C = C + sC0*dtdx
#endif
#if NENER>0
           do irad = 1,nener
              e(irad) = e(irad) + se0(irad)*dtdx
           end do
#endif
           ! Right state
           qp(i,j,k,ir,1) = r - drx
           qp(i,j,k,iu,1) = u - dux
           qp(i,j,k,iv,1) = v - dvx
           qp(i,j,k,iw,1) = w - dwx
           qp(i,j,k,ip,1) = p - dpx
           if(qp(i,j,k,ir,1)<smallr)qp(i,j,k,ir,1)=q(i,j,k,ir)
#ifdef MHD
           qp(i,j,k,iA,1) = A
           qp(i,j,k,iB,1) = B - dBx
           qp(i,j,k,iC,1) = C - dCx
#endif
#if NENER>0
           do irad=1,nener
              qp(i,j,k,ie+irad,1) = e(irad) - dex(irad)
           end do
#endif
           ! Left state
           qm(i,j,k,ir,1) = r + drx
           qm(i,j,k,iu,1) = u + dux
           qm(i,j,k,iv,1) = v + dvx
           qm(i,j,k,iw,1) = w + dwx
           qm(i,j,k,ip,1) = p + dpx
           if(qm(i,j,k,ir,1)<smallr)qm(i,j,k,ir,1)=q(i,j,k,ir)
#ifdef MHD
           qm(i,j,k,iA,1) = A
           qm(i,j,k,iB,1) = B + dBx
           qm(i,j,k,iC,1) = C + dCx
#endif
#if NENER>0
           do irad=1,nener
              qm(i,j,k,ie+irad,1) = e(irad) + dex(irad)
           end do
#endif
        end do
     end do
  end do

#if NVAR>5+NENER
  ! Passive scalars
  do n = ie+nener+1, nprim
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              r   = q(i,j,k,n)          ! Cell centered values
              u   = q(i,j,k,iu)
              drx = half*dq(i,j,k,n,1)  ! TVD slopes
              sr0 = -u*drx              ! Source terms
              r   = r + sr0*dtdx
              qp(i,j,k,n,1) = r - drx   ! Right state
              qm(i,j,k,n,1) = r + drx   ! Left state
           end do
        end do
     end do
  end do
#endif

end subroutine trace1d
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>1
subroutine trace2d(q,dq,qm,qp, &
#ifdef MHD
       & bf,dbf,Ez,qRT,qRB,qLT,qLB, &
#endif
       & dx,dy,dt,iu1,iu2,ju1,ju2,ku1,ku2,gamma,gamma_rad,smallr,smallc)
  use amr_parameters, only:dp, ndim
  use hydro_parameters, only: nprim, nener, ie
  use const
  implicit none

  real(dp)::dx, dy, dt
  integer::iu1,iu2,ju1,ju2,ku1,ku2
  real(dp)::gamma,smallr,smallc
  real(dp),dimension(1:nener+1)::gamma_rad
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim)::q
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:ndim)::dq
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:ndim)::qm
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:ndim)::qp
#ifdef MHD
  real(dp),dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf
  real(dp),dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3,1:2)::dbf
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2)::Ez
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:3)::qRT
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:3)::qRB
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:3)::qLT
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:3)::qLB
#endif

  ! declare local variables
  integer::i, j, k, n
  integer::ilo,ihi,jlo,jhi,klo,khi
  integer::ir, iu, iv, iw, ip
  real(dp)::dtdx
  real(dp)::r, u, v, w, p
  real(dp)::drx, dux, dvx, dwx, dpx
  real(dp)::dry, duy, dvy, dwy, dpy
  real(dp)::sr0, su0, sv0, sw0, sp0
#ifdef MHD
  integer::iA, iB, iC
  real(dp)::A, B, C
  real(dp)::dAx, dBx, dCx
  real(dp)::dAy, dBy, dCy
  real(dp)::sA0, sB0, sC0
  real(dp)::AL, AR, BL, BR
  real(dp)::dALy, dARy, dBLx, dBRx
  real(dp)::sAL0, sAR0, sBL0, sBR0
  real(dp)::ELL, ELR, ERL, ERR
#endif
#if NENER>0
  integer::irad
  real(dp),dimension(1:nener)::e, dex, dey, se0
#endif
  
  dtdx = dt/dx
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4; ip=5

#ifdef MHD
  iA=6; iB=7; iC=8
  DO k = klo, ku2
     DO j = jlo, ju2
        DO i = ilo, iu2
           u = 0.25*(q(i-1,j-1,k,iu)+q(i-1,j,k,iu)+q(i,j-1,k,iu)+q(i,j,k,iu))
           v = 0.25*(q(i-1,j-1,k,iv)+q(i-1,j,k,iv)+q(i,j-1,k,iv)+q(i,j,k,iv))
           A = 0.5*(bf(i,j-1,k,1)+bf(i,j,k,1))
           B = 0.5*(bf(i-1,j,k,2)+bf(i,j,k,2))
           Ez(i,j,k)=u*B-v*A
        END DO
     END DO
  END DO
#endif
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
#ifdef MHD
           ! Face centered variables
           AL =  bf(i  ,j  ,k,1)
           AR =  bf(i+1,j  ,k,1)
           BL =  bf(i  ,j  ,k,2)
           BR =  bf(i  ,j+1,k,2)

           ! Face centered TVD slopes in transverse direction
           dALy = half*dbf(i  ,j  ,k,1,1)
           dARy = half*dbf(i+1,j  ,k,1,1)
           dBLx = half*dbf(i  ,j  ,k,2,1)
           dBRx = half*dbf(i  ,j+1,k,2,1)

           ! Edge centered electric field Ez = uB-vA
           ELL = Ez(i  ,j  ,k)
           ELR = Ez(i  ,j+1,k)
           ERL = Ez(i+1,j  ,k)
           ERR = Ez(i+1,j+1,k)

           ! Face-centered predicted states
           sAL0 = +(ELR-ELL)
           sAR0 = +(ERR-ERL)
           sBL0 = -(ERL-ELL)
           sBR0 = -(ERR-ELR)
           
           AL = AL + sAL0*dtdx*half
           AR = AR + sAR0*dtdx*half
           BL = BL + sBL0*dtdx*half
           BR = BR + sBR0*dtdx*half
#endif
           ! Cell centered values
           r = q(i,j,k,ir)
           u = q(i,j,k,iu)
           v = q(i,j,k,iv)
           w = q(i,j,k,iw)
           p = q(i,j,k,ip)
#ifdef MHD
           A = q(i,j,k,iA)
           B = q(i,j,k,iB)
           C = q(i,j,k,iC)
#endif
#if NENER>0
           do irad=1,nener
              e(irad) = q(i,j,k,ie+irad)
           end do
#endif
           ! TVD slopes in all directions
           drx = half*dq(i,j,k,ir,1)
           dux = half*dq(i,j,k,iu,1)
           dvx = half*dq(i,j,k,iv,1)
           dwx = half*dq(i,j,k,iw,1)
           dpx = half*dq(i,j,k,ip,1)
#ifdef MHD
           dBx = half*dq(i,j,k,iB,1)
           dCx = half*dq(i,j,k,iC,1)
#endif
#if NENER>0
           do irad=1,nener
              dex(irad) = half*dq(i,j,k,ie+irad,1)
           end do
#endif
           dry = half*dq(i,j,k,ir,2)
           duy = half*dq(i,j,k,iu,2)
           dvy = half*dq(i,j,k,iv,2)
           dwy = half*dq(i,j,k,iw,2)
           dpy = half*dq(i,j,k,ip,2)
#ifdef MHD
           dAy = half*dq(i,j,k,iA,2)
           dCy = half*dq(i,j,k,iC,2)
#endif
#if NENER>0
           do irad=1,nener
              dey(irad) = half*dq(i,j,k,ie+irad,2)
           end do
#endif
           ! source terms (with transverse derivatives)
           sr0 = -u*drx-v*dry - (dux+dvy)*r
           sp0 = -u*dpx-v*dpy - (dux+dvy)*gamma*p
           su0 = -u*dux-v*duy - (dpx    )/r
           sv0 = -u*dvx-v*dvy - (dpy    )/r
           sw0 = -u*dwx-v*dwy
#ifdef MHD
           su0 = su0 + (-(B*dBx+C*dCx)/r) + (B*dAy/r)
           sv0 = sv0 + (A*dBx/r) + (-(A*dAy+C*dCy)/r)
           sw0 = sw0 + (A*dCx/r) + (B*dCy/r)
           sC0 = (-u*dCx-C*dux+A*dwx) + (-v*dCy-C*dvy+B*dwy)
#endif
#if NENER>0
           do irad=1,nener
              su0 = su0 - (dex(irad))/r
              sv0 = sv0 - (dey(irad))/r
              se0(irad) = -u*dex(irad)-v*dey(irad) &
                   & - (dux+dvy)*gamma_rad(irad)*e(irad)
           end do
#endif
           ! Cell-centered predicted states
           r = r + sr0*dtdx
           u = u + su0*dtdx
           v = v + sv0*dtdx
           w = w + sw0*dtdx
           p = p + sp0*dtdx
#ifdef MHD
           A = half*(AL+AR)
           B = half*(BL+BR)
           C = C + sC0*dtdx
#endif
#if NENER>0
           do irad=1,nener
              e(irad) = e(irad) + se0(irad)*dtdx
           end do
#endif
           ! Right state at left interface
           qp(i,j,k,ir,1) = r - drx
           qp(i,j,k,iu,1) = u - dux
           qp(i,j,k,iv,1) = v - dvx
           qp(i,j,k,iw,1) = w - dwx
           qp(i,j,k,ip,1) = p - dpx
#ifdef MHD
           qp(i,j,k,iA,1) = AL
           qp(i,j,k,iB,1) = B - dBx
           qp(i,j,k,iC,1) = C - dCx
#endif
           if(qp(i,j,k,ir,1)<smallr)qp(i,j,k,ir,1)=q(i,j,k,ir)
#if NENER>0
           do irad=1,nener
              qp(i,j,k,ie+irad,1) = e(irad) - dex(irad)
           end do
#endif
           ! Left state at right interface
           qm(i,j,k,ir,1) = r + drx
           qm(i,j,k,iu,1) = u + dux
           qm(i,j,k,iv,1) = v + dvx
           qm(i,j,k,iw,1) = w + dwx
           qm(i,j,k,ip,1) = p + dpx
#ifdef MHD
           qm(i,j,k,iA,1) = AR
           qm(i,j,k,iB,1) = B + dBx
           qm(i,j,k,iC,1) = C + dCx
#endif
           if(qm(i,j,k,ir,1)<smallr)qm(i,j,k,ir,1)=q(i,j,k,ir)
#if NENER>0
           do irad=1,nener
              qm(i,j,k,ie+irad,1) = e(irad) + dex(irad)
           end do
#endif
           ! Top state at bottom interface
           qp(i,j,k,ir,2) = r - dry
           qp(i,j,k,iu,2) = u - duy
           qp(i,j,k,iv,2) = v - dvy
           qp(i,j,k,iw,2) = w - dwy
           qp(i,j,k,ip,2) = p - dpy
#ifdef MHD
           qp(i,j,k,iA,2) = A - dAy
           qp(i,j,k,iB,2) = BL
           qp(i,j,k,iC,2) = C - dCy
#endif
           if(qp(i,j,k,ir,2)<smallr)qp(i,j,k,ir,2)=q(i,j,k,ir)
#if NENER>0
           do irad=1,nener
              qp(i,j,k,ie+irad,2) = e(irad) - dey(irad)
           end do
#endif
           ! Bottom state at top interface
           qm(i,j,k,ir,2) = r + dry
           qm(i,j,k,iu,2) = u + duy
           qm(i,j,k,iv,2) = v + dvy
           qm(i,j,k,iw,2) = w + dwy
           qm(i,j,k,ip,2) = p + dpy
#ifdef MHD
           qm(i,j,k,iA,2) = A + dAy
           qm(i,j,k,iB,2) = BR
           qm(i,j,k,iC,2) = C + dCy
#endif
           if(qm(i,j,k,ir,2)<smallr)qm(i,j,k,ir,2)=q(i,j,k,ir)
#if NENER>0
           do irad=1,nener
              qm(i,j,k,ie+irad,2) = e(irad) + dey(irad)
           end do
#endif
#ifdef MHD
           ! Edge averaged right-top corner state (RT->LL)
           qRT(i,j,k,ir,3) = r + (+drx+dry)
           qRT(i,j,k,iu,3) = u + (+dux+duy)
           qRT(i,j,k,iv,3) = v + (+dvx+dvy)
           qRT(i,j,k,iw,3) = w + (+dwx+dwy)
           qRT(i,j,k,ip,3) = p + (+dpx+dpy)
           qRT(i,j,k,iC,3) = C + (+dCx+dCy)
           qRT(i,j,k,iA,3) = AR+ (   +dARy)
           qRT(i,j,k,iB,3) = BR+ (+dBRx   )
           if (qRT(i,j,k,ir,3)<smallr) qRT(i,j,k,ir,3)=r
#if NENER>0
           do irad=1,nener
              qRT(i,j,k,iC+irad,3) = e(irad) + (+dex(irad)+dey(irad))
           end do
#endif
           ! Edge averaged right-bottom corner state (RB->LR)
           qRB(i,j,k,ir,3) = r + (+drx-dry)
           qRB(i,j,k,iu,3) = u + (+dux-duy)
           qRB(i,j,k,iv,3) = v + (+dvx-dvy)
           qRB(i,j,k,iw,3) = w + (+dwx-dwy)
           qRB(i,j,k,ip,3) = p + (+dpx-dpy)
           qRB(i,j,k,iC,3) = C + (+dCx-dCy)
           qRB(i,j,k,iA,3) = AR+ (   -dARy)
           qRB(i,j,k,iB,3) = BL+ (+dBLx   )
           if (qRB(i,j,k,ir,3)<smallr) qRB(i,j,k,ir,3)=r
#if NENER>0
           do irad=1,nener
              qRB(i,j,k,iC+irad,3) = e(irad) + (+dex(irad)-dey(irad))
           end do
#endif
           ! Edge averaged left-top corner state (LT->RL)
           qLT(i,j,k,ir,3) = r + (-drx+dry)
           qLT(i,j,k,iu,3) = u + (-dux+duy)
           qLT(i,j,k,iv,3) = v + (-dvx+dvy)
           qLT(i,j,k,iw,3) = w + (-dwx+dwy)
           qLT(i,j,k,ip,3) = p + (-dpx+dpy)
           qLT(i,j,k,iC,3) = C + (-dCx+dCy)
           qLT(i,j,k,iA,3) = AL+ (   +dALy)
           qLT(i,j,k,iB,3) = BR+ (-dBRx   )
           if (qLT(i,j,k,ir,3)<smallr) qLT(i,j,k,ir,3)=r
#if NENER>0
           do irad=1,nener
              qLT(i,j,k,iC+irad,3) = e(irad) + (-dex(irad)+dey(irad))
           end do
#endif
           ! Edge averaged left-bottom corner state (LB->RR)
           qLB(i,j,k,ir,3) = r + (-drx-dry)
           qLB(i,j,k,iu,3) = u + (-dux-duy)
           qLB(i,j,k,iv,3) = v + (-dvx-dvy)
           qLB(i,j,k,iw,3) = w + (-dwx-dwy)
           qLB(i,j,k,ip,3) = p + (-dpx-dpy)
           qLB(i,j,k,iC,3) = C + (-dCx-dCy)
           qLB(i,j,k,iA,3) = AL+ (   -dALy)
           qLB(i,j,k,iB,3) = BL+ (-dBLx   )
           if (qLB(i,j,k,ir,3)<smallr) qLB(i,j,k,ir,3)=r
#if NENER>0
           do irad=1,nener
              qLB(i,j,k,iC+irad,3) = e(irad) + (-dex(irad)-dey(irad))
           end do
#endif
#endif
        end do
     end do
  end do

#if NVAR>5+NENER
  ! Passive scalars
  do n = ie+nener+1, nprim
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              r   = q(i,j,k,n)          ! Cell centered values
              u   = q(i,j,k,iu)
              v   = q(i,j,k,iv)
              drx = half*dq(i,j,k,n,1)  ! TVD slopes
              dry = half*dq(i,j,k,n,2)
              sr0 = -u*drx-v*dry        ! Source terms
              r   = r + sr0*dtdx
              qp(i,j,k,n,1) = r - drx   ! Right state
              qm(i,j,k,n,1) = r + drx   ! Left state
              qp(i,j,k,n,2) = r - dry   ! Top state
              qm(i,j,k,n,2) = r + dry   ! Bottom state
           end do
        end do
     end do
  end do
#endif

end subroutine trace2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>2
subroutine trace3d(q,dq,qm,qp, &
#ifdef MHD
     & bf,dbf,Ex,Ey,Ez,qRT,qRB,qLT,qLB, &
     & induction, &
#endif
     & dx,dy,dz,dt,iu1,iu2,ju1,ju2,ku1,ku2,gamma,gamma_rad,smallr,smallc)
  use amr_parameters, only: dp, ndim
  use hydro_parameters, only: nprim, nener, ie
  use const
  implicit none

  real(dp)::dx, dy, dz, dt
  integer::iu1,iu2,ju1,ju2,ku1,ku2
  real(dp)::gamma,smallr,smallc
  real(dp),dimension(1:nener+1)::gamma_rad
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim)::q
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:ndim)::dq
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:ndim)::qm
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:ndim)::qp
#ifdef MHD
  logical::induction
  real(dp),dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf
  real(dp),dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3,1:2)::dbf
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2)::Ex
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2)::Ey
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2)::Ez
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:3)::qRT
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:3)::qRB
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:3)::qLT
  REAL(dp),DIMENSION(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:3)::qLB
#endif

  ! declare local variables
  integer::i, j, k, n
  integer::ilo,ihi,jlo,jhi,klo,khi
  integer::ir, iu, iv, iw, ip
  real(dp)::dtdx
  real(dp)::r, u, v, w, p
  real(dp)::drx, dux, dvx, dwx, dpx
  real(dp)::dry, duy, dvy, dwy, dpy
  real(dp)::drz, duz, dvz, dwz, dpz
  real(dp)::sr0, su0, sv0, sw0, sp0
#ifdef MHD
  integer::iA, iB, iC
  real(dp)::A, B, C
  real(dp)::dAx, dBx, dCx
  real(dp)::dAy, dBy, dCy
  real(dp)::dAz, dBz, dCz
  real(dp)::sA0, sB0, sC0
  real(dp)::AL, AR, BL, BR, CL, CR
  REAL(dp)::dALy, dARy, dALz, dARz
  REAL(dp)::dBLx, dBRx, dBLz, dBRz
  REAL(dp)::dCLx, dCRx, dCLy, dCRy
  real(dp)::sAL0, sAR0, sBL0, sBR0, sCL0, sCR0
  REAL(dp)::ELL, ELR, ERL, ERR
  REAL(dp)::FLL, FLR, FRL, FRR
  REAL(dp)::GLL, GLR, GRL, GRR
#endif
#if NENER>0
  integer::irad
  real(dp),dimension(1:nener)::e, dex, dey, dez, se0
#endif
  
  dtdx = dt/dx
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4; ip=5

#ifdef MHD
  iA=6; iB=7; iC=8
  DO k = klo, ku2
     DO j = jlo, ju2
        DO i = ilo, iu2
           v = 0.25*(q(i,j-1,k-1,iv)+q(i,j-1,k,iv)+q(i,j,k-1,iv)+q(i,j,k,iv))
           w = 0.25*(q(i,j-1,k-1,iw)+q(i,j-1,k,iw)+q(i,j,k-1,iw)+q(i,j,k,iw))
           B = 0.5*(bf(i,j,k-1,2)+bf(i,j,k,2))
           C = 0.5*(bf(i,j-1,k,3)+bf(i,j,k,3))
           Ex(i,j,k) = v*C-w*B

           u = 0.25*(q(i-1,j,k-1,iu)+q(i-1,j,k,iu)+q(i,j,k-1,iu)+q(i,j,k,iu))
           w = 0.25*(q(i-1,j,k-1,iw)+q(i-1,j,k,iw)+q(i,j,k-1,iw)+q(i,j,k,iw))
           A = 0.5*(bf(i,j,k-1,1)+bf(i,j,k,1))
           C = 0.5*(bf(i-1,j,k,3)+bf(i,j,k,3))
           Ey(i,j,k) = w*A-u*C

           u = 0.25*(q(i-1,j-1,k,iu)+q(i-1,j,k,iu)+q(i,j-1,k,iu)+q(i,j,k,iu))
           v = 0.25*(q(i-1,j-1,k,iv)+q(i-1,j,k,iv)+q(i,j-1,k,iv)+q(i,j,k,iv))
           A = 0.5*(bf(i,j-1,k,1)+bf(i,j,k,1))
           B = 0.5*(bf(i-1,j,k,2)+bf(i,j,k,2))
           Ez(i,j,k) = u*B-v*A
        END DO
     END DO
  END DO
#endif
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
#ifdef MHD
           ! Face centered variables
           AL = bf(i  ,j  ,k  ,1)
           AR = bf(i+1,j  ,k  ,1)
           BL = bf(i  ,j  ,k  ,2)
           BR = bf(i  ,j+1,k  ,2)
           CL = bf(i  ,j  ,k  ,3)
           CR = bf(i  ,j  ,k+1,3)

           ! Face centered TVD slopes in transverse direction
           dALy = half*dbf(i  ,j  ,k  ,1,1)
           dARy = half*dbf(i+1,j  ,k  ,1,1)
           dALz = half*dbf(i  ,j  ,k  ,1,2)
           dARz = half*dbf(i+1,j  ,k  ,1,2)

           dBLx = half*dbf(i  ,j  ,k  ,2,1)
           dBRx = half*dbf(i  ,j+1,k  ,2,1)
           dBLz = half*dbf(i  ,j  ,k  ,2,2)
           dBRz = half*dbf(i  ,j+1,k  ,2,2)

           dCLx = half*dbf(i  ,j  ,k  ,3,1)
           dCRx = half*dbf(i  ,j  ,k+1,3,1)
           dCLy = half*dbf(i  ,j  ,k  ,3,2)
           dCRy = half*dbf(i  ,j  ,k+1,3,2)

           ! Edge centered electric field Ex = vC-wB
           ELL = Ex(i,j  ,k  )
           ELR = Ex(i,j  ,k+1)
           ERL = Ex(i,j+1,k  )
           ERR = Ex(i,j+1,k+1)

           ! Edge centered electric field Ey = wA-uC
           FLL = Ey(i  ,j,k  )
           FLR = Ey(i  ,j,k+1)
           FRL = Ey(i+1,j,k  )
           FRR = Ey(i+1,j,k+1)

           ! Edge centered electric field Ez = uB-vA
           GLL = Ez(i  ,j  ,k)
           GLR = Ez(i  ,j+1,k)
           GRL = Ez(i+1,j  ,k)
           GRR = Ez(i+1,j+1,k)

           ! Face-centered predicted states
           sAL0 = + (GLR-GLL) - (FLR-FLL)
           sAR0 = + (GRR-GRL) - (FRR-FRL)
           sBL0 = - (GRL-GLL) + (ELR-ELL)
           sBR0 = - (GRR-GLR) + (ERR-ERL)
           sCL0 = + (FRL-FLL) - (ERL-ELL)
           sCR0 = + (FRR-FLR) - (ERR-ELR)

           AL = AL + sAL0*dtdx*half
           AR = AR + sAR0*dtdx*half
           BL = BL + sBL0*dtdx*half
           BR = BR + sBR0*dtdx*half
           CL = CL + sCL0*dtdx*half
           CR = CR + sCR0*dtdx*half
#endif
           ! Cell centered values
           r = q(i,j,k,ir)
           u = q(i,j,k,iu)
           v = q(i,j,k,iv)
           w = q(i,j,k,iw)
           p = q(i,j,k,ip)
#ifdef MHD
           A = q(i,j,k,iA)
           B = q(i,j,k,iB)
           C = q(i,j,k,iC)
#endif
#if NENER>0
           do irad=1,nener
              e(irad) = q(i,j,k,ie+irad)
           end do
#endif
           ! TVD slopes in all 3 directions
           drx = half*dq(i,j,k,ir,1)
           dux = half*dq(i,j,k,iu,1)
           dvx = half*dq(i,j,k,iv,1)
           dwx = half*dq(i,j,k,iw,1)
           dpx = half*dq(i,j,k,ip,1)
#ifdef MHD
           dBx = half*dq(i,j,k,iB,1)
           dCx = half*dq(i,j,k,iC,1)
#endif
#if NENER>0
           do irad=1,nener
              dex(irad) = half*dq(i,j,k,ie+irad,1)
           end do
#endif
           dry = half*dq(i,j,k,ir,2)
           duy = half*dq(i,j,k,iu,2)
           dvy = half*dq(i,j,k,iv,2)
           dwy = half*dq(i,j,k,iw,2)
           dpy = half*dq(i,j,k,ip,2)
#ifdef MHD
           dAy = half*dq(i,j,k,iA,2)
           dCy = half*dq(i,j,k,iC,2)
#endif
#if NENER>0
           do irad=1,nener
              dey(irad) = half*dq(i,j,k,ie+irad,2)
           end do
#endif
           drz = half*dq(i,j,k,ir,3)
           duz = half*dq(i,j,k,iu,3)
           dvz = half*dq(i,j,k,iv,3)
           dwz = half*dq(i,j,k,iw,3)
           dpz = half*dq(i,j,k,ip,3)
#ifdef MHD
           dAz = half*dq(i,j,k,iA,3)
           dBz = half*dq(i,j,k,iB,3)
#endif
#if NENER>0
           do irad=1,nener
              dez(irad) = half*dq(i,j,k,ie+irad,3)
           end do
#endif
           ! Source terms (including transverse derivatives)
           sr0 = -u*drx-v*dry-w*drz - (dux+dvy+dwz)*r
           su0 = -u*dux-v*duy-w*duz - (dpx        )/r
           sv0 = -u*dvx-v*dvy-w*dvz - (dpy        )/r
           sw0 = -u*dwx-v*dwy-w*dwz - (dpz        )/r
           sp0 = -u*dpx-v*dpy-w*dpz - (dux+dvy+dwz)*gamma*p
#ifdef MHD
           su0 = su0 + (-(B*dBx+C*dCx)/r) + (B*dAy/r) + (C*dAz/r)
           sv0 = sv0 + (A*dBx/r) + (-(A*dAy+C*dCy)/r) + (C*dBz/r)
           sw0 = sw0 + (A*dCx/r) + (B*dCy/r) + (-(A*dAz+B*dBz)/r)
           if(induction)then
              su0=0
              sv0=0
              sw0=0
           endif
#endif
#if NENER>0
           do irad=1,nener
              su0 = su0 - (dex(irad))/r
              sv0 = sv0 - (dey(irad))/r
              sw0 = sw0 - (dez(irad))/r
              se0(irad) = -u*dex(irad)-v*dey(irad)-w*dez(irad) & 
                   & - (dux+dvy+dwz)*gamma_rad(irad)*e(irad)
           end do
#endif
           ! Cell-centered predicted states
           r = r + sr0*dtdx
           u = u + su0*dtdx
           v = v + sv0*dtdx
           w = w + sw0*dtdx
           p = p + sp0*dtdx
#ifdef MHD
           A = 0.5*(AL+AR)
           B = 0.5*(BL+BR)
           C = 0.5*(CL+CR)
#endif
#if NENER>0
           do irad=1,nener
              e(irad)=e(irad)+se0(irad)*dtdx
           end do
#endif
           ! Right state at left interface
           qp(i,j,k,ir,1) = r - drx
           qp(i,j,k,iu,1) = u - dux
           qp(i,j,k,iv,1) = v - dvx
           qp(i,j,k,iw,1) = w - dwx
           qp(i,j,k,ip,1) = p - dpx
           if(qp(i,j,k,ir,1)<smallr)qp(i,j,k,ir,1)=q(i,j,k,ir)
#ifdef MHD
           qp(i,j,k,iA,1) = AL
           qp(i,j,k,iB,1) = B - dBx
           qp(i,j,k,iC,1) = C - dCx
#endif
#if NENER>0
           do irad=1,nener
              qp(i,j,k,ie+irad,1) = e(irad) - dex(irad)
           end do
#endif
           ! Left state at left interface
           qm(i,j,k,ir,1) = r + drx
           qm(i,j,k,iu,1) = u + dux
           qm(i,j,k,iv,1) = v + dvx
           qm(i,j,k,iw,1) = w + dwx
           qm(i,j,k,ip,1) = p + dpx
           if(qm(i,j,k,ir,1)<smallr)qm(i,j,k,ir,1)=q(i,j,k,ir)
#ifdef MHD
           qm(i,j,k,iA,1) = AR
           qm(i,j,k,iB,1) = B + dBx
           qm(i,j,k,iC,1) = C + dCx
#endif
#if NENER>0
           do irad=1,nener
              qm(i,j,k,ie+irad,1) = e(irad) + dex(irad)
           end do
#endif
           ! Top state at bottom interface
           qp(i,j,k,ir,2) = r - dry
           qp(i,j,k,iu,2) = u - duy
           qp(i,j,k,iv,2) = v - dvy
           qp(i,j,k,iw,2) = w - dwy
           qp(i,j,k,ip,2) = p - dpy
           if(qp(i,j,k,ir,2)<smallr)qp(i,j,k,ir,2)=q(i,j,k,ir)
#ifdef MHD
           qp(i,j,k,iA,2) = A - dAy
           qp(i,j,k,iB,2) = BL
           qp(i,j,k,iC,2) = C - dCy
#endif
#if NENER>0
           do irad=1,nener
              qp(i,j,k,ie+irad,2) = e(irad) - dey(irad)
           end do
#endif
           ! Bottom state at top interface
           qm(i,j,k,ir,2) = r + dry
           qm(i,j,k,iu,2) = u + duy
           qm(i,j,k,iv,2) = v + dvy
           qm(i,j,k,iw,2) = w + dwy
           qm(i,j,k,ip,2) = p + dpy
           if(qm(i,j,k,ir,2)<smallr)qm(i,j,k,ir,2)=q(i,j,k,ir)
#ifdef MHD
           qm(i,j,k,iA,2) = A + dAy
           qm(i,j,k,iB,2) = BR
           qm(i,j,k,iC,2) = C + dCy
#endif
#if NENER>0
           do irad=1,nener
              qm(i,j,k,ie+irad,2) = e(irad) + dey(irad)
           end do
#endif
           ! Back state at front interface
           qp(i,j,k,ir,3) = r - drz
           qp(i,j,k,iu,3) = u - duz
           qp(i,j,k,iv,3) = v - dvz
           qp(i,j,k,iw,3) = w - dwz
           qp(i,j,k,ip,3) = p - dpz
           if(qp(i,j,k,ir,3)<smallr)qp(i,j,k,ir,3)=q(i,j,k,ir)
#ifdef MHD
           qp(i,j,k,iA,3) = A - dAz
           qp(i,j,k,iB,3) = B - dBz
           qp(i,j,k,iC,3) = CL
#endif
#if NENER>0
           do irad=1,nener
              qp(i,j,k,ie+irad,3) = e(irad) - dez(irad)
           end do
#endif
           ! Front state at back interface
           qm(i,j,k,ir,3) = r + drz
           qm(i,j,k,iu,3) = u + duz
           qm(i,j,k,iv,3) = v + dvz
           qm(i,j,k,iw,3) = w + dwz
           qm(i,j,k,ip,3) = p + dpz
           if(qm(i,j,k,ir,3)<smallr)qm(i,j,k,ir,3)=q(i,j,k,ir)
#ifdef MHD
           qm(i,j,k,iA,3) = A + dAz
           qm(i,j,k,iB,3) = B + dBz
           qm(i,j,k,iC,3) = CR
#endif
#if NENER>0
           do irad=1,nener
              qm(i,j,k,ie+irad,3) = e(irad) + dez(irad)
           end do
#endif
#ifdef MHD
           ! X-edge averaged right-top corner state (RT->LL)
           qRT(i,j,k,ir,1) = r + (+dry+drz)
           qRT(i,j,k,iu,1) = u + (+duy+duz)
           qRT(i,j,k,iv,1) = v + (+dvy+dvz)
           qRT(i,j,k,iw,1) = w + (+dwy+dwz)
           qRT(i,j,k,ip,1) = p + (+dpy+dpz)
           qRT(i,j,k,iA,1) = A + (+dAy+dAz)
           qRT(i,j,k,iB,1) = BR+ (   +dBRz)
           qRT(i,j,k,iC,1) = CR+ (+dCRy   )
           if (qRT(i,j,k,ir,1)<smallr) qRT(i,j,k,ir,1)=r
#if NENER>0
           do irad=1,nener
              qRT(i,j,k,iC+irad,1) = e(irad) + (+dey(irad)+dez(irad))
           end do
#endif
           ! X-edge averaged right-bottom corner state (RB->LR)
           qRB(i,j,k,ir,1) = r + (+dry-drz)
           qRB(i,j,k,iu,1) = u + (+duy-duz)
           qRB(i,j,k,iv,1) = v + (+dvy-dvz)
           qRB(i,j,k,iw,1) = w + (+dwy-dwz)
           qRB(i,j,k,ip,1) = p + (+dpy-dpz)
           qRB(i,j,k,iA,1) = A + (+dAy-dAz)
           qRB(i,j,k,iB,1) = BR+ (   -dBRz)
           qRB(i,j,k,iC,1) = CL+ (+dCLy   )
           if (qRB(i,j,k,ir,1)<smallr) qRB(i,j,k,ir,1)=r
#if NENER>0
           do irad=1,nener
              qRB(i,j,k,iC+irad,1) = e(irad) + (+dey(irad)-dez(irad))
           end do
#endif
           ! X-edge averaged left-top corner state (LT->RL)
           qLT(i,j,k,ir,1) = r + (-dry+drz)
           qLT(i,j,k,iu,1) = u + (-duy+duz)
           qLT(i,j,k,iv,1) = v + (-dvy+dvz)
           qLT(i,j,k,iw,1) = w + (-dwy+dwz)
           qLT(i,j,k,ip,1) = p + (-dpy+dpz)
           qLT(i,j,k,iA,1) = A + (-dAy+dAz)
           qLT(i,j,k,iB,1) = BL+ (   +dBLz)
           qLT(i,j,k,iC,1) = CR+ (-dCRy   )
           if (qLT(i,j,k,ir,1)<smallr) qLT(i,j,k,ir,1)=r
#if NENER>0
           do irad=1,nener
              qLT(i,j,k,iC+irad,1) = e(irad) + (-dey(irad)+dez(irad))
           end do
#endif
           ! X-edge averaged left-bottom corner state (LB->RR)
           qLB(i,j,k,ir,1) = r + (-dry-drz)
           qLB(i,j,k,iu,1) = u + (-duy-duz)
           qLB(i,j,k,iv,1) = v + (-dvy-dvz)
           qLB(i,j,k,iw,1) = w + (-dwy-dwz)
           qLB(i,j,k,ip,1) = p + (-dpy-dpz)
           qLB(i,j,k,iA,1) = A + (-dAy-dAz)
           qLB(i,j,k,iB,1) = BL+ (   -dBLz)
           qLB(i,j,k,iC,1) = CL+ (-dCLy   )
           if (qLB(i,j,k,ir,1)<smallr) qLB(i,j,k,ir,1)=r
#if NENER>0
           do irad=1,nener
              qLB(i,j,k,iC+irad,1) = e(irad) + (-dey(irad)-dez(irad))
           end do
#endif
           ! Y-edge averaged right-top corner state (RT->LL)
           qRT(i,j,k,ir,2) = r + (+drx+drz)
           qRT(i,j,k,iu,2) = u + (+dux+duz)
           qRT(i,j,k,iv,2) = v + (+dvx+dvz)
           qRT(i,j,k,iw,2) = w + (+dwx+dwz)
           qRT(i,j,k,ip,2) = p + (+dpx+dpz)
           qRT(i,j,k,iA,2) = AR+ (   +dARz)
           qRT(i,j,k,iB,2) = B + (+dBx+dBz)
           qRT(i,j,k,iC,2) = CR+ (+dCRx   )
           if (qRT(i,j,k,ir,2)<smallr) qRT(i,j,k,ir,2)=r
#if NENER>0
           do irad=1,nener
              qRT(i,j,k,iC+irad,2) = e(irad) + (+dex(irad)+dez(irad))
           end do
#endif
           ! Y-edge averaged right-bottom corner state (RB->LR)
           qRB(i,j,k,ir,2) = r + (+drx-drz)
           qRB(i,j,k,iu,2) = u + (+dux-duz)
           qRB(i,j,k,iv,2) = v + (+dvx-dvz)
           qRB(i,j,k,iw,2) = w + (+dwx-dwz)
           qRB(i,j,k,ip,2) = p + (+dpx-dpz)
           qRB(i,j,k,iA,2) = AR+ (   -dARz)
           qRB(i,j,k,iB,2) = B + (+dBx-dBz)
           qRB(i,j,k,iC,2) = CL+ (+dCLx   )
           if (qRB(i,j,k,ir,2)<smallr) qRB(i,j,k,ir,2)=r
#if NENER>0
           do irad=1,nener
              qRB(i,j,k,iC+irad,2) = e(irad) + (+dex(irad)-dez(irad))
           end do
#endif
           ! Y-edge averaged left-top corner state (LT->RL)
           qLT(i,j,k,ir,2) = r + (-drx+drz)
           qLT(i,j,k,iu,2) = u + (-dux+duz)
           qLT(i,j,k,iv,2) = v + (-dvx+dvz)
           qLT(i,j,k,iw,2) = w + (-dwx+dwz)
           qLT(i,j,k,ip,2) = p + (-dpx+dpz)
           qLT(i,j,k,iA,2) = AL+ (   +dALz)
           qLT(i,j,k,iB,2) = B + (-dBx+dBz)
           qLT(i,j,k,iC,2) = CR+ (-dCRx   )
           if (qLT(i,j,k,ir,2)<smallr) qLT(i,j,k,ir,2)=r
#if NENER>0
           do irad=1,nener
              qLT(i,j,k,iC+irad,2) = e(irad) + (-dex(irad)+dez(irad))
           end do
#endif
           ! Y-edge averaged left-bottom corner state (LB->RR)
           qLB(i,j,k,ir,2) = r + (-drx-drz)
           qLB(i,j,k,iu,2) = u + (-dux-duz)
           qLB(i,j,k,iv,2) = v + (-dvx-dvz)
           qLB(i,j,k,iw,2) = w + (-dwx-dwz)
           qLB(i,j,k,ip,2) = p + (-dpx-dpz)
           qLB(i,j,k,iA,2) = AL+ (   -dALz)
           qLB(i,j,k,iB,2) = B + (-dBx-dBz)
           qLB(i,j,k,iC,2) = CL+ (-dCLx   )
           if (qLB(i,j,k,ir,2)<smallr) qLB(i,j,k,ir,2)=r
#if NENER>0
           do irad=1,nener
              qLB(i,j,k,iC+irad,2) = e(irad) + (-dex(irad)-dez(irad))
           end do
#endif
           ! Z-edge averaged right-top corner state (RT->LL)
           qRT(i,j,k,ir,3) = r + (+drx+dry)
           qRT(i,j,k,iu,3) = u + (+dux+duy)
           qRT(i,j,k,iv,3) = v + (+dvx+dvy)
           qRT(i,j,k,iw,3) = w + (+dwx+dwy)
           qRT(i,j,k,ip,3) = p + (+dpx+dpy)
           qRT(i,j,k,iA,3) = AR+ (   +dARy)
           qRT(i,j,k,iB,3) = BR+ (+dBRx   )
           qRT(i,j,k,iC,3) = C + (+dCx+dCy)
           if (qRT(i,j,k,ir,3)<smallr) qRT(i,j,k,ir,3)=r
#if NENER>0
           do irad=1,nener
              qRT(i,j,k,iC+irad,3) = e(irad) + (+dex(irad)+dey(irad))
           end do
#endif
           ! Z-edge averaged right-bottom corner state (RB->LR)
           qRB(i,j,k,ir,3) = r + (+drx-dry)
           qRB(i,j,k,iu,3) = u + (+dux-duy)
           qRB(i,j,k,iv,3) = v + (+dvx-dvy)
           qRB(i,j,k,iw,3) = w + (+dwx-dwy)
           qRB(i,j,k,ip,3) = p + (+dpx-dpy)
           qRB(i,j,k,iA,3) = AR+ (   -dARy)
           qRB(i,j,k,iB,3) = BL+ (+dBLx   )
           qRB(i,j,k,iC,3) = C + (+dCx-dCy)
           if (qRB(i,j,k,ir,3)<smallr) qRB(i,j,k,ir,3)=r
#if NENER>0
           do irad=1,nener
              qRB(i,j,k,iC+irad,3) = e(irad) + (+dex(irad)-dey(irad))
           end do
#endif
           ! Z-edge averaged left-top corner state (LT->RL)
           qLT(i,j,k,ir,3) = r + (-drx+dry)
           qLT(i,j,k,iu,3) = u + (-dux+duy)
           qLT(i,j,k,iv,3) = v + (-dvx+dvy)
           qLT(i,j,k,iw,3) = w + (-dwx+dwy)
           qLT(i,j,k,ip,3) = p + (-dpx+dpy)
           qLT(i,j,k,iA,3) = AL+ (   +dALy)
           qLT(i,j,k,iB,3) = BR+ (-dBRx   )
           qLT(i,j,k,iC,3) = C + (-dCx+dCy)
           if (qLT(i,j,k,ir,3)<smallr) qLT(i,j,k,ir,3)=r
#if NENER>0
           do irad=1,nener
              qLT(i,j,k,iC+irad,3) = e(irad) + (-dex(irad)+dey(irad))
           end do
#endif
           ! Z-edge averaged left-bottom corner state (LB->RR)
           qLB(i,j,k,ir,3) = r + (-drx-dry)
           qLB(i,j,k,iu,3) = u + (-dux-duy)
           qLB(i,j,k,iv,3) = v + (-dvx-dvy)
           qLB(i,j,k,iw,3) = w + (-dwx-dwy)
           qLB(i,j,k,ip,3) = p + (-dpx-dpy)
           qLB(i,j,k,iA,3) = AL+ (   -dALy)
           qLB(i,j,k,iB,3) = BL+ (-dBLx   )
           qLB(i,j,k,iC,3) = C + (-dCx-dCy)
           if (qLB(i,j,k,ir,3)<smallr) qLB(i,j,k,ir,3)=r
#if NENER>0
           do irad=1,nener
              qLB(i,j,k,iC+irad,3) = e(irad) + (-dex(irad)-dey(irad))
           end do
#endif
#endif
        end do
     end do
  end do

#if NVAR>5+NENER
  ! Passive scalars
  do n = ie+nener+1, nprim
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              r   = q(i,j,k,n)         ! Cell centered values
              u   = q(i,j,k,iu)
              v   = q(i,j,k,iv)
              w   = q(i,j,k,iw)
              drx = half*dq(i,j,k,n,1) ! TVD slopes
              dry = half*dq(i,j,k,n,2)
              drz = half*dq(i,j,k,n,3)
              sr0 = -u*drx-v*dry-w*drz ! Source terms
              r   = r + sr0*dtdx
              qp(i,j,k,n,1) = r - drx  ! Right state
              qm(i,j,k,n,1) = r + drx  ! Left state
              qp(i,j,k,n,2) = r - dry  ! Bottom state
              qm(i,j,k,n,2) = r + dry  ! Upper state
              qp(i,j,k,n,3) = r - drz  ! Front state
              qm(i,j,k,n,3) = r + drz  ! Back state
           end do
        end do
     end do
  end do
#endif

end subroutine trace3d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpflxm(qm,im1,im2,jm1,jm2,km1,km2, &
     &             qp,ip1,ip2,jp1,jp2,kp1,kp2, &
     &                ilo,ihi,jlo,jhi,klo,khi, &
#ifdef MHD
     &             bn,bt1,bt2, &
#endif
     &             ln,lt1,lt2,flx,tmp, &
     &             gamma,gamma_rad,smallr,smallc,riemann)
  use amr_parameters, only: dp,ndim
  use hydro_parameters
  use const
  implicit none
  ! routine arguments
  integer::ln,lt1,lt2
#ifdef MHD
  integer::bn,bt1,bt2
#endif
  integer::im1,im2,jm1,jm2,km1,km2
  integer::ip1,ip2,jp1,jp2,kp1,kp2
  integer::ilo,ihi,jlo,jhi,klo,khi
  real(dp)::gamma,smallr,smallc
  real(dp),dimension(1:nener+1)::gamma_rad
  integer::riemann
  real(dp),dimension(im1:im2,jm1:jm2,km1:km2,1:nprim,1:ndim)::qm
  real(dp),dimension(ip1:ip2,jp1:jp2,kp1:kp2,1:nprim,1:ndim)::qp
  real(dp),dimension(ip1:ip2,jp1:jp2,kp1:kp2,1:nprim)::flx
  real(dp),dimension(ip1:ip2,jp1:jp2,kp1:kp2,1:2)::tmp
  ! local variables
  integer ::i,j,k,xdim
  real(dp)::entho
  real(dp),dimension(1:nprim)::qleft,qright
  real(dp),dimension(1:nprim+1)::fgdnv
#if NVAR>5
  integer::n
#endif
  
  entho=one/(gamma-one)
  xdim=ln-1

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi

           ! Left state
           qleft (1) = qm(i,j,k,1  ,xdim) ! Mass density
           qleft (2) = qm(i,j,k,ln ,xdim) ! Normal velocity
           qleft (3) = qm(i,j,k,lt1,xdim) ! Tangential velocity 1
           qleft (4) = qm(i,j,k,lt2,xdim) ! Tangential velocity 2
           qleft (5) = qm(i,j,k,5  ,xdim) ! Pressure
#ifdef MHD
           qleft (6) = qm(i,j,k,bn ,xdim) ! Normal magnetic field
           qleft (7) = qm(i,j,k,bt1,xdim) ! Tangential magnetic field 1
           qleft (8) = qm(i,j,k,bt2,xdim) ! Tangential magnetic field 2
#endif
           ! Left state
           qright(1) = qp(i,j,k,1  ,xdim) ! Mass density
           qright(2) = qp(i,j,k,ln ,xdim) ! Normal velocity
           qright(3) = qp(i,j,k,lt1,xdim) ! Tangential velocity 1
           qright(4) = qp(i,j,k,lt2,xdim) ! Tangential velocity 2
           qright(5) = qp(i,j,k,5  ,xdim) ! Pressure
#ifdef MHD
           qright(6) = qp(i,j,k,bn ,xdim) ! Normal magnetic field
           qright(7) = qp(i,j,k,bt1,xdim) ! Tangential magnetic field 1
           qright(8) = qp(i,j,k,bt2,xdim) ! Tangential magnetic field 2
#endif
#if NVAR>5
           ! Other advected quantities
           do n = ie+1, nprim
              qleft (n) = qm(i,j,k,n,xdim)
              qright(n) = qp(i,j,k,n,xdim)
           end do
#endif
#ifndef MHD
           ! Solve hydro Riemann problem
           if(riemann.eq.solver_llf)then
              call riemann_llf(qleft,qright,fgdnv,gamma,gamma_rad,smallr,smallc)
           else if (riemann.eq.solver_hllc)then
              call riemann_hllc(qleft,qright,fgdnv,gamma,gamma_rad,smallr,smallc)
           else if (riemann.eq.solver_hll)then
              call riemann_hll(qleft,qright,fgdnv,gamma,gamma_rad,smallr,smallc)
#endif
#ifdef MHD
              ! Solve MHD Riemann problem
           if (riemann.eq.solver_llf)then
              call riemann_llf_mhd(qleft,qright,fgdnv,real(1.0,kind=dp),gamma,gamma_rad,smallr,smallc)
           else if (riemann.eq.solver_hlld)then
              call riemann_hlld(qleft,qright,fgdnv,gamma,gamma_rad,smallr,smallc)
           else if (riemann.eq.solver_hll)then
              call riemann_hll_mhd(qleft,qright,fgdnv,gamma,gamma_rad,smallr,smallc)
           else if (riemann.eq.solver_roe)then
              call riemann_roe_mhd(qleft,qright,fgdnv,real(1.0,kind=dp),gamma,gamma_rad,smallr,smallc)
           else if (riemann.eq.solver_upwind)then
              call riemann_upwind_mhd(qleft,qright,fgdnv,real(1.0,kind=dp),gamma,gamma_rad,smallr,smallc)
#endif
           else
              write(*,*)'unknown Riemann solver'
              stop
           end if
           
           ! Compute fluxes
           flx(i,j,k,1  ) = fgdnv(1) ! Mass density
           flx(i,j,k,ln ) = fgdnv(2) ! Normal momentum
           flx(i,j,k,lt1) = fgdnv(3) ! Transverse momentum 1
           flx(i,j,k,lt2) = fgdnv(4) ! Transverse momentum 2
           flx(i,j,k,5  ) = fgdnv(5) ! Total energy
#ifdef MHD
           flx(i,j,k,bn ) = fgdnv(6) ! Normal magnetic field
           flx(i,j,k,bt1) = fgdnv(7) ! Transverse magnetic field 1
           flx(i,j,k,bt2) = fgdnv(8) ! Transverse magnetic field 2
#endif
#if NVAR>5
           ! Other advected quantities
           do n = ie+1, nprim
              flx(i,j,k,n) = fgdnv(n)
           end do
#endif
           ! Normal velocity
           tmp(i,j,k,1) = half*(qleft(2)+qright(2))
           ! Internal energy flux
           tmp(i,j,k,2) = fgdnv(nprim+1)

        end do
     end do
  end do
  
end subroutine cmpflxm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#ifdef MHD
subroutine cmp_mag_flx(qRT,irt1,irt2,jrt1,jrt2,krt1,krt2, &
       &               qRB,irb1,irb2,jrb1,jrb2,krb1,krb2, &
       &               qLT,ilt1,ilt2,jlt1,jlt2,klt1,klt2, &
       &               qLB,ilb1,ilb2,jlb1,jlb2,klb1,klb2, &
       &                   ilo ,ihi ,jlo ,jhi ,klo ,khi , &
       &                   lp1 ,lp2 ,lor ,bp1 ,bp2 ,bor , &
       &               emf,gamma,gamma_rad,smallr,smallc,riemann2d)
  ! 2D Riemann solver to compute EMF at cell edges
  use amr_parameters, only: dp,ndim
  use hydro_parameters
  use const
  implicit none
  ! routine arguments
  ! indices of the 2 planar velocity lp1 and lp2
  ! and the orthogonal one, lor.
  ! Idem for the magnetic field.
  integer::lp1,lp2,lor,bp1,bp2,bor
  integer::irt1,irt2,jrt1,jrt2,krt1,krt2
  integer::irb1,irb2,jrb1,jrb2,krb1,krb2
  integer::ilt1,ilt2,jlt1,jlt2,klt1,klt2
  integer::ilb1,ilb2,jlb1,jlb2,klb1,klb2
  integer::ilo,ihi,jlo,jhi,klo,khi
  real(dp)::gamma,smallr,smallc
  real(dp),dimension(1:nener+1)::gamma_rad
  integer::riemann2d
  real(dp),dimension(irt1:irt2,jrt1:jrt2,krt1:krt2,1:nprim,1:3)::qRT
  real(dp),dimension(irb1:irb2,jrb1:jrb2,krb1:krb2,1:nprim,1:3)::qRB
  real(dp),dimension(ilt1:ilt2,jlt1:jlt2,klt1:klt2,1:nprim,1:3)::qLT
  real(dp),dimension(ilb1:ilb2,jlb1:jlb2,klb1:klb2,1:nprim,1:3)::qLB
  real(dp),dimension(ilb1:ilb2,jlb1:jlb2,klb1:klb2)::emf
  ! local variables
  integer::i, j, k, xdim
  real(dp),dimension(1:nprim)::qLL,qRL,qLR,qRR
  real(dp)::E

  xdim = lor - 1

  DO k = klo, khi
     DO j = jlo, jhi
        DO i = ilo, ihi

           ! Density
           qLL(1) = qRT(i,j,k,1,xdim)
           qRL(1) = qLT(i,j,k,1,xdim)
           qLR(1) = qRB(i,j,k,1,xdim)
           qRR(1) = qLB(i,j,k,1,xdim)

           ! First parallel velocity
           qLL(2) = qRT(i,j,k,lp1,xdim)
           qRL(2) = qLT(i,j,k,lp1,xdim)
           qLR(2) = qRB(i,j,k,lp1,xdim)
           qRR(2) = qLB(i,j,k,lp1,xdim)

           ! Second parallel velocity
           qLL(3) = qRT(i,j,k,lp2,xdim)
           qRL(3) = qLT(i,j,k,lp2,xdim)
           qLR(3) = qRB(i,j,k,lp2,xdim)
           qRR(3) = qLB(i,j,k,lp2,xdim)

           ! Orthogonal velocity
           qLL(4) = qRT(i,j,k,lor,xdim)
           qRL(4) = qLT(i,j,k,lor,xdim)
           qLR(4) = qRB(i,j,k,lor,xdim)
           qRR(4) = qLB(i,j,k,lor,xdim)

           ! Pressure
           qLL(5) = qRT(i,j,k,5,xdim)
           qRL(5) = qLT(i,j,k,5,xdim)
           qLR(5) = qRB(i,j,k,5,xdim)
           qRR(5) = qLB(i,j,k,5,xdim)

           ! First parallel magnetic field (enforce continuity)
           qLL(6) = half*(qRT(i,j,k,bp1,xdim)+qLT(i,j,k,bp1,xdim))
           qRL(6) = half*(qRT(i,j,k,bp1,xdim)+qLT(i,j,k,bp1,xdim))
           qLR(6) = half*(qRB(i,j,k,bp1,xdim)+qLB(i,j,k,bp1,xdim))
           qRR(6) = half*(qRB(i,j,k,bp1,xdim)+qLB(i,j,k,bp1,xdim))

           ! Second parallel magnetic field (enforce continuity)
           qLL(7) = half*(qRT(i,j,k,bp2,xdim)+qRB(i,j,k,bp2,xdim))
           qRL(7) = half*(qLT(i,j,k,bp2,xdim)+qLB(i,j,k,bp2,xdim))
           qLR(7) = half*(qRT(i,j,k,bp2,xdim)+qRB(i,j,k,bp2,xdim))
           qRR(7) = half*(qLT(i,j,k,bp2,xdim)+qLB(i,j,k,bp2,xdim))

           ! Orthogonal magnetic Field
           qLL(8) = qRT(i,j,k,bor,xdim)
           qRL(8) = qLT(i,j,k,bor,xdim)
           qLR(8) = qRB(i,j,k,bor,xdim)
           qRR(8) = qLB(i,j,k,bor,xdim)

#if NENER>0
           ! Non-thermal energies
           do irad = 1,nener
              qLL(8+irad) = qRT(i,j,k,8+irad,xdim)
              qRL(8+irad) = qLT(i,j,k,8+irad,xdim)
              qLR(8+irad) = qRB(i,j,k,8+irad,xdim)
              qRR(8+irad) = qLB(i,j,k,8+irad,xdim)
           end do
#endif
           ! Solve 2D Riemann problem
           if(riemann2d.eq.solver2d_hlld)then
              call riemann2d_hlld(qLL,qLR,qRL,qRR,E,gamma,gamma_rad,smallr,smallc)
           else if(riemann2d.eq.solver2d_hllf)then
              call riemann2d_hll_fast(qLL,qLR,qRL,qRR,E,gamma,gamma_rad,smallr,smallc)
           else if(riemann2d.eq.solver2d_hlla)then
              call riemann2d_hll_alfven(qLL,qLR,qRL,qRR,E,gamma,gamma_rad,smallr,smallc)
           else
              call riemann2d_simple(qLL,qLR,qRL,qRR,E,gamma,gamma_rad,smallr,smallc,riemann2d)
           endif

           ! Store EMF
           emf(i,j,k) = E

        end do
     end do
  end do

end subroutine cmp_mag_flx
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine ctoprim(uin,q,c,gravin, &
#ifdef MHD
     & bin,bf, &
#endif
     & dt,iu1,iu2,ju1,ju2,ku1,ku2,gamma,gamma_rad,smallr,smallc)
  use amr_parameters, only: dp, ndim
  use hydro_parameters, only: nvar, nprim, nener, ie
  use const
  implicit none

  real(dp)::dt
  integer::iu1,iu2,ju1,ju2,ku1,ku2
  real(dp)::smallr,smallc,gamma
  real(dp),dimension(1:nener+1)::gamma_rad
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
#ifdef MHD
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:6)::bin
  real(dp),dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf
#endif
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim)::q
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2)::c  

  real(dp)::eint, smalle, dtxhalf, oneoverrho
  real(dp)::ekin, erad, emag
  integer ::i, j, k
#if NENER>0
  integer ::irad
#endif
#if NVAR>5+NENER
  integer ::n
#endif

  smalle=smallc**2/gamma/(gamma-one)
  dtxhalf=dt*half

#ifdef MHD
  ! Compute face-centered manetic field from left and right face fields
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2+1
           if(i>iu1.and.i<=iu2)then
              bf(i,j,k,1) = half*(bin(i,j,k,1)+bin(i-1,j,k,4))
           else if(i==iu1)then
              bf(i,j,k,1) = bin(i,j,k,1)
           else
              bf(i,j,k,1) = bin(i-1,j,k,4)
           endif
        end do
     end do
  end do
  do k = ku1, ku2
     do j = ju1, ju2+1
        do i = iu1, iu2
           if(j>ju1.and.j<=ju2)then
              bf(i,j,k,2) = half*(bin(i,j,k,2)+bin(i,j-1,k,5))
           else if(j==ju1)then
              bf(i,j,k,2) = bin(i,j,k,2)
           else
              bf(i,j,k,2) = bin(i,j-1,k,5)
           endif
        end do
     end do
  end do
  do k = ku1, ku2+1
     do j = ju1, ju2
        do i = iu1, iu2
           if(k>ku1.and.k<=ku2)then
              bf(i,j,k,3) = half*(bin(i,j,k,3)+bin(i,j,k-1,6))
           else if(k==ku1)then
              bf(i,j,k,3) = bin(i,j,k,3)
           else
              bf(i,j,k,3) = bin(i,j,k-1,6)
           endif
        end do
     end do
  end do
#endif
  ! Convert to primitive variable
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2

           ! Compute density
           q(i,j,k,1) = max(uin(i,j,k,1),smallr)
           
           ! Compute velocities
           oneoverrho = one/q(i,j,k,1)
           q(i,j,k,2) = uin(i,j,k,2)*oneoverrho
           q(i,j,k,3) = uin(i,j,k,3)*oneoverrho
           q(i,j,k,4) = uin(i,j,k,4)*oneoverrho

           ! Compute specific kinetic energy
           ekin =        half*q(i,j,k,1)*q(i,j,k,2)**2
           ekin = ekin + half*q(i,j,k,1)*q(i,j,k,3)**2
           ekin = ekin + half*q(i,j,k,1)*q(i,j,k,4)**2

           emag = zero
#ifdef MHD
           ! Compute cell centered magnetic field
           q(i,j,k,6) = half*(bin(i,j,k,1)+bin(i,j,k,4))
           q(i,j,k,7) = half*(bin(i,j,k,2)+bin(i,j,k,5))
           q(i,j,k,8) = half*(bin(i,j,k,3)+bin(i,j,k,6))

           ! Compute magnetic energy
           emag = emag + half*q(i,j,k,6)**2
           emag = emag + half*q(i,j,k,7)**2
           emag = emag + half*q(i,j,k,8)**2
#endif
           ! Compute non-thermal pressure
           erad = zero
#if NENER>0
           do irad = 1,nener
              q(i,j,k,ie+irad) = (gamma_rad(irad)-one)*uin(i,j,k,5+irad)
              erad = erad + uin(i,j,k,5+irad)
           enddo
#endif
           ! Compute thermal pressure
           eint = MAX(uin(i,j,k,5)-ekin-erad-emag,smalle)
           q(i,j,k,5) = (gamma-one)*eint

           ! Compute thermal sound speed
           c(i,j,k)=gamma*q(i,j,k,5)
#if NENER>0
           do irad=1,nener
              c(i,j,k)=c(i,j,k)+gamma_rad(irad)*q(i,j,k,ie+irad)
           enddo
#endif
           c(i,j,k)=sqrt(c(i,j,k)*oneoverrho)
           
           ! Gravity predictor step
           q(i,j,k,2) = q(i,j,k,2) + gravin(i,j,k,1)*dtxhalf
#if NDIM>1
           q(i,j,k,3) = q(i,j,k,3) + gravin(i,j,k,2)*dtxhalf
#endif
#if NDIM>2
           q(i,j,k,4) = q(i,j,k,4) + gravin(i,j,k,3)*dtxhalf
#endif
        end do
     end do
  end do

#if NVAR>5+NENER
  ! Passive scalar
  do n = 1, nvar-5-nener
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              oneoverrho = one/q(i,j,k,1)
              q(i,j,k,ie+nener+n) = uin(i,j,k,5+nener+n)*oneoverrho
           end do
        end do
     end do
  end do
#endif
 
end subroutine ctoprim
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine uslope(q,dq, &
#ifdef MHD
     & bf,dbf, &
#endif
     & dx,dt,iu1,iu2,ju1,ju2,ku1,ku2,slope_type,slope_mag_type)
  use amr_parameters, only: dp, ndim
  use hydro_parameters, only: nprim
  use const
  implicit none
  ! routine arguments
  real(dp)::dx,dt
  integer::iu1,iu2,ju1,ju2,ku1,ku2
  integer::slope_type,slope_mag_type
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim)::q
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim,1:ndim)::dq
#ifdef MHD
  real(dp),dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf
  real(dp),dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3,1:2)::dbf
#endif
  ! local arrays
  integer::i, j, k,  n
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop
#if NDIM==2
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
#endif
#if NDIM==3
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
#endif
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff
  integer::ilo,ihi,jlo,jhi,klo,khi
  
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

#if NDIM==1
  if(slope_type==0)then
     dq = zero
  else
     do n = 1, nprim
        do k = klo, khi
        do j = jlo, jhi
        do i = ilo, ihi
           if(slope_type==1.or.slope_type==2.or.slope_type==3)then  ! minmod or average
              dlft = MIN(slope_type,2)*(q(i  ,j,k,n) - q(i-1,j,k,n))
              drgt = MIN(slope_type,2)*(q(i+1,j,k,n) - q(i  ,j,k,n))
              dcen = half*(dlft+drgt)/MIN(slope_type,2)
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dq(i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
           else if(slope_type==4)then ! superbee
              dcen = q(i,j,k,2)*dt/dx
              dlft = two/(one+dcen)*(q(i,j,k,n)-q(i-1,j,k,n))
              drgt = two/(one-dcen)*(q(i+1,j,k,n)-q(i,j,k,n))
              dcen = half*(q(i+1,j,k,n)-q(i-1,j,k,n))
              dsgn = sign(one, dlft)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dq(i,j,k,n,1) = dsgn*dlim !min(dlim,abs(dcen))
           else if(slope_type==5)then ! ultrabee
              if(n==1)then
                 dcen = q(i,j,k,2)*dt/dx
                 if(dcen>=0)then
                    dlft = two/(zero+dcen+1d-10)*(q(i,j,k,n)-q(i-1,j,k,n))
                    drgt = two/(one -dcen      )*(q(i+1,j,k,n)-q(i,j,k,n))
                 else
                    dlft = two/(one +dcen      )*(q(i,j,k,n)-q(i-1,j,k,n))
                    drgt = two/(zero-dcen+1d-10)*(q(i+1,j,k,n)-q(i,j,k,n))
                 endif
                 dsgn = sign(one, dlft)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 dcen = half*(q(i+1,j,k,n)-q(i-1,j,k,n))
                 if((dlft*drgt)<=zero)dlim=zero
                 dq(i,j,k,n,1) = dsgn*dlim !min(dlim,abs(dcen))
              else
                 dq(i,j,k,n,1) = 0.0
              end if
           else if(slope_type==6)then ! unstable
              if(n==1)then
                 dlft = (q(i,j,k,n)-q(i-1,j,k,n))
                 drgt = (q(i+1,j,k,n)-q(i,j,k,n))
                 slop = 0.5*(dlft+drgt)
                 dlim = slop
                 dq(i,j,k,n,1) = dlim
              else
                 dq(i,j,k,n,1) = 0.0
              end if
           else
              write(*,*)'Unknown slope type',dx,dt
              stop
           end if
        end do
        end do
        end do
     end do
  end if
#endif

#if NDIM==2              
  if(slope_type==0)then
     dq = zero
  else if(slope_type==1.or.slope_type==2)then  ! minmod or average
     do n = 1, nprim
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 dlft = slope_type*(q(i  ,j,k,n) - q(i-1,j,k,n))
                 drgt = slope_type*(q(i+1,j,k,n) - q(i  ,j,k,n))
                 dcen = half*(dlft+drgt)/slope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dq(i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 ! slopes in second coordinate direction
                 dlft = slope_type*(q(i,j  ,k,n) - q(i,j-1,k,n))
                 drgt = slope_type*(q(i,j+1,k,n) - q(i,j  ,k,n))
                 dcen = half*(dlft+drgt)/slope_type
                 dsgn = sign(one,dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dq(i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
              end do
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 2d unsplit slope
     do n = 1, nprim
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 dfll = q(i-1,j-1,k,n)-q(i,j,k,n)
                 dflm = q(i-1,j  ,k,n)-q(i,j,k,n)
                 dflr = q(i-1,j+1,k,n)-q(i,j,k,n)
                 dfml = q(i  ,j-1,k,n)-q(i,j,k,n)
                 dfmm = q(i  ,j  ,k,n)-q(i,j,k,n)
                 dfmr = q(i  ,j+1,k,n)-q(i,j,k,n)
                 dfrl = q(i+1,j-1,k,n)-q(i,j,k,n)
                 dfrm = q(i+1,j  ,k,n)-q(i,j,k,n)
                 dfrr = q(i+1,j+1,k,n)-q(i,j,k,n)
                 
                 vmin = min(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                 vmax = max(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                 
                 dfx  = half*(q(i+1,j,k,n)-q(i-1,j,k,n))
                 dfy  = half*(q(i,j+1,k,n)-q(i,j-1,k,n))
                 dff  = half*(abs(dfx)+abs(dfy))
                 
                 if(dff>zero)then
                    slop = min(one,min(abs(vmin),abs(vmax))/dff)
                 else
                    slop = one
                 endif
                 
                 dlim = slop
                 
                 dq(i,j,k,n,1) = dlim*dfx
                 dq(i,j,k,n,2) = dlim*dfy
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type',dx,dt
     stop
  endif

#ifdef MHD
  if (slope_mag_type==0) then
     dbf = zero
  else if (slope_mag_type==1 .or. slope_mag_type==2.or. slope_mag_type==3) then
     ! Face-centered Bx slope along direction Y
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi+1
              dlft = MIN(slope_mag_type,2)*(bf(i,j  ,k,1) - bf(i,j-1,k,1))
              drgt = MIN(slope_mag_type,2)*(bf(i,j+1,k,1) - bf(i,j  ,k,1))
              dcen = half*(dlft+drgt)/MIN(slope_mag_type,2)
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dbf(i,j,k,1,1) = dsgn*min(dlim,abs(dcen))
           enddo
        end do
     end do
     ! Face-centered By slope along direction X
     do k = klo, khi
        do j = jlo, jhi+1
           do i = ilo, ihi
              dlft = MIN(slope_mag_type,2)*(bf(i  ,j,k,2) - bf(i-1,j,k,2))
              drgt = MIN(slope_mag_type,2)*(bf(i+1,j,k,2) - bf(i  ,j,k,2))
              dcen = half*(dlft+drgt)/MIN(slope_mag_type,2)
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dbf(i,j,k,2,1) = dsgn*min(dlim,abs(dcen))
           enddo
        end do
     end do
  else
     write(*,*)'Unknown slope type mag',dx,dt
     stop
  endif
#endif
#endif

#if NDIM==3
  if(slope_type==0)then
     dq = zero
  else if(slope_type==1)then  ! minmod
     do n = 1, nprim
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 dlft = q(i  ,j,k,n) - q(i-1,j,k,n)
                 drgt = q(i+1,j,k,n) - q(i  ,j,k,n)
                 if((dlft*drgt)<=zero) then
                    dq(i,j,k,n,1) = zero
                 else if(dlft>0) then
                    dq(i,j,k,n,1) = min(dlft,drgt)
                 else
                    dq(i,j,k,n,1) = max(dlft,drgt)
                 end if
                 ! slopes in second coordinate direction
                 dlft = q(i,j  ,k,n) - q(i,j-1,k,n)
                 drgt = q(i,j+1,k,n) - q(i,j  ,k,n)
                 if((dlft*drgt)<=zero) then
                    dq(i,j,k,n,2) = zero
                 else if(dlft>0) then
                    dq(i,j,k,n,2) = min(dlft,drgt)
                 else
                    dq(i,j,k,n,2) = max(dlft,drgt)
                 end if
                 ! slopes in third coordinate direction
                 dlft = q(i,j,k  ,n) - q(i,j,k-1,n)
                 drgt = q(i,j,k+1,n) - q(i,j,k  ,n)
                 if((dlft*drgt)<=zero) then
                    dq(i,j,k,n,3) = zero
                 else if(dlft>0) then
                    dq(i,j,k,n,3) = min(dlft,drgt)
                 else
                    dq(i,j,k,n,3) = max(dlft,drgt)
                 end if
              end do
           end do
        end do
     end do
  else if(slope_type==2)then ! moncen
     do n = 1, nprim
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 dlft = slope_type*(q(i  ,j,k,n) - q(i-1,j,k,n))
                 drgt = slope_type*(q(i+1,j,k,n) - q(i  ,j,k,n))
                 dcen = half*(dlft+drgt)/slope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dq(i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 ! slopes in second coordinate direction
                 dlft = slope_type*(q(i,j  ,k,n) - q(i,j-1,k,n))
                 drgt = slope_type*(q(i,j+1,k,n) - q(i,j  ,k,n))
                 dcen = half*(dlft+drgt)/slope_type
                 dsgn = sign(one,dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dq(i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 ! slopes in third coordinate direction
                 dlft = slope_type*(q(i,j,k  ,n) - q(i,j,k-1,n))
                 drgt = slope_type*(q(i,j,k+1,n) - q(i,j,k  ,n))
                 dcen = half*(dlft+drgt)/slope_type
                 dsgn = sign(one,dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dq(i,j,k,n,3) = dsgn*min(dlim,abs(dcen))
              end do
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 3d unsplit slope
     do n = 1, nprim
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 dflll = q(i-1,j-1,k-1,n)-q(i,j,k,n)
                 dflml = q(i-1,j  ,k-1,n)-q(i,j,k,n)
                 dflrl = q(i-1,j+1,k-1,n)-q(i,j,k,n)
                 dfmll = q(i  ,j-1,k-1,n)-q(i,j,k,n)
                 dfmml = q(i  ,j  ,k-1,n)-q(i,j,k,n)
                 dfmrl = q(i  ,j+1,k-1,n)-q(i,j,k,n)
                 dfrll = q(i+1,j-1,k-1,n)-q(i,j,k,n)
                 dfrml = q(i+1,j  ,k-1,n)-q(i,j,k,n)
                 dfrrl = q(i+1,j+1,k-1,n)-q(i,j,k,n)
                 
                 dfllm = q(i-1,j-1,k  ,n)-q(i,j,k,n)
                 dflmm = q(i-1,j  ,k  ,n)-q(i,j,k,n)
                 dflrm = q(i-1,j+1,k  ,n)-q(i,j,k,n)
                 dfmlm = q(i  ,j-1,k  ,n)-q(i,j,k,n)
                 dfmmm = q(i  ,j  ,k  ,n)-q(i,j,k,n)
                 dfmrm = q(i  ,j+1,k  ,n)-q(i,j,k,n)
                 dfrlm = q(i+1,j-1,k  ,n)-q(i,j,k,n)
                 dfrmm = q(i+1,j  ,k  ,n)-q(i,j,k,n)
                 dfrrm = q(i+1,j+1,k  ,n)-q(i,j,k,n)
                 
                 dfllr = q(i-1,j-1,k+1,n)-q(i,j,k,n)
                 dflmr = q(i-1,j  ,k+1,n)-q(i,j,k,n)
                 dflrr = q(i-1,j+1,k+1,n)-q(i,j,k,n)
                 dfmlr = q(i  ,j-1,k+1,n)-q(i,j,k,n)
                 dfmmr = q(i  ,j  ,k+1,n)-q(i,j,k,n)
                 dfmrr = q(i  ,j+1,k+1,n)-q(i,j,k,n)
                 dfrlr = q(i+1,j-1,k+1,n)-q(i,j,k,n)
                 dfrmr = q(i+1,j  ,k+1,n)-q(i,j,k,n)
                 dfrrr = q(i+1,j+1,k+1,n)-q(i,j,k,n)
                 
                 vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                      &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                      &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                 vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                      &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                      &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                 
                 dfx  = half*(q(i+1,j,k,n)-q(i-1,j,k,n))
                 dfy  = half*(q(i,j+1,k,n)-q(i,j-1,k,n))
                 dfz  = half*(q(i,j,k+1,n)-q(i,j,k-1,n))
                 dff  = half*(abs(dfx)+abs(dfy)+abs(dfz))
                 
                 if(dff>zero)then
                    slop = min(one,min(abs(vmin),abs(vmax))/dff)
                 else
                    slop = one
                 endif
                 
                 dlim = slop
                 
                 dq(i,j,k,n,1) = dlim*dfx
                 dq(i,j,k,n,2) = dlim*dfy
                 dq(i,j,k,n,3) = dlim*dfz
                 
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type',dx,dt
     stop
  endif

#ifdef MHD
  if(slope_mag_type==0)then
     dbf = zero
  else if(slope_mag_type==1.or.slope_mag_type==2.or.slope_mag_type==3)then  ! minmod or moncen
     ! Face-centered Bx slopes along direction Y and Z
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi+1
              ! Slopes along first coordinate direction
              dlft = MIN(slope_mag_type,2)*(bf(i,j  ,k,1) - bf(i,j-1,k,1))
              drgt = MIN(slope_mag_type,2)*(bf(i,j+1,k,1) - bf(i,j  ,k,1))
              dcen = half*(dlft+drgt)/MIN(slope_mag_type,2)
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dbf(i,j,k,1,1) = dsgn*min(dlim,abs(dcen))
              ! Slopes along second coordinate direction
              dlft = MIN(slope_mag_type,2)*(bf(i,j,k  ,1) - bf(i,j,k-1,1))
              drgt = MIN(slope_mag_type,2)*(bf(i,j,k+1,1) - bf(i,j,k  ,1))
              dcen = half*(dlft+drgt)/MIN(slope_mag_type,2)
              dsgn = sign(one,dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dbf(i,j,k,1,2) = dsgn*min(dlim,abs(dcen))
           end do
        end do
     end do
     ! Face-centered By slopes along direction X and Z
     do k = klo, khi
        do j = jlo, jhi+1
           do i = ilo, ihi
              ! Slopes along first coordinate direction
              dlft = MIN(slope_mag_type,2)*(bf(i  ,j,k,2) - bf(i-1,j,k,2))
              drgt = MIN(slope_mag_type,2)*(bf(i+1,j,k,2) - bf(i  ,j,k,2))
              dcen = half*(dlft+drgt)/MIN(slope_mag_type,2)
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dbf(i,j,k,2,1) = dsgn*min(dlim,abs(dcen))
              ! Slopes along second coordinate direction
              dlft = MIN(slope_mag_type,2)*(bf(i,j,k  ,2) - bf(i,j,k-1,2))
              drgt = MIN(slope_mag_type,2)*(bf(i,j,k+1,2) - bf(i,j,k  ,2))
              dcen = half*(dlft+drgt)/MIN(slope_mag_type,2)
              dsgn = sign(one,dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dbf(i,j,k,2,2) = dsgn*min(dlim,abs(dcen))
           end do
        end do
     end do
     ! Face-cemtered Bz slopes along direction X and Y
     do k = klo, khi+1
        do j = jlo, jhi
           do i = ilo, ihi
              ! Slopes along first coordinate direction
              dlft = MIN(slope_mag_type,2)*(bf(i  ,j,k,3) - bf(i-1,j,k,3))
              drgt = MIN(slope_mag_type,2)*(bf(i+1,j,k,3) - bf(i  ,j,k,3))
              dcen = half*(dlft+drgt)/MIN(slope_mag_type,2)
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dbf(i,j,k,3,1) = dsgn*min(dlim,abs(dcen))
              ! Slopes along second coordinate direction
              dlft = MIN(slope_mag_type,2)*(bf(i,j  ,k,3) - bf(i,j-1,k,3))
              drgt = MIN(slope_mag_type,2)*(bf(i,j+1,k,3) - bf(i,j  ,k,3))
              dcen = half*(dlft+drgt)/MIN(slope_mag_type,2)
              dsgn = sign(one,dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dbf(i,j,k,3,2) = dsgn*min(dlim,abs(dcen))
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type mag',dx,dt
     stop
  endif
#endif
#endif
  
end subroutine uslope
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdivu(q,div,dx,dy,dz,&
     & iu1,iu2,ju1,ju2,ku1,ku2,&
     & if1,if2,jf1,jf2,kf1,kf2)
  use amr_parameters, only: dp, ndim
  use hydro_parameters, only: nvar, nprim
  use const
  implicit none
  ! This routine computes the velocity divergence at the corners of the cells.
  real(dp)::dx, dy, dz
  integer::iu1,iu2,ju1,ju2,ku1,ku2
  integer::if1,if2,jf1,jf2,kf1,kf2
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nprim)::q
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2)::div

  integer::i, j, k
  real(dp)::factorx, factory, factorz
  real(dp)::ux, vy, wz

  factorx=half**(ndim-1)/dx
  factory=half**(ndim-1)/dy
  factorz=half**(ndim-1)/dz

  do k = kf1, kf2
     do j = jf1, jf2
        do i = if1, if2
           ux=zero; vy=zero; wz=zero
           ux=ux+factorx*(q(i,j,k,2) - q(i-1,j,k,2))
#if NDIM>1
           ux=ux+factorx*(q(i  ,j-1,k,2) - q(i-1,j-1,k,2))
           vy=vy+factory*(q(i  ,j  ,k,3) - q(i  ,j-1,k,3)+&
                &         q(i-1,j  ,k,3) - q(i-1,j-1,k,3))
#endif
#if NDIM>2
           ux=ux+factorx*(q(i  ,j  ,k-1,2) - q(i-1,j  ,k-1,2)+&
                &         q(i  ,j-1,k-1,2) - q(i-1,j-1,k-1,2))
           vy=vy+factory*(q(i  ,j  ,k-1,3) - q(i  ,j-1,k-1,3)+&
                &         q(i-1,j  ,k-1,3) - q(i-1,j-1,k-1,3))
           wz=wz+factorz*(q(i  ,j  ,k  ,4) - q(i  ,j  ,k-1,4)+&
                &         q(i  ,j-1,k  ,4) - q(i  ,j-1,k-1,4)+&
                &         q(i-1,j  ,k  ,4) - q(i-1,j  ,k-1,4)+&
                &         q(i-1,j-1,k  ,4) - q(i-1,j-1,k-1,4))
#endif
           div(i,j,k) = ux + vy + wz

        end do
     end do
  end do

end subroutine cmpdivu
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdiff(uin,flux,div,dt,&
#ifdef MHD
     & bin, &
#endif
     & iu1,iu2,ju1,ju2,ku1,ku2,&
     & if1,if2,jf1,jf2,kf1,kf2,difmag)
  use amr_parameters, only: dp, ndim
  use hydro_parameters, only: nvar, nprim, ie
  use const
  implicit none
  ! Add diffusive flux where flow is compressing
  real(dp)::dt,difmag
  integer::iu1,iu2,ju1,ju2,ku1,ku2
  integer::if1,if2,jf1,jf2,kf1,kf2
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
#ifdef MHD
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:6)::bin
#endif
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2,1:nprim,1:ndim)::flux
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2)::div

  integer:: i, j, k, n
  real(dp)::factor
  real(dp)::div1

  factor=half**(ndim-1)

  do n = 1, nprim
     do k = kf1, MAX(kf1,ku2-2)
        do j = jf1, MAX(jf1, ju2-2)
           do i = if1, if2
              div1 = factor*div(i,j,k)
#if NDIM>1
              div1 = div1 + factor*div(i,j+1,k)
#endif
#if NDIM>2
              div1 = div1 + factor*(div(i,j,k+1)+div(i,j+1,k+1))
#endif
              div1 = difmag*min(zero,div1)
              if(n.LE.5)then
                 flux(i,j,k,n,1) = flux(i,j,k,n,1) + dt*div1*(uin(i,j,k,n) - uin(i-1,j,k,n))
#ifdef MHD
              else if(n.LE.8)then
                 flux(i,j,k,n,1) = flux(i,j,k,n,1) + dt*div1*(bin(i,j,k,n-5) - bin(i-1,j,k,n-5))
#endif
              else
                 flux(i,j,k,n,1) = flux(i,j,k,n,1) + dt*div1*(uin(i,j,k,n+5-ie) - uin(i-1,j,k,n+5-ie))
              endif
           end do
        end do
     end do

#if NDIM>1
     do k = kf1, MAX(kf1,ku2-2)
        do j = jf1, jf2
           do i = iu1+2, iu2-2
              div1 = zero
              div1 = div1 + factor*(div(i,j,k ) + div(i+1,j,k))
#if NDIM>2
              div1 = div1 + factor*(div(i,j,k+1) + div(i+1,j,k+1))
#endif
              div1 = difmag*min(zero,div1)
              if(n.LE.5)then
                 flux(i,j,k,n,2) = flux(i,j,k,n,2) + dt*div1*(uin(i,j,k,n) - uin(i,j-1,k,n))
#ifdef MHD
              else if(n.LE.8)then
                 flux(i,j,k,n,2) = flux(i,j,k,n,2) + dt*div1*(bin(i,j,k,n-5) - bin(i,j-1,k,n-5))
#endif
              else
                 flux(i,j,k,n,2) = flux(i,j,k,n,2) + dt*div1*(uin(i,j,k,n+5-ie) - uin(i,j-1,k,n+5-ie))
              end if
           end do
        end do
     end do
#endif

#if NDIM>2
     do k = kf1, kf2
        do j = ju1+2, ju2-2
           do i = iu1+2, iu2-2
              div1 = factor*(div(i,j  ,k) + div(i+1,j  ,k) + div(i,j+1,k) + div(i+1,j+1,k))
              div1 = difmag*min(zero,div1)
              if(n.LE.5)then
                 flux(i,j,k,n,3) = flux(i,j,k,n,3) + dt*div1*(uin(i,j,k,n) - uin(i,j,k-1,n))
#ifdef MHD
              else if(n.LE.8)then
                 flux(i,j,k,n,3) = flux(i,j,k,n,3) + dt*div1*(bin(i,j,k,n-5) - bin(i,j,k-1,n-5))
#endif
              else
                 flux(i,j,k,n,3) = flux(i,j,k,n,3) + dt*div1*(uin(i,j,k,n+5-ie) - uin(i,j,k-1,n+5-ie))
              end if
           end do
        end do
     end do
#endif

  end do

end subroutine cmpdiff
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpcurrent(bf,emfx,emfy,emfz,dx,dy,dz,dt,iu1,iu2,ju1,ju2,ku1,ku2,if1,if2,jf1,jf2,kf1,kf2,etamag)
  use amr_parameters,ONLY:dp
  implicit none
  integer::iu1,iu2,ju1,ju2,ku1,ku2
  integer::if1,if2,jf1,jf2,kf1,kf2
  real(dp),dimension(iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2)::emfx
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2)::emfy
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2)::emfz
  real(dp)::dx,dy,dz,dt,etamag
  ! Add to EMF -eta J where J = nabla x B
  real(dp)::dBx_arete_dy,dBx_arete_dz
  real(dp)::dBy_arete_dx,dBy_arete_dz
  real(dp)::dBz_arete_dx,dBz_arete_dy
  real(dp)::fact
  integer::i,j,k
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=iu1+2; ihi=iu2-2
  jlo=ju1+2; jhi=ju2-2
  klo=ku1+2; khi=ku2-2
  fact=etamag*dt/dx/dx

  ! Edge along x-axis
  do k=kf1,kf2
     do j=jf1,jf2
        do i=ilo,ihi
           dBz_arete_dy=(bf(i,j,k,3)-bf(i,j-1,k,3))
           dBy_arete_dz=(bf(i,j,k,2)-bf(i,j,k-1,2))
           emfx(i,j,k)=emfx(i,j,k)-(dBz_arete_dy-dBy_arete_dz)*fact
        enddo
     enddo
  enddo

  ! Edge along y-axis
  do k=kf1,kf2
     do j=jlo,jhi
        do i=if1,if2
           dBx_arete_dz=(bf(i,j,k,1)-bf(i,j,k-1,1))
           dBz_arete_dx=(bf(i,j,k,3)-bf(i-1,j,k,3))
           emfy(i,j,k)=emfy(i,j,k)-(dBx_arete_dz-dBz_arete_dx)*fact
        enddo
     enddo
  enddo

  ! Edge along z-axis
  do k=klo,khi
     do j=jf1,jf2
        do i=if1,if2
           dBy_arete_dx=(bf(i,j,k,2)-bf(i-1,j,k,2))
           dBx_arete_dy=(bf(i,j,k,1)-bf(i,j-1,k,1))
           emfz(i,j,k)=emfz(i,j,k)-(dBy_arete_dx-dBx_arete_dy)*fact
        enddo
     enddo
  enddo

end subroutine cmpcurrent
!###########################################################
!###########################################################
!###########################################################
!###########################################################
