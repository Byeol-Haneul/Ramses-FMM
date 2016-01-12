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
subroutine unsplit(uin,gravin,flux,tmp,dx,dy,dz,dt)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin 
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin 

  ! Output fluxes
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2,1:2   ,1:ndim)::tmp 

  ! Primitive variables
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::qin 
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2       ),save::cin

  ! Slopes
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::dq

  ! Left and right state arrays
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qm
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qp
  
  ! Intermediate fluxes
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::fx
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:2   ),save::tx

  ! Velocity divergence
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2)::divu

  ! Local scalar variables
  integer::i,j,k,l,ivar
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  ! Translate to primative variables, compute sound speeds  
  call ctoprim(uin,qin,cin,gravin,dt)

  ! Compute TVD slopes
  call uslope(qin,dq,dx,dt)

  ! Compute 3D traced-states in all three directions
#if NDIM==1
  call trace1d(qin,dq,qm,qp,dx      ,dt)
#endif
#if NDIM==2
  call trace2d(qin,dq,qm,qp,dx,dy   ,dt)
#endif
#if NDIM==3
  call trace3d(qin,dq,qm,qp,dx,dy,dz,dt)
#endif

  ! Solve for 1D flux in X direction
  call cmpflxm(qm,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          if1  ,if2  ,jlo  ,jhi  ,klo  ,khi  , 2,3,4,fx,tx)
  ! Save flux in output array
  do i=if1,if2
  do j=jlo,jhi
  do k=klo,khi
     do ivar=1,nvar
        flux(i,j,k,ivar,1)=fx(i,j,k,ivar)*dt/dx
     end do
     do ivar=1,2
        tmp (i,j,k,ivar,1)=tx(i,j,k,ivar)*dt/dx
     end do
  end do
  end do
  end do

  ! Solve for 1D flux in Y direction
#if NDIM>1
  call cmpflxm(qm,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jf1  ,jf2  ,klo  ,khi  , 3,2,4,fx,tx)
  ! Save flux in output array
  do i=ilo,ihi
  do j=jf1,jf2
  do k=klo,khi
     do ivar=1,nvar
        flux(i,j,k,ivar,2)=fx(i,j,k,ivar)*dt/dy
     end do
     do ivar=1,2
        tmp (i,j,k,ivar,2)=tx(i,j,k,ivar)*dt/dy
     end do
  end do
  end do
  end do
#endif

  ! Solve for 1D flux in Z direction
#if NDIM>2
  call cmpflxm(qm,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jlo  ,jhi  ,kf1  ,kf2  , 4,2,3,fx,tx)
  ! Save flux in output array
  do i=ilo,ihi
  do j=jlo,jhi
  do k=kf1,kf2
     do ivar=1,nvar
        flux(i,j,k,ivar,3)=fx(i,j,k,ivar)*dt/dz
     end do
     do ivar=1,2
        tmp (i,j,k,ivar,3)=tx(i,j,k,ivar)*dt/dz
     end do
  end do
  end do
  end do
#endif

  if(difmag>0.0)then
     call cmpdivu(qin,divu,dx,dy,dz)
     call consup(uin,flux,divu,dt)
  endif

end subroutine unsplit
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine trace1d(q,dq,qm,qp,dx,dt)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  real(dp)::dx, dt

  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  ! Local variables
  integer ::i, j, k,  n
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, ip, irad
  real(dp)::dtdx
  real(dp)::r, u, p, a
  real(dp)::drx, dux, dpx, dax
  real(dp)::sr0, su0, sp0, sa0
#if NENER>0
  real(dp),dimension(1:nener)::e, dex, se0
#endif
  
  dtdx = dt/dx

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; ip=3

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           
           ! Cell centered values
           r   =  q(i,j,k,ir)
           u   =  q(i,j,k,iu)
           p   =  q(i,j,k,ip)
#if NENER>0
           do irad=1,nener
              e(irad) = q(i,j,k,ip+irad)
           end do
#endif
           ! TVD slopes in X direction
           drx = dq(i,j,k,ir,1)
           dux = dq(i,j,k,iu,1)
           dpx = dq(i,j,k,ip,1)
#if NENER>0
           do irad=1,nener
              dex(irad) = dq(i,j,k,ip+irad,1)
           end do
#endif
           
           ! Source terms (including transverse derivatives)
           sr0 = -u*drx - (dux)*r
           sp0 = -u*dpx - (dux)*gamma*p
           su0 = -u*dux - (dpx)/r
#if NENER>0
           do irad=1,nener
              su0 = su0 - (dex(irad))/r
              se0(irad) = -u*dex(irad) &
                   & - (dux)*gamma_rad(irad)*e(irad)
           end do
#endif
           
           ! Right state
           qp(i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
           qp(i,j,k,iu,1) = u - half*dux + su0*dtdx*half
           qp(i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
           !              qp(i,j,k,ir,1) = max(smallr, qp(i,j,k,ir,1))
           if(qp(i,j,k,ir,1)<smallr)qp(i,j,k,ir,1)=r
#if NENER>0
           do irad=1,nener
              qp(i,j,k,ip+irad,1) = e(irad) - half*dex(irad) + se0(irad)*dtdx*half
           end do
#endif
           
           ! Left state
           qm(i,j,k,ir,1) = r + half*drx + sr0*dtdx*half
           qm(i,j,k,iu,1) = u + half*dux + su0*dtdx*half
           qm(i,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
           !              qm(i,j,k,ir,1) = max(smallr, qm(i,j,k,ir,1))
           if(qm(i,j,k,ir,1)<smallr)qm(i,j,k,ir,1)=r
#if NENER>0
           do irad=1,nener
              qm(i,j,k,ip+irad,1) = e(irad) + half*dex(irad) + se0(irad)*dtdx*half
           end do
#endif
           
        end do
     end do
  end do

#if NVAR > NDIM + 2 + NENER
  ! Passive scalars
  do n = ndim+nener+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              a   = q(i,j,k,n)       ! Cell centered values
              u   = q(i,j,k,iu)
              dax = dq(i,j,k,n,1)    ! TVD slopes
              sa0 = -u*dax             ! Source terms
              qp(i,j,k,n,1) = a - half*dax + sa0*dtdx*half   ! Right state
              qm(i,j,k,n,1) = a + half*dax + sa0*dtdx*half   ! Left state
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
subroutine trace2d(q,dq,qm,qp,dx,dy,dt)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  real(dp)::dx, dy, dt

  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  ! declare local variables
  integer ::i, j, k,  n
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, iv, ip, irad
  real(dp)::dtdx, dtdy
  real(dp)::r, u, v, p, a
  real(dp)::drx, dux, dvx, dpx, dax
  real(dp)::dry, duy, dvy, dpy, day
  real(dp)::sr0, su0, sv0, sp0, sa0
#if NENER>0
  real(dp),dimension(1:nener)::e, dex, dey, se0
#endif
  
  dtdx = dt/dx
  dtdy = dt/dy
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; ip=4

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           
           ! Cell centered values
           r   =  q(i,j,k,ir)
           u   =  q(i,j,k,iu)
           v   =  q(i,j,k,iv)
           p   =  q(i,j,k,ip)
#if NENER>0
           do irad=1,nener
              e(irad) = q(i,j,k,ip+irad)
           end do
#endif
           
           ! TVD slopes in all directions
           drx = dq(i,j,k,ir,1)
           dux = dq(i,j,k,iu,1)
           dvx = dq(i,j,k,iv,1)
           dpx = dq(i,j,k,ip,1)
#if NENER>0
           do irad=1,nener
              dex(irad) = dq(i,j,k,ip+irad,1)
           end do
#endif
           
           dry = dq(i,j,k,ir,2)
           duy = dq(i,j,k,iu,2)
           dvy = dq(i,j,k,iv,2)
           dpy = dq(i,j,k,ip,2)
#if NENER>0
           do irad=1,nener
              dey(irad) = dq(i,j,k,ip+irad,2)
           end do
#endif
           
           ! source terms (with transverse derivatives)
           sr0 = -u*drx-v*dry - (dux+dvy)*r
           sp0 = -u*dpx-v*dpy - (dux+dvy)*gamma*p
           su0 = -u*dux-v*duy - (dpx    )/r
           sv0 = -u*dvx-v*dvy - (dpy    )/r
#if NENER>0
           do irad=1,nener
              su0 = su0 - (dex(irad))/r
              sv0 = sv0 - (dey(irad))/r
              se0(irad) = -u*dex(irad)-v*dey(irad) &
                   & - (dux+dvy)*gamma_rad(irad)*e(irad)
           end do
#endif
           
           ! Right state at left interface
           qp(i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
           qp(i,j,k,iu,1) = u - half*dux + su0*dtdx*half
           qp(i,j,k,iv,1) = v - half*dvx + sv0*dtdx*half
           qp(i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
           !              qp(i,j,k,ir,1) = max(smallr, qp(i,j,k,ir,1))
           if(qp(i,j,k,ir,1)<smallr)qp(i,j,k,ir,1)=r
#if NENER>0
           do irad=1,nener
              qp(i,j,k,ip+irad,1) = e(irad) - half*dex(irad) + se0(irad)*dtdx*half
           end do
#endif
           
           ! Left state at right interface
           qm(i,j,k,ir,1) = r + half*drx + sr0*dtdx*half
           qm(i,j,k,iu,1) = u + half*dux + su0*dtdx*half
           qm(i,j,k,iv,1) = v + half*dvx + sv0*dtdx*half
           qm(i,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
           !              qm(i,j,k,ir,1) = max(smallr, qm(i,j,k,ir,1))
           if(qm(i,j,k,ir,1)<smallr)qm(i,j,k,ir,1)=r
#if NENER>0
           do irad=1,nener
              qm(i,j,k,ip+irad,1) = e(irad) + half*dex(irad) + se0(irad)*dtdx*half
           end do
#endif
           
           ! Top state at bottom interface
           qp(i,j,k,ir,2) = r - half*dry + sr0*dtdy*half
           qp(i,j,k,iu,2) = u - half*duy + su0*dtdy*half
           qp(i,j,k,iv,2) = v - half*dvy + sv0*dtdy*half
           qp(i,j,k,ip,2) = p - half*dpy + sp0*dtdy*half
           !              qp(i,j,k,ir,2) = max(smallr, qp(i,j,k,ir,2))
           if(qp(i,j,k,ir,2)<smallr)qp(i,j,k,ir,2)=r
#if NENER>0
           do irad=1,nener
              qp(i,j,k,ip+irad,2) = e(irad) - half*dey(irad) + se0(irad)*dtdy*half
           end do
#endif
           
           ! Bottom state at top interface
           qm(i,j,k,ir,2) = r + half*dry + sr0*dtdy*half
           qm(i,j,k,iu,2) = u + half*duy + su0*dtdy*half
           qm(i,j,k,iv,2) = v + half*dvy + sv0*dtdy*half
           qm(i,j,k,ip,2) = p + half*dpy + sp0*dtdy*half
           !              qm(i,j,k,ir,2) = max(smallr, qm(i,j,k,ir,2))
           if(qm(i,j,k,ir,2)<smallr)qm(i,j,k,ir,2)=r
#if NENER>0
           do irad=1,nener
              qm(i,j,k,ip+irad,2) = e(irad) + half*dey(irad) + se0(irad)*dtdy*half
           end do
#endif
           
        end do
     end do
  end do

#if NVAR > NDIM + 2 + NENER
  ! passive scalars
  do n = ndim+nener+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              a   = q(i,j,k,n)       ! Cell centered values
              u   = q(i,j,k,iu)
              v   = q(i,j,k,iv)
              dax = dq(i,j,k,n,1)    ! TVD slopes
              day = dq(i,j,k,n,2)
              sa0 = -u*dax-v*day       ! Source terms
              qp(i,j,k,n,1) = a - half*dax + sa0*dtdx*half   ! Right state
              qm(i,j,k,n,1) = a + half*dax + sa0*dtdx*half   ! Left state
              qp(i,j,k,n,2) = a - half*day + sa0*dtdy*half   ! Top state
              qm(i,j,k,n,2) = a + half*day + sa0*dtdy*half   ! Bottom state
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
subroutine trace3d(q,dq,qm,qp,dx,dy,dz,dt)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  real(dp)::dx, dy, dz, dt

  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  ! declare local variables
  integer ::i, j, k,  n
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, iv, iw, ip, irad
  real(dp)::dtdx, dtdy, dtdz
  real(dp)::r, u, v, w, p, a
  real(dp)::drx, dux, dvx, dwx, dpx, dax
  real(dp)::dry, duy, dvy, dwy, dpy, day
  real(dp)::drz, duz, dvz, dwz, dpz, daz
  real(dp)::sr0, su0, sv0, sw0, sp0, sa0
#if NENER>0
  real(dp),dimension(1:nener)::e, dex, dey, dez, se0
#endif
  
  dtdx = dt/dx
  dtdy = dt/dy
  dtdz = dt/dz
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4; ip=5

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           
           ! Cell centered values
           r   =  q(i,j,k,ir)
           u   =  q(i,j,k,iu)
           v   =  q(i,j,k,iv)
           w   =  q(i,j,k,iw)
           p   =  q(i,j,k,ip)
#if NENER>0
           do irad=1,nener
              e(irad) = q(i,j,k,ip+irad)
           end do
#endif
           
           ! TVD slopes in all 3 directions
           drx = dq(i,j,k,ir,1)
           dpx = dq(i,j,k,ip,1)
           dux = dq(i,j,k,iu,1)
           dvx = dq(i,j,k,iv,1)
           dwx = dq(i,j,k,iw,1)
#if NENER>0
           do irad=1,nener
              dex(irad) = dq(i,j,k,ip+irad,1)
           end do
#endif
           
           dry = dq(i,j,k,ir,2)
           dpy = dq(i,j,k,ip,2)
           duy = dq(i,j,k,iu,2)
           dvy = dq(i,j,k,iv,2)
           dwy = dq(i,j,k,iw,2)
#if NENER>0
           do irad=1,nener
              dey(irad) = dq(i,j,k,ip+irad,2)
           end do
#endif
           
           drz = dq(i,j,k,ir,3)
           dpz = dq(i,j,k,ip,3)
           duz = dq(i,j,k,iu,3)
           dvz = dq(i,j,k,iv,3)
           dwz = dq(i,j,k,iw,3)
#if NENER>0
           do irad=1,nener
              dez(irad) = dq(i,j,k,ip+irad,3)
           end do
#endif
           
           ! Source terms (including transverse derivatives)
           sr0 = -u*drx-v*dry-w*drz - (dux+dvy+dwz)*r
           sp0 = -u*dpx-v*dpy-w*dpz - (dux+dvy+dwz)*gamma*p
           su0 = -u*dux-v*duy-w*duz - (dpx        )/r
           sv0 = -u*dvx-v*dvy-w*dvz - (dpy        )/r
           sw0 = -u*dwx-v*dwy-w*dwz - (dpz        )/r
#if NENER>0
           do irad=1,nener
              su0 = su0 - (dex(irad))/r
              sv0 = sv0 - (dey(irad))/r
              sw0 = sw0 - (dez(irad))/r
              se0(irad) = -u*dex(irad)-v*dey(irad)-w*dez(irad) & 
                   & - (dux+dvy+dwz)*gamma_rad(irad)*e(irad)
           end do
#endif
           
           ! Right state at left interface
           qp(i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
           qp(i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
           qp(i,j,k,iu,1) = u - half*dux + su0*dtdx*half
           qp(i,j,k,iv,1) = v - half*dvx + sv0*dtdx*half
           qp(i,j,k,iw,1) = w - half*dwx + sw0*dtdx*half
           !              qp(i,j,k,ir,1) = max(smallr, qp(i,j,k,ir,1))
           if(qp(i,j,k,ir,1)<smallr)qp(i,j,k,ir,1)=r
#if NENER>0
           do irad=1,nener
              qp(i,j,k,ip+irad,1) = e(irad) - half*dex(irad) + se0(irad)*dtdx*half
           end do
#endif
           
           ! Left state at left interface
           qm(i,j,k,ir,1) = r + half*drx + sr0*dtdx*half
           qm(i,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
           qm(i,j,k,iu,1) = u + half*dux + su0*dtdx*half
           qm(i,j,k,iv,1) = v + half*dvx + sv0*dtdx*half
           qm(i,j,k,iw,1) = w + half*dwx + sw0*dtdx*half
           !              qm(i,j,k,ir,1) = max(smallr, qm(i,j,k,ir,1))
           if(qm(i,j,k,ir,1)<smallr)qm(i,j,k,ir,1)=r
#if NENER>0
           do irad=1,nener
              qm(i,j,k,ip+irad,1) = e(irad) + half*dex(irad) + se0(irad)*dtdx*half
           end do
#endif
           
           ! Top state at bottom interface
           qp(i,j,k,ir,2) = r - half*dry + sr0*dtdy*half
           qp(i,j,k,ip,2) = p - half*dpy + sp0*dtdy*half
           qp(i,j,k,iu,2) = u - half*duy + su0*dtdy*half
           qp(i,j,k,iv,2) = v - half*dvy + sv0*dtdy*half
           qp(i,j,k,iw,2) = w - half*dwy + sw0*dtdy*half
           !              qp(i,j,k,ir,2) = max(smallr, qp(i,j,k,ir,2))
           if(qp(i,j,k,ir,2)<smallr)qp(i,j,k,ir,2)=r
#if NENER>0
           do irad=1,nener
              qp(i,j,k,ip+irad,2) = e(irad) - half*dey(irad) + se0(irad)*dtdy*half
           end do
#endif
           
           ! Bottom state at top interface
           qm(i,j,k,ir,2) = r + half*dry + sr0*dtdy*half
           qm(i,j,k,ip,2) = p + half*dpy + sp0*dtdy*half
           qm(i,j,k,iu,2) = u + half*duy + su0*dtdy*half
           qm(i,j,k,iv,2) = v + half*dvy + sv0*dtdy*half
           qm(i,j,k,iw,2) = w + half*dwy + sw0*dtdy*half
           !              qm(i,j,k,ir,2) = max(smallr, qm(i,j,k,ir,2))
           if(qm(i,j,k,ir,2)<smallr)qm(i,j,k,ir,2)=r
#if NENER>0
           do irad=1,nener
              qm(i,j,k,ip+irad,2) = e(irad) + half*dey(irad) + se0(irad)*dtdy*half
           end do
#endif
           
           ! Back state at front interface
           qp(i,j,k,ir,3) = r - half*drz + sr0*dtdz*half
           qp(i,j,k,ip,3) = p - half*dpz + sp0*dtdz*half
           qp(i,j,k,iu,3) = u - half*duz + su0*dtdz*half
           qp(i,j,k,iv,3) = v - half*dvz + sv0*dtdz*half
           qp(i,j,k,iw,3) = w - half*dwz + sw0*dtdz*half
           !              qp(i,j,k,ir,3) = max(smallr, qp(i,j,k,ir,3))
           if(qp(i,j,k,ir,3)<smallr)qp(i,j,k,ir,3)=r
#if NENER>0
           do irad=1,nener
              qp(i,j,k,ip+irad,3) = e(irad) - half*dez(irad) + se0(irad)*dtdz*half
           end do
#endif
           
           ! Front state at back interface
           qm(i,j,k,ir,3) = r + half*drz + sr0*dtdz*half
           qm(i,j,k,ip,3) = p + half*dpz + sp0*dtdz*half
           qm(i,j,k,iu,3) = u + half*duz + su0*dtdz*half
           qm(i,j,k,iv,3) = v + half*dvz + sv0*dtdz*half
           qm(i,j,k,iw,3) = w + half*dwz + sw0*dtdz*half
           !              qm(i,j,k,ir,3) = max(smallr, qm(i,j,k,ir,3))
           if(qm(i,j,k,ir,3)<smallr)qm(i,j,k,ir,3)=r
#if NENER>0
           do irad=1,nener
              qm(i,j,k,ip+irad,3) = e(irad) + half*dez(irad) + se0(irad)*dtdz*half
           end do
#endif
           
        end do
     end do
  end do

#if NVAR > NDIM + 2 + NENER
  ! Passive scalars
  do n = ndim+nener+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              a   = q(i,j,k,n)       ! Cell centered values
              u   = q(i,j,k,iu)
              v   = q(i,j,k,iv)
              w   = q(i,j,k,iw)
              dax = dq(i,j,k,n,1)    ! TVD slopes
              day = dq(i,j,k,n,2)
              daz = dq(i,j,k,n,3)
              sa0 = -u*dax-v*day-w*daz     ! Source terms
              qp(i,j,k,n,1) = a - half*dax + sa0*dtdx*half  ! Right state
              qm(i,j,k,n,1) = a + half*dax + sa0*dtdx*half  ! Left state
              qp(i,j,k,n,2) = a - half*day + sa0*dtdy*half  ! Bottom state
              qm(i,j,k,n,2) = a + half*day + sa0*dtdy*half  ! Upper state
              qp(i,j,k,n,3) = a - half*daz + sa0*dtdz*half  ! Front state
              qm(i,j,k,n,3) = a + half*daz + sa0*dtdz*half  ! Back state
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
     &                ilo,ihi,jlo,jhi,klo,khi, ln,lt1,lt2, &
     &            flx,tmp)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ln,lt1,lt2
  integer ::im1,im2,jm1,jm2,km1,km2
  integer ::ip1,ip2,jp1,jp2,kp1,kp2
  integer ::ilo,ihi,jlo,jhi,klo,khi
  real(dp),dimension(im1:im2,jm1:jm2,km1:km2,1:nvar,1:ndim)::qm
  real(dp),dimension(ip1:ip2,jp1:jp2,kp1:kp2,1:nvar,1:ndim)::qp 
  real(dp),dimension(ip1:ip2,jp1:jp2,kp1:kp2,1:nvar)::flx
  real(dp),dimension(ip1:ip2,jp1:jp2,kp1:kp2,1:2)::tmp
  
  ! local variables
  integer ::i, j, k, n,  idim, xdim
  real(dp)::entho
  real(dp),dimension(1:nvar),save::qleft,qright
  real(dp),dimension(1:nvar+1),save::fgdnv

  entho=one/(gamma-one)
  xdim=ln-1

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           
           ! Mass density
           qleft (1) = qm(i,j,k,1,xdim)
           qright(1) = qp(i,j,k,1,xdim)
           
           ! Normal velocity
           qleft (2) = qm(i,j,k,ln,xdim)
           qright(2) = qp(i,j,k,ln,xdim)
           
           ! Pressure
           qleft (3) = qm(i,j,k,ndim+2,xdim)
           qright(3) = qp(i,j,k,ndim+2,xdim)
           
           ! Tangential velocity 1
#if NDIM>1
           qleft (4) = qm(i,j,k,lt1,xdim)
           qright(4) = qp(i,j,k,lt1,xdim)
#endif
           ! Tangential velocity 2
#if NDIM>2
           qleft (5) = qm(i,j,k,lt2,xdim)
           qright(5) = qp(i,j,k,lt2,xdim)
#endif           
#if NVAR > NDIM + 2
           ! Other advected quantities
           do n = ndim+3, nvar
              qleft (n) = qm(i,j,k,n,xdim)
              qright(n) = qp(i,j,k,n,xdim)
           end do
#endif          
           ! Solve Riemann problem
           if(riemann.eq.'llf')then
              call riemann_llf     (qleft,qright,fgdnv)
           else if (riemann.eq.'hllc')then
              call riemann_hllc    (qleft,qright,fgdnv)
           else if (riemann.eq.'hll')then
              call riemann_hll     (qleft,qright,fgdnv)
           else
              write(*,*)'unknown Riemann solver'
              stop
           end if
           
           ! Compute fluxes
           
           ! Mass density
           flx(i,j,k,1) = fgdnv(1)
           
           ! Normal momentum
           flx(i,j,k,ln) = fgdnv(2)

           ! Transverse momentum 1
#if NDIM>1
           flx(i,j,k,lt1) = fgdnv(4)
#endif
           ! Transverse momentum 2
#if NDIM>2
           flx(i,j,k,lt2) = fgdnv(5)
#endif           
           ! Total energy
           flx(i,j,k,ndim+2) = fgdnv(3)

#if NVAR > NDIM + 2
           ! Other advected quantities
           do n = ndim+3, nvar
              flx(i,j,k,n) = fgdnv(n)
           end do
#endif
           ! Normal velocity
           tmp(i,j,k,1) = half*(qleft(2)+qright(2))
           ! Internal energy flux
           tmp(i,j,k,2) = fgdnv(nvar+1)

        end do
     end do
  end do
  
end subroutine cmpflxm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine ctoprim(uin,q,c,gravin,dt)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  real(dp)::dt
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2)::c  

  integer ::i, j, k,  n, idim, irad
  real(dp)::eint, smalle, dtxhalf, oneoverrho
  real(dp)::eken, erad

  smalle = smallc**2/gamma/(gamma-one)
  dtxhalf = dt*half

  ! Convert to primitive variable
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2

           ! Compute density
           q(i,j,k,1) = max(uin(i,j,k,1),smallr)
           
           ! Compute velocities
           oneoverrho = one/q(i,j,k,1)
           q(i,j,k,2) = uin(i,j,k,2)*oneoverrho
#if NDIM>1
           q(i,j,k,3) = uin(i,j,k,3)*oneoverrho
#endif
#if NDIM>2
           q(i,j,k,4) = uin(i,j,k,4)*oneoverrho
#endif
           
           ! Compute specific kinetic energy
           eken = half*q(i,j,k,2)*q(i,j,k,2)
#if NDIM>1
           eken = eken + half*q(i,j,k,3)*q(i,j,k,3)
#endif
#if NDIM>2
           eken = eken + half*q(i,j,k,4)*q(i,j,k,4)
#endif
           ! Compute non-thermal pressure
           erad = zero
#if NENER>0
           do irad = 1,nener
              q(i,j,k,ndim+2+irad) = (gamma_rad(irad)-one)*uin(i,j,k,ndim+2+irad)
              erad = erad+uin(i,j,k,ndim+2+irad)*oneoverrho
           enddo
#endif
           ! Compute thermal pressure
           eint = MAX(uin(i,j,k,ndim+2)*oneoverrho-eken-erad,smalle)
           q(i,j,k,ndim+2) = (gamma-one)*q(i,j,k,1)*eint
           
           ! Compute sound speed
           c(i,j,k)=gamma*q(i,j,k,ndim+2)
#if NENER>0
           do irad=1,nener
              c(i,j,k)=c(i,j,k)+gamma_rad(irad)*q(i,j,k,ndim+2+irad)
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

#if NVAR > NDIM + 2 + NENER
  ! Passive scalar
  do n = ndim+nener+3, nvar
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              oneoverrho = one/q(i,j,k,1)
              q(i,j,k,n) = uin(i,j,k,n)*oneoverrho
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
subroutine uslope(q,dq,dx,dt)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  real(dp)::dx,dt
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq

  ! local arrays
  integer::i, j, k,  n
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff
  integer::ilo,ihi,jlo,jhi,klo,khi
  
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  if(slope_type==0)then
     dq=zero
     return
  end if

#if NDIM==1
  do n = 1, nvar
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
                 write(*,*)'Unknown slope type'
                 stop
              end if
           end do
        end do
     end do     
  end do
#endif

#if NDIM==2              
  if(slope_type==1.or.slope_type==2)then  ! minmod or average
     do n = 1, nvar
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
     do n = 1, nvar
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
     write(*,*)'Unknown slope type'
     stop
  endif
#endif

#if NDIM==3
  if(slope_type==1)then  ! minmod
     do n = 1, nvar
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
     do n = 1, nvar
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
     do n = 1, nvar
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
     write(*,*)'Unknown slope type'
     stop
  endif     
#endif
  
end subroutine uslope
