!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hydro_refine(r,ug,um,ud,ok)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t
  use hydro_parameters, only: nvar,nener
  use const
  implicit none
  ! dummy arguments
  type(run_t)::r
  real(dp)::ug(1:nvar)
  real(dp)::um(1:nvar)
  real(dp)::ud(1:nvar)
  logical ::ok
  
  integer::idim
#if NENER>0
  integer::irad
#endif
  real(dp)::eking,ekinm,ekind
  real(dp)::dg,dm,dd,pg,pm,pd,vg,vm,vd,cg,cm,cd,error
  
  ! Convert to primitive variables
  ug(1) = max(ug(1),r%smallr)
  um(1) = max(um(1),r%smallr)
  ud(1) = max(ud(1),r%smallr)

  ! Velocity
  do idim = 1,ndim
     ug(idim+1) = ug(idim+1)/ug(1)
     um(idim+1) = um(idim+1)/um(1)
     ud(idim+1) = ud(idim+1)/ud(1)
  end do

  ! Pressure
  eking = zero
  ekinm = zero
  ekind = zero
  do idim = 1,ndim
     eking = eking + half*ug(1)*ug(idim+1)**2
     ekinm = ekinm + half*um(1)*um(idim+1)**2
     ekind = ekind + half*ud(1)*ud(idim+1)**2
  end do
#if NENER>0
  do irad = 1,nener
     eking = eking + ug(ndim+2+irad)
     ekinm = ekinm + um(ndim+2+irad)
     ekind = ekind + ud(ndim+2+irad)
  end do
#endif
  ug(ndim+2) = (r%gamma-one)*(ug(ndim+2)-eking)
  um(ndim+2) = (r%gamma-one)*(um(ndim+2)-ekinm)
  ud(ndim+2) = (r%gamma-one)*(ud(ndim+2)-ekind)

  ! Passive scalars
#if NVAR>NDIM + 2
  do idim = ndim+3,nvar
     ug(idim) = ug(idim)/ug(1)
     um(idim) = um(idim)/um(1)
     ud(idim) = ud(idim)/ud(1)
  end do
#endif

  ! Compute errors
  if(r%err_grad_d >= 0.)then
     dg=ug(1); dm=um(1); dd=ud(1)
     error=2.0d0*MAX( &
          & ABS((dd-dm)/(dd+dm+r%floor_d)) , &
          & ABS((dm-dg)/(dm+dg+r%floor_d)) )
     ok = ok .or. error > r%err_grad_d
  end if

  if(r%err_grad_p >= 0.)then
     pg=ug(ndim+2); pm=um(ndim+2); pd=ud(ndim+2)
     error=2.0d0*MAX( &
          & ABS((pd-pm)/(pd+pm+r%floor_p)), &
          & ABS((pm-pg)/(pm+pg+r%floor_p)) )
     ok = ok .or. error > r%err_grad_p
  end if

  if(r%err_grad_u >= 0.)then
     do idim = 1,ndim
        vg=ug(idim+1); vm=um(idim+1); vd=ud(idim+1)
        cg=sqrt(max(r%gamma*ug(ndim+2)/ug(1),r%floor_u**2))
        cm=sqrt(max(r%gamma*um(ndim+2)/um(1),r%floor_u**2))
        cd=sqrt(max(r%gamma*ud(ndim+2)/ud(1),r%floor_u**2))
        error=2.0d0*MAX( &
             & ABS((vd-vm)/(cd+cm+ABS(vd)+ABS(vm)+r%floor_u)) , &
             & ABS((vm-vg)/(cm+cg+ABS(vm)+ABS(vg)+r%floor_u)) )
        ok = ok .or. error > r%err_grad_u
     end do
  end if

end subroutine hydro_refine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdt(r,uu,gg,dx,dt)
  use amr_parameters, only: dp, ndim
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t
  use const
  implicit none
  type(run_t)::r
  real(dp)::dx,dt
  real(dp),dimension(1:nvar)::uu
  real(dp),dimension(1:ndim)::gg
  
  real(dp)::dtcell,smallp
  integer::idim
#if NENER>0
  integer::irad
#endif
  
  smallp = r%smallc**2/r%gamma

  ! Convert to primitive variables
  uu(1)=max(uu(1),r%smallr)
  ! Velocity
  do idim = 1,ndim
     uu(idim+1) = uu(idim+1)/uu(1)
  end do
  ! Internal energy
  do idim = 1,ndim
     uu(ndim+2) = uu(ndim+2)-half*uu(1)*uu(idim+1)**2
  end do
#if NENER>0
  do irad = 1,nener
     uu(ndim+2) = uu(ndim+2)-uu(ndim+2+irad)
  end do
#endif

  ! Debug
  if(r%debug)then
     if(uu(ndim+2).le.0.0D0.or.uu(1).le.r%smallr)then
        write(*,*)'stop in cmpdt'
        write(*,*)'dx   =',dx
        write(*,*)'rho  =',uu(1)
        write(*,*)'P    =',uu(ndim+2)
        write(*,*)'vel  =',uu(2:ndim+1)
        stop
     end if
  end if

  ! Compute pressure
  uu(ndim+2) = max((r%gamma-one)*uu(ndim+2),uu(1)*smallp)
#if NENER>0
  do irad = 1,nener
     uu(ndim+2+irad) = (r%gamma_rad(irad)-one)*uu(ndim+2+irad)
  end do
#endif

  ! Compute sound speed
  uu(ndim+2) = r%gamma*uu(ndim+2)
#if NENER>0
  do irad = 1,nener
     uu(ndim+2) = uu(ndim+2) + r%gamma_rad(irad)*uu(ndim+2+irad)
  end do
#endif  
  uu(ndim+2)=sqrt(uu(ndim+2)/uu(1))

  ! Compute wave speed
  uu(ndim+2)=dble(ndim)*uu(ndim+2)
  do idim = 1,ndim
     uu(ndim+2)=uu(ndim+2)+abs(uu(idim+1))
  end do

  ! Compute gravity strength ratio
  uu(1)=zero
  do idim = 1,ndim
     uu(1)=uu(1)+abs(gg(idim))
  end do
  uu(1)=uu(1)*dx/uu(ndim+2)**2
  uu(1)=MAX(uu(1),0.0001_dp)

  ! Compute maximum time step for each authorized cell
  dt = r%courant_factor*dx/r%smallc
  dtcell = dx/uu(ndim+2)*(sqrt(one+two*r%courant_factor*uu(1))-one)/uu(1)
  dt = min(dt,dtcell)

end subroutine cmpdt
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_llf(qleft,qright,fgdnv,gamma,gamma_rad,smallr,smallc)
  use amr_parameters, only: dp, ndim
  use hydro_parameters, only: nvar, nener
  use const
  implicit none

  ! dummy arguments
  real(dp)::gamma,smallr,smallc
  real(dp),dimension(1:nener)::gamma_rad
  real(dp),dimension(1:nvar)::qleft,qright
  real(dp),dimension(1:nvar+1)::fgdnv

  real(dp),dimension(1:nvar+1)::fleft,fright
  real(dp),dimension(1:nvar+1)::uleft,uright
  real(dp)::cmax
  integer::n
  real(dp)::smallp, entho
  real(dp)::rl   ,ul   ,pl   ,cl
  real(dp)::rr   ,ur   ,pr   ,cr   

  ! Constants
  smallp = smallc**2/gamma
  entho = one/(gamma-one)

  !===========================
  ! Compute maximum wave speed
  !===========================
  ! Left states
  rl = max(qleft(1),smallr)
  ul =     qleft(2)
  pl = max(qleft(3),rl*smallp)
  cl = gamma*pl
#if NENER>0
  do n = 1,nener
     cl = cl + gamma_rad(n)*qleft(ndim+2+n)
  end do
#endif
  cl = sqrt(cl/rl)
  ! Right states
  rr = max(qright(1),smallr)
  ur =     qright(2)
  pr = max(qright(3),rr*smallp)
  cr = gamma*pr
#if NENER>0
  do n = 1,nener
     cr = cr + gamma_rad(n)*qright(ndim+2+n)
  end do
#endif
  cr = sqrt(cr/rr)
  ! Local max. wave speed
  cmax = max(abs(ul)+cl,abs(ur)+cr)

  !===============================
  ! Compute conservative variables  
  !===============================
  ! Mass density
  uleft (1) = qleft (1)
  uright(1) = qright(1)
  ! Normal momentum
  uleft (2) = qleft (1)*qleft (2)
  uright(2) = qright(1)*qright(2)
  ! Total energy
  uleft (3) = qleft (3)*entho + half*qleft (1)*qleft (2)**2
  uright(3) = qright(3)*entho + half*qright(1)*qright(2)**2
#if NDIM>1
  uleft (3) = uleft (3)       + half*qleft (1)*qleft (4)**2
  uright(3) = uright(3)       + half*qright(1)*qright(4)**2
#endif
#if NDIM>2
  uleft (3) = uleft (3)       + half*qleft (1)*qleft (5)**2
  uright(3) = uright(3)       + half*qright(1)*qright(5)**2
#endif
#if NENER>0
  do n = 1,nener
     uleft (3) = uleft (3) + qleft (ndim+2+n)/(gamma_rad(n)-one)
     uright(3) = uright(3) + qright(ndim+2+n)/(gamma_rad(n)-one)
  end do
#endif
  ! Transverse velocities
#if NDIM>1
  do n = 4, ndim+2
     uleft (n) = qleft (1)*qleft (n)
     uright(n) = qright(1)*qright(n)
  end do
#endif
  ! Non-thermal energies
#if NENER>0
  do n = 1, nener
     uleft (ndim+2+n) = qleft (ndim+2+n)/(gamma_rad(n)-one)
     uright(ndim+2+n) = qright(ndim+2+n)/(gamma_rad(n)-one)
  end do
#endif
  ! Other passively advected quantities
#if NVAR>2+NDIM+NENER
  do n = 3+ndim+nener, nvar
     uleft (n) = qleft (1)*qleft (n)
     uright(n) = qright(1)*qright(n)
  end do
#endif
  ! Thermal energy
  uleft (nvar+1) = qleft (3)*entho
  uright(nvar+1) = qright(3)*entho

  !==============================
  ! Compute left and right fluxes  
  !==============================
  ! Mass density
  fleft (1) = qleft (2)*uleft (1)
  fright(1) = qright(2)*uright(1)
  ! Normal momentum
  fleft (2) = qleft (2)*uleft (2) + qleft (3)
  fright(2) = qright(2)*uright(2) + qright(3)
#if NENER>0
  do n = 1,nener
     fleft (2) = fleft (2) + qleft (ndim+2+n)
     fright(2) = fright(2) + qright(ndim+2+n)
  end do
#endif
  ! Total energy
  fleft (3) = qleft (2)*(uleft (3)+qleft (3))
  fright(3) = qright(2)*(uright(3)+qright(3))
#if NENER>0
  do n = 1,nener
     fleft (3) = fleft (3) + qleft (2)*qleft (ndim+2+n)
     fright(3) = fright(3) + qright(2)*qright(ndim+2+n)
  end do
#endif
  ! Other passively advected quantities
  do n = 4, nvar+1
     fleft (n) = qleft (2)*uleft (n)
     fright(n) = qright(2)*uright(n)
  end do

  !=============================
  ! Compute Lax-Friedrich fluxes
  !=============================
  do n = 1, nvar+1
     fgdnv(n) = half*(fleft(n)+fright(n)-cmax*(uright(n)-uleft(n)))
  end do

end subroutine riemann_llf
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_hll(qleft,qright,fgdnv,gamma,gamma_rad,smallr,smallc)
  use amr_parameters, only: dp, ndim
  use hydro_parameters, only: nvar, nener
  use const
  implicit none
  ! 1D HLL Riemann solver

  ! dummy arguments
  real(dp)::gamma,smallr,smallc
  real(dp),dimension(1:nener)::gamma_rad
  real(dp),dimension(1:nvar)::qleft,qright
  real(dp),dimension(1:nvar+1)::fgdnv

  real(dp),dimension(1:nvar+1)::fleft,fright
  real(dp),dimension(1:nvar+1)::uleft,uright
  real(dp)::SL,SR
  integer::n
  real(dp)::smallp, entho
  real(dp)::rl   ,ul   ,pl   ,cl
  real(dp)::rr   ,ur   ,pr   ,cr   

  ! Constants
  smallp = smallc**2/gamma
  entho = one/(gamma-one)

  !===========================
  ! Compute maximum wave speed
  !===========================
  ! Left states
  rl = max(qleft (1),smallr)
  ul =     qleft (2)
  pl = max(qleft (3),rl*smallp)
  cl = gamma*pl
#if NENER>0
  do n = 1,nener
     cl = cl + gamma_rad(n)*qleft(ndim+2+n)
  end do
#endif
  cl = sqrt(cl/rl)
  ! Right states
  rr = max(qright(1),smallr)
  ur =     qright(2)
  pr = max(qright(3),rr*smallp)
  cr = gamma*pr
#if NENER>0
  do n = 1,nener
     cr = cr + gamma_rad(n)*qright(ndim+2+n)
  end do
#endif
  cr = sqrt(cr/rr)
  ! Left and right max. wave speed
  SL=min(min(ul,ur)-max(cl,cr),zero)
  SR=max(max(ul,ur)+max(cl,cr),zero)

  !===============================
  ! Compute conservative variables
  !===============================
  ! Mass density
  uleft (1) = qleft (1)
  uright(1) = qright(1)
  ! Normal momentum
  uleft (2) = qleft (1)*qleft (2)
  uright(2) = qright(1)*qright(2)
  ! Total energy
  uleft (3) = qleft (3)*entho + half*qleft (1)*qleft (2)**2
  uright(3) = qright(3)*entho + half*qright(1)*qright(2)**2
#if NDIM>1
  uleft (3) = uleft (3)       + half*qleft (1)*qleft (4)**2
  uright(3) = uright(3)       + half*qright(1)*qright(4)**2
#endif
#if NDIM>2
  uleft (3) = uleft (3)       + half*qleft (1)*qleft (5)**2
  uright(3) = uright(3)       + half*qright(1)*qright(5)**2
#endif
#if NENER>0
  do n = 1,nener
     uleft (3) = uleft (3) + qleft (ndim+2+n)/(gamma_rad(n)-one)
     uright(3) = uright(3) + qright(ndim+2+n)/(gamma_rad(n)-one)
  end do
#endif
  ! Transverse velocities
#if NDIM>1
  do n = 4, ndim+2
     uleft (n) = qleft (1)*qleft (n)
     uright(n) = qright(1)*qright(n)
  end do
#endif
  ! Non-thermal energies
#if NENER>0
  do n = 1, nener
     uleft (ndim+2+n) = qleft (ndim+2+n)/(gamma_rad(n)-one)
     uright(ndim+2+n) = qright(ndim+2+n)/(gamma_rad(n)-one)
  end do
#endif
  ! Other passively advected quantities
#if NVAR>2+NDIM+NENER
  do n = 3+ndim+nener, nvar
     uleft (n) = qleft (1)*qleft (n)
     uright(n) = qright(1)*qright(n)
  end do
#endif
  ! Thermal energy
  uleft (nvar+1) = qleft (3)*entho
  uright(nvar+1) = qright(3)*entho

  !==============================
  ! Compute left and right fluxes
  !==============================
  ! Mass density
  fleft (1) = uleft (2)
  fright(1) = uright(2)
  ! Normal momentum
  fleft (2) = qleft (3)+uleft (2)*qleft (2)
  fright(2) = qright(3)+uright(2)*qright(2)
#if NENER>0
  do n = 1,nener
     fleft (2) = fleft (2) + qleft (ndim+2+n)
     fright(2) = fright(2) + qright(ndim+2+n)
  end do
#endif
  ! Total energy
  fleft (3) = qleft (2)*(uleft (3)+qleft (3))
  fright(3) = qright(2)*(uright(3)+qright(3))
#if NENER>0
  do n = 1,nener
     fleft (3) = fleft (3) + qleft (2)*qleft (ndim+2+n)
     fright(3) = fright(3) + qright(2)*qright(ndim+2+n)
  end do
#endif
  ! Other advected quantities
  do n = 4, nvar+1
     fleft (n) = qleft (2)*uleft (n)
     fright(n) = qright(2)*uright(n)
  end do

  !===================
  ! Compute HLL fluxes
  !===================
  do n = 1, nvar+1
     fgdnv(n) = (SR*fleft(n)-SL*fright(n) &
          & + SR*SL*(uright(n)-uleft(n)))/(SR-SL)
  end do

end subroutine riemann_hll
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_hllc(qleft,qright,fgdnv,gamma,gamma_rad,smallr,smallc)
  use amr_parameters, only: dp, ndim
  use hydro_parameters, only: nvar, nener
  use const
  implicit none
  ! HLLC Riemann solver (Toro)

  ! dummy arguments
  real(dp)::gamma,smallr,smallc
  real(dp),dimension(1:nener)::gamma_rad
  real(dp),dimension(1:nvar)::qleft,qright
  real(dp),dimension(1:nvar+1)::fgdnv

  REAL(dp)::SL,SR
  REAL(dp)::entho
  REAL(dp)::rl,pl,ul,ecinl,etotl,el,ptotl
  REAL(dp)::rr,pr,ur,ecinr,etotr,er,ptotr
  REAL(dp)::cfastl,rcl,rstarl,estarl
  REAL(dp)::cfastr,rcr,rstarr,estarr
  REAL(dp)::etotstarl,etotstarr
  REAL(dp)::ustar,ptotstar
  REAL(dp)::ro,uo,ptoto,etoto,eo
  REAL(dp)::smallp
#if NENER>0
  INTEGER::irad
  REAL(dp),dimension(1:nener)::eradl,eradr,erado
  REAL(dp),dimension(1:nener)::eradstarl,eradstarr
#endif
  INTEGER::ivar

  ! constants
  smallp = smallc**2/gamma
  entho = one/(gamma-one)

  ! Left variables
  rl=max(qleft (1),smallr)
  Pl=max(qleft (3),rl*smallp)
  ul=    qleft (2)
  
  el=Pl*entho
  ecinl = half*rl*ul*ul
#if NDIM>1
  ecinl=ecinl+half*rl*qleft(4)**2
#endif
#if NDIM>2
  ecinl=ecinl+half*rl*qleft(5)**2
#endif
  etotl = el+ecinl
#if NENER>0
  do irad=1,nener
     eradl(irad)=qleft(2+ndim+irad)/(gamma_rad(irad)-one)
     etotl=etotl+eradl(irad)
  end do
#endif
  Ptotl = Pl
#if NENER>0
  do irad=1,nener
     Ptotl=Ptotl+qleft(2+ndim+irad)
  end do
#endif
  
  ! Right variables
  rr=max(qright(1),smallr)
  Pr=max(qright(3),rr*smallp)
  ur=    qright(2)
  
  er=Pr*entho
  ecinr = half*rr*ur*ur
#if NDIM>1
  ecinr=ecinr+half*rr*qright(4)**2
#endif
#if NDIM>2
  ecinr=ecinr+half*rr*qright(5)**2
#endif
  etotr = er+ecinr
#if NENER>0
  do irad=1,nener
     eradr(irad)=qright(2+ndim+irad)/(gamma_rad(irad)-one)
     etotr=etotr+eradr(irad)
  end do
#endif
  Ptotr = Pr
#if NENER>0
  do irad=1,nener
     Ptotr=Ptotr+qright(2+ndim+irad)
  end do
#endif
  
  ! Find the largest eigenvalues in the normal direction to the interface
  cfastl=gamma*Pl
#if NENER>0
  do irad = 1,nener
     cfastl = cfastl + gamma_rad(irad)*qleft(ndim+2+irad)
  end do
#endif
  cfastl=sqrt(max(cfastl/rl,smallc**2))
  
  cfastr=gamma*Pr
#if NENER>0
  do irad = 1,nener
     cfastr = cfastr + gamma_rad(irad)*qright(ndim+2+irad)
  end do
#endif
  cfastr=sqrt(max(cfastr/rr,smallc**2))
  
  ! Compute HLL wave speed
  SL=min(ul,ur)-max(cfastl,cfastr)
  SR=max(ul,ur)+max(cfastl,cfastr)
  
  ! Compute lagrangian sound speed
  rcl=rl*(ul-SL)
  rcr=rr*(SR-ur)
  
  ! Compute acoustic star state
  ustar   =(rcr*ur   +rcl*ul   +  (Ptotl-Ptotr))/(rcr+rcl)
  Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)
  
  ! Left star region variables
  rstarl=rl*(SL-ul)/(SL-ustar)
  etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar)/(SL-ustar)
  estarl=el*(SL-ul)/(SL-ustar)
#if NENER>0
  do irad = 1,nener
     eradstarl(irad)=eradl(irad)*(SL-ul)/(SL-ustar)
  end do
#endif
  
  ! Right star region variables
  rstarr=rr*(SR-ur)/(SR-ustar)
  etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar)/(SR-ustar)
  estarr=er*(SR-ur)/(SR-ustar)
#if NENER>0
  do irad = 1,nener
     eradstarr(irad)=eradr(irad)*(SR-ur)/(SR-ustar)
  end do
#endif
  
  ! Sample the solution at x/t=0
  if(SL>0d0)then
     ro=rl
     uo=ul
     Ptoto=Ptotl
     etoto=etotl
     eo=el
#if NENER>0
     do irad = 1,nener
        erado(irad)=eradl(irad)
     end do
#endif
  else if(ustar>0d0)then
     ro=rstarl
     uo=ustar
     Ptoto=Ptotstar
     etoto=etotstarl
     eo=estarl
#if NENER>0
     do irad = 1,nener
        erado(irad)=eradstarl(irad)
     end do
#endif
  else if (SR>0d0)then
     ro=rstarr
     uo=ustar
     Ptoto=Ptotstar
     etoto=etotstarr
     eo=estarr
#if NENER>0
     do irad = 1,nener
        erado(irad)=eradstarr(irad)
     end do
#endif
  else
     ro=rr
     uo=ur
     Ptoto=Ptotr
     etoto=etotr
     eo=er
#if NENER>0
     do irad = 1,nener
        erado(irad)=eradr(irad)
     end do
#endif
  end if
  
  !=========================
  ! Compute the Godunov flux
  !=========================
  fgdnv(1) = ro*uo
  fgdnv(2) = ro*uo*uo+Ptoto
  fgdnv(3) = (etoto+Ptoto)*uo
  ! Transverse velocities
#if NDIM>1
  do ivar = 4,ndim+2
     if(ustar>0)then
        fgdnv(ivar) = ro*uo*qleft (ivar)
     else
        fgdnv(ivar) = ro*uo*qright(ivar)
     endif
  end do
#endif
  ! Non-thermal energies
#if NENER>0
  do irad = 1,nener
     fgdnv(ndim+2+irad) = uo*erado(irad)
  end do
#endif
  ! Other passively advected quantities
#if NVAR>2+NDIM+NENER
  do ivar = 3+ndim+nener,nvar
     if(ustar>0)then
        fgdnv(ivar) = ro*uo*qleft (ivar)
     else
        fgdnv(ivar) = ro*uo*qright(ivar)
     endif
  end do
#endif
  ! Thermal energy
  fgdnv(nvar+1) = uo*eo
  
end subroutine riemann_hllc
!###########################################################
!###########################################################
!###########################################################
!###########################################################
  


