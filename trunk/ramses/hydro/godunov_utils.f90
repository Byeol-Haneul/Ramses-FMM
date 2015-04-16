!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdt(uu,gg,dx,dt,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  integer::ncell
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp),dimension(1:nvector,1:ndim)::gg
  
  real(dp)::dtcell,smallp
  integer::k,idim,irad
  
  smallp = smallc**2/gamma

  ! Convert to primitive variables
  do k = 1,ncell
     uu(k,1)=max(uu(k,1),smallr)
  end do
  ! Velocity
  do idim = 1,ndim
     do k = 1, ncell
        uu(k,idim+1) = uu(k,idim+1)/uu(k,1)
     end do
  end do
  ! Internal energy
  do idim = 1,ndim
     do k = 1, ncell
        uu(k,ndim+2) = uu(k,ndim+2)-half*uu(k,1)*uu(k,idim+1)**2
     end do
  end do
#if NENER>0
  do irad = 1,nener
     do k = 1, ncell
        uu(k,ndim+2) = uu(k,ndim+2)-uu(k,ndim+2+irad)
     end do
  end do
#endif

  ! Debug
  if(debug)then
     do k = 1, ncell 
        if(uu(k,ndim+2).le.0.or.uu(k,1).le.smallr)then
           write(*,*)'stop in cmpdt'
           write(*,*)'dx   =',dx
           write(*,*)'ncell=',ncell
           write(*,*)'rho  =',uu(k,1)
           write(*,*)'P    =',uu(k,ndim+2)
           write(*,*)'vel  =',uu(k,2:ndim+1)
           stop
        end if
     end do
  end if

  ! Compute pressure
  do k = 1, ncell
     uu(k,ndim+2) = max((gamma-one)*uu(k,ndim+2),uu(k,1)*smallp)
  end do
#if NENER>0
  do irad = 1,nener
     do k = 1, ncell
        uu(k,ndim+2+irad) = (gamma_rad(irad)-one)*uu(k,ndim+2+irad)
     end do
  end do
#endif

  ! Compute sound speed
  do k = 1, ncell
     uu(k,ndim+2) = gamma*uu(k,ndim+2)
  end do
#if NENER>0
  do irad = 1,nener
     do k = 1, ncell
        uu(k,ndim+2) = uu(k,ndim+2) + gamma_rad(irad)*uu(k,ndim+2+irad)
     end do
  end do
#endif  
  do k = 1, ncell 
     uu(k,ndim+2)=sqrt(uu(k,ndim+2)/uu(k,1))
  end do

  ! Compute wave speed
  do k = 1, ncell
     uu(k,ndim+2)=dble(ndim)*uu(k,ndim+2)
  end do
  do idim = 1,ndim
     do k = 1, ncell 
        uu(k,ndim+2)=uu(k,ndim+2)+abs(uu(k,idim+1))
     end do
  end do

  ! Compute gravity strength ratio
  do k = 1, ncell
     uu(k,1)=zero
  end do
  do idim = 1,ndim
     do k = 1, ncell 
        uu(k,1)=uu(k,1)+abs(gg(k,idim))
     end do
  end do
  do k = 1, ncell
     uu(k,1)=uu(k,1)*dx/uu(k,ndim+2)**2
     uu(k,1)=MAX(uu(k,1),0.0001_dp)
  end do

  ! Compute maximum time step for each authorized cell
  dt = courant_factor*dx/smallc
  do k = 1,ncell
     dtcell = dx/uu(k,ndim+2)*(sqrt(one+two*courant_factor*uu(k,1))-one)/uu(k,1)
     dt = min(dt,dtcell)
  end do

end subroutine cmpdt
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hydro_refine(ug,um,ud,ok,nn)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  ! dummy arguments
  integer nn
  real(dp)::ug(1:nvector,1:nvar)
  real(dp)::um(1:nvector,1:nvar)
  real(dp)::ud(1:nvector,1:nvar)
  logical ::ok(1:nvector)
  
  integer::k,idim,irad
  real(dp),dimension(1:nvector),save::eking,ekinm,ekind
  real(dp)::dg,dm,dd,pg,pm,pd,vg,vm,vd,cg,cm,cd,error
  
  ! Convert to primitive variables
  do k = 1,nn
     ug(k,1) = max(ug(k,1),smallr)
     um(k,1) = max(um(k,1),smallr)
     ud(k,1) = max(ud(k,1),smallr)
  end do
  ! Velocity
  do idim = 1,ndim
     do k = 1,nn
        ug(k,idim+1) = ug(k,idim+1)/ug(k,1)
        um(k,idim+1) = um(k,idim+1)/um(k,1)
        ud(k,idim+1) = ud(k,idim+1)/ud(k,1)
     end do
  end do
  ! Pressure
  do k = 1,nn
     eking(k) = zero
     ekinm(k) = zero
     ekind(k) = zero
  end do
  do idim = 1,ndim
     do k = 1,nn
        eking(k) = eking(k) + half*ug(k,1)*ug(k,idim+1)**2
        ekinm(k) = ekinm(k) + half*um(k,1)*um(k,idim+1)**2
        ekind(k) = ekind(k) + half*ud(k,1)*ud(k,idim+1)**2
     end do
  end do
#if NENER>0
  do irad = 1,nener
     do k = 1, nn
        eking(k) = eking(k) + ug(k,ndim+2+irad)
        ekinm(k) = ekinm(k) + um(k,ndim+2+irad)
        ekind(k) = ekind(k) + ud(k,ndim+2+irad)
     end do
  end do
#endif
  do k = 1,nn
     ug(k,ndim+2) = (gamma-one)*(ug(k,ndim+2)-eking(k))
     um(k,ndim+2) = (gamma-one)*(um(k,ndim+2)-ekinm(k))
     ud(k,ndim+2) = (gamma-one)*(ud(k,ndim+2)-ekind(k))
  end do
  ! Passive scalars
#if NVAR > NDIM + 2
  do idim = ndim+3,nvar
     do k = 1,nn
        ug(k,idim) = ug(k,idim)/ug(k,1)
        um(k,idim) = um(k,idim)/um(k,1)
        ud(k,idim) = ud(k,idim)/ud(k,1)
     end do
  end do
#endif

  ! Compute errors
  if(err_grad_d >= 0.)then
     do k=1,nn
        dg=ug(k,1); dm=um(k,1); dd=ud(k,1)
        error=2.0d0*MAX( &
             & ABS((dd-dm)/(dd+dm+floor_d)) , &
             & ABS((dm-dg)/(dm+dg+floor_d)) )
        ok(k) = ok(k) .or. error > err_grad_d
     end do
  end if

  if(err_grad_p >= 0.)then
     do k=1,nn
        pg=ug(k,ndim+2); pm=um(k,ndim+2); pd=ud(k,ndim+2)
        error=2.0d0*MAX( &
             & ABS((pd-pm)/(pd+pm+floor_p)), &
             & ABS((pm-pg)/(pm+pg+floor_p)) )
        ok(k) = ok(k) .or. error > err_grad_p
     end do
  end if

  if(err_grad_u >= 0.)then
     do idim = 1,ndim
        do k=1,nn
           vg=ug(k,idim+1); vm=um(k,idim+1); vd=ud(k,idim+1)
           cg=sqrt(max(gamma*ug(k,ndim+2)/ug(k,1),floor_u**2))
           cm=sqrt(max(gamma*um(k,ndim+2)/um(k,1),floor_u**2))
           cd=sqrt(max(gamma*ud(k,ndim+2)/ud(k,1),floor_u**2))
           error=2.0d0*MAX( &
                & ABS((vd-vm)/(cd+cm+ABS(vd)+ABS(vm)+floor_u)) , &
                & ABS((vm-vg)/(cm+cg+ABS(vm)+ABS(vg)+floor_u)) )
           ok(k) = ok(k) .or. error > err_grad_u
        end do
     end do
  end if

end subroutine hydro_refine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_llf(qleft,qright,fgdnv,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  ! dummy arguments
  integer::ngrid
  real(dp),dimension(1:nvector,1:nvar)::qleft,qright
  real(dp),dimension(1:nvector,1:nvar+1)::fgdnv

  ! local arrays
  real(dp),dimension(1:nvector,1:nvar+1),save::fleft,fright
  real(dp),dimension(1:nvector,1:nvar+1),save::uleft,uright
  real(dp),dimension(1:nvector),save::cmax

  ! local variables
  integer::i,n
  real(dp)::smallp, entho
  real(dp)::rl   ,ul   ,pl   ,cl
  real(dp)::rr   ,ur   ,pr   ,cr   

  ! Constants
  smallp = smallc**2/gamma
  entho = one/(gamma-one)

  !===========================
  ! Compute maximum wave speed
  !===========================
  do i = 1,ngrid
     ! Left states
     rl = max(qleft (i,1),smallr)
     ul =     qleft (i,2)
     pl = max(qleft (i,3),rl*smallp)
     cl = gamma*pl
#if NENER>0
     do n = 1,nener
        cl = cl + gamma_rad(n)*qleft(i,ndim+2+n)
     end do
#endif
     cl = sqrt(cl/rl)
     ! Right states
     rr = max(qright(i,1),smallr)
     ur =     qright(i,2)
     pr = max(qright(i,3),rr*smallp)
     cr = gamma*pr
#if NENER>0
     do n = 1,nener
        cr = cr + gamma_rad(n)*qright(i,ndim+2+n)
     end do
#endif
     cr = sqrt(cr/rr)
     ! Local max. wave speed
     cmax(i) = max(abs(ul)+cl,abs(ur)+cr)
  end do

  !===============================
  ! Compute conservative variables  
  !===============================
  do i = 1, ngrid 
     ! Mass density
     uleft (i,1) = qleft (i,1)
     uright(i,1) = qright(i,1)
     ! Normal momentum
     uleft (i,2) = qleft (i,1)*qleft (i,2)
     uright(i,2) = qright(i,1)*qright(i,2)
     ! Total energy
     uleft (i,3) = qleft (i,3)*entho + half*qleft (i,1)*qleft (i,2)**2
     uright(i,3) = qright(i,3)*entho + half*qright(i,1)*qright(i,2)**2
#if NDIM>1
     uleft (i,3) = uleft (i,3)       + half*qleft (i,1)*qleft (i,4)**2
     uright(i,3) = uright(i,3)       + half*qright(i,1)*qright(i,4)**2
#endif
#if NDIM>2
     uleft (i,3) = uleft (i,3)       + half*qleft (i,1)*qleft (i,5)**2
     uright(i,3) = uright(i,3)       + half*qright(i,1)*qright(i,5)**2
#endif
#if NENER>0
     do n = 1,nener
        uleft (i,3) = uleft (i,3) + qleft (i,ndim+2+n)/(gamma_rad(n)-one)
        uright(i,3) = uright(i,3) + qright(i,ndim+2+n)/(gamma_rad(n)-one)
     end do
#endif
  end do
  ! Transverse velocities
#if NDIM > 1
  do n = 4, ndim+2
     do i = 1, ngrid
        uleft (i,n) = qleft (i,1)*qleft (i,n)
        uright(i,n) = qright(i,1)*qright(i,n)
     end do
  end do
#endif
  ! Non-thermal energies
#if NENER>0
  do n = 1, nener
     do i = 1, ngrid
        uleft (i,ndim+2+n) = qleft (i,ndim+2+n)/(gamma_rad(n)-one)
        uright(i,ndim+2+n) = qright(i,ndim+2+n)/(gamma_rad(n)-one)
     end do
  end do
#endif
  ! Other passively advected quantities
#if NVAR > 2+NDIM+NENER
  do n = 3+ndim+nener, nvar
     do i = 1, ngrid
        uleft (i,n) = qleft (i,1)*qleft (i,n)
        uright(i,n) = qright(i,1)*qright(i,n)
     end do
  end do
#endif
  ! Thermal energy
  do i = 1, ngrid
     uleft (i,nvar+1) = qleft (i,3)*entho
     uright(i,nvar+1) = qright(i,3)*entho
  end do

  !==============================
  ! Compute left and right fluxes  
  !==============================
  do i = 1, ngrid 
     ! Mass density
     fleft (i,1) = qleft (i,2)*uleft (i,1)
     fright(i,1) = qright(i,2)*uright(i,1)
     ! Normal momentum
     fleft (i,2) = qleft (i,2)*uleft (i,2) + qleft (i,3)
     fright(i,2) = qright(i,2)*uright(i,2) + qright(i,3)
#if NENER>0
     do n = 1,nener
        fleft (i,2) = fleft (i,2) + qleft (i,ndim+2+n)
        fright(i,2) = fright(i,2) + qright(i,ndim+2+n)
     end do
#endif
     ! Total energy
     fleft (i,3) = qleft (i,2)*(uleft (i,3)+qleft (i,3))
     fright(i,3) = qright(i,2)*(uright(i,3)+qright(i,3))
#if NENER>0
     do n = 1,nener
        fleft (i,3) = fleft (i,3) + qleft (i,2)*qleft (i,ndim+2+n)
        fright(i,3) = fright(i,3) + qright(i,2)*qright(i,ndim+2+n)
     end do
#endif
  end do
  ! Other passively advected quantities
  do n = 4, nvar+1
     do i = 1, ngrid
        fleft (i,n) = qleft (i,2)*uleft (i,n)
        fright(i,n) = qright(i,2)*uright(i,n)
     end do
  end do

  !=============================
  ! Compute Lax-Friedrich fluxes
  !=============================
  do n = 1, nvar+1
     do i = 1, ngrid 
        fgdnv(i,n) = half*(fleft(i,n)+fright(i,n)-cmax(i)*(uright(i,n)-uleft(i,n)))
     end do
  end do

end subroutine riemann_llf
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_hll(qleft,qright,fgdnv,ngrid)
  USE amr_parameters
  USE const
  USE hydro_parameters
  ! 1D HLL Riemann solver
  IMPLICIT NONE
  integer::ngrid
  real(dp),dimension(1:nvector,1:nvar)::qleft,qright
  real(dp),dimension(1:nvector,1:nvar+1)::fgdnv

  real(dp),dimension(1:nvector,1:nvar+1),save::fleft,fright
  real(dp),dimension(1:nvector,1:nvar+1),save::uleft,uright
  real(dp),dimension(1:nvector),save::SL,SR
  integer::i,n
  real(dp)::smallp, entho
  real(dp)::rl   ,ul   ,pl   ,cl
  real(dp)::rr   ,ur   ,pr   ,cr   

  ! Constants
  smallp = smallc**2/gamma
  entho = one/(gamma-one)

  !===========================
  ! Compute maximum wave speed
  !===========================
  do i=1,ngrid
     ! Left states
     rl = max(qleft (i,1),smallr)
     ul =     qleft (i,2)
     pl = max(qleft (i,3),rl*smallp)
     cl = gamma*pl
#if NENER>0
     do n = 1,nener
        cl = cl + gamma_rad(n)*qleft(i,ndim+2+n)
     end do
#endif
     cl = sqrt(cl/rl)
     ! Right states
     rr = max(qright(i,1),smallr)
     ur =     qright(i,2)
     pr = max(qright(i,3),rr*smallp)
     cr = gamma*pr
#if NENER>0
     do n = 1,nener
        cr = cr + gamma_rad(n)*qright(i,ndim+2+n)
     end do
#endif
     cr = sqrt(cr/rr)
     ! Left and right max. wave speed
     SL(i)=min(min(ul,ur)-max(cl,cr),zero)
     SR(i)=max(max(ul,ur)+max(cl,cr),zero)
  end do

  !===============================
  ! Compute conservative variables
  !===============================
  do i = 1, ngrid
     ! Mass density
     uleft (i,1) = qleft (i,1)
     uright(i,1) = qright(i,1)
     ! Normal momentum
     uleft (i,2) = qleft (i,1)*qleft (i,2)
     uright(i,2) = qright(i,1)*qright(i,2)
     ! Total energy
     uleft (i,3) = qleft (i,3)*entho + half*qleft (i,1)*qleft (i,2)**2
     uright(i,3) = qright(i,3)*entho + half*qright(i,1)*qright(i,2)**2
#if NDIM>1
     uleft (i,3) = uleft (i,3)       + half*qleft (i,1)*qleft (i,4)**2
     uright(i,3) = uright(i,3)       + half*qright(i,1)*qright(i,4)**2
#endif
#if NDIM>2
     uleft (i,3) = uleft (i,3)       + half*qleft (i,1)*qleft (i,5)**2
     uright(i,3) = uright(i,3)       + half*qright(i,1)*qright(i,5)**2
#endif
#if NENER>0
     do n = 1,nener
        uleft (i,3) = uleft (i,3) + qleft (i,ndim+2+n)/(gamma_rad(n)-one)
        uright(i,3) = uright(i,3) + qright(i,ndim+2+n)/(gamma_rad(n)-one)
     end do
#endif
  end do
  ! Transverse velocities
#if NDIM > 1
  do n = 4, ndim+2
     do i = 1, ngrid
        uleft (i,n) = qleft (i,1)*qleft (i,n)
        uright(i,n) = qright(i,1)*qright(i,n)
     end do
  end do
#endif
  ! Non-thermal energies
#if NENER>0
  do n = 1, nener
     do i = 1, ngrid
        uleft (i,ndim+2+n) = qleft (i,ndim+2+n)/(gamma_rad(n)-one)
        uright(i,ndim+2+n) = qright(i,ndim+2+n)/(gamma_rad(n)-one)
     end do
  end do
#endif
  ! Other passively advected quantities
#if NVAR > 2+NDIM+NENER
  do n = 3+ndim+nener, nvar
     do i = 1, ngrid
        uleft (i,n) = qleft (i,1)*qleft (i,n)
        uright(i,n) = qright(i,1)*qright(i,n)
     end do
  end do
#endif
  ! Thermal energy
  do i = 1, ngrid
     uleft (i,nvar+1) = qleft (i,3)*entho
     uright(i,nvar+1) = qright(i,3)*entho
  end do

  !==============================
  ! Compute left and right fluxes
  !==============================
  do i = 1, ngrid
     ! Mass density
     fleft (i,1) = uleft (i,2)
     fright(i,1) = uright(i,2)
     ! Normal momentum
     fleft (i,2) = qleft (i,3)+uleft (i,2)*qleft (i,2)
     fright(i,2) = qright(i,3)+uright(i,2)*qright(i,2)
#if NENER>0
     do n = 1,nener
        fleft (i,2) = fleft (i,2) + qleft (i,ndim+2+n)
        fright(i,2) = fright(i,2) + qright(i,ndim+2+n)
     end do
#endif
     ! Total energy
     fleft (i,3) = qleft (i,2)*(uleft (i,3)+qleft (i,3))
     fright(i,3) = qright(i,2)*(uright(i,3)+qright(i,3))
#if NENER>0
     do n = 1,nener
        fleft (i,3) = fleft (i,3) + qleft (i,2)*qleft (i,ndim+2+n)
        fright(i,3) = fright(i,3) + qright(i,2)*qright(i,ndim+2+n)
     end do
#endif
  end do
  ! Other advected quantities
  do n = 4, nvar+1
     do i = 1, ngrid
        fleft (i,n) = qleft (i,2)*uleft (i,n)
        fright(i,n) = qright(i,2)*uright(i,n)
     end do
  end do

  !===================
  ! Compute HLL fluxes
  !===================
  do n = 1, nvar+1
     do i = 1, ngrid
        fgdnv(i,n) = (SR(i)*fleft(i,n)-SL(i)*fright(i,n) &
             & + SR(i)*SL(i)*(uright(i,n)-uleft(i,n)))/(SR(i)-SL(i))
     end do
  end do

end subroutine riemann_hll
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_hllc(qleft,qright,fgdnv,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  ! HLLC Riemann solver (Toro)
  integer::ngrid
  real(dp),dimension(1:nvector,1:nvar)::qleft,qright
  real(dp),dimension(1:nvector,1:nvar+1)::fgdnv

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
  REAL(dp),dimension(1:nener)::eradl,eradr,erado
  REAL(dp),dimension(1:nener)::eradstarl,eradstarr
#endif
  INTEGER ::irad, ivar, i

  ! constants
  smallp = smallc**2/gamma
  entho = one/(gamma-one)

  do i=1,ngrid

     ! Left variables
     rl=max(qleft (i,1),smallr)
     Pl=max(qleft (i,3),rl*smallp)
     ul=    qleft (i,2)

     el=Pl*entho
     ecinl = half*rl*ul*ul
#if NDIM>1
     ecinl=ecinl+half*rl*qleft(i,4)**2
#endif
#if NDIM>2
     ecinl=ecinl+half*rl*qleft(i,5)**2
#endif
     etotl = el+ecinl
#if NENER>0
     do irad=1,nener
        eradl(irad)=qleft(i,2+ndim+irad)/(gamma_rad(irad)-one)
        etotl=etotl+eradl(irad)
     end do
#endif
     Ptotl = Pl
#if NENER>0
     do irad=1,nener
        Ptotl=Ptotl+qleft(i,2+ndim+irad)
     end do
#endif

     ! Right variables
     rr=max(qright(i,1),smallr)
     Pr=max(qright(i,3),rr*smallp)
     ur=    qright(i,2)

     er=Pr*entho
     ecinr = half*rr*ur*ur
#if NDIM>1
     ecinr=ecinr+half*rr*qright(i,4)**2
#endif
#if NDIM>2
     ecinr=ecinr+half*rr*qright(i,5)**2
#endif
     etotr = er+ecinr
#if NENER>0
     do irad=1,nener
        eradr(irad)=qright(i,2+ndim+irad)/(gamma_rad(irad)-one)
        etotr=etotr+eradr(irad)
     end do
#endif
     Ptotr = Pr
#if NENER>0
     do irad=1,nener
        Ptotr=Ptotr+qright(i,2+ndim+irad)
     end do
#endif

     ! Find the largest eigenvalues in the normal direction to the interface
     cfastl=gamma*Pl
#if NENER>0
     do irad = 1,nener
        cfastl = cfastl + gamma_rad(irad)*qleft(i,ndim+2+irad)
     end do
#endif
     cfastl=sqrt(max(cfastl/rl,smallc**2))

     cfastr=gamma*Pr
#if NENER>0
     do irad = 1,nener
        cfastr = cfastr + gamma_rad(irad)*qright(i,ndim+2+irad)
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
     fgdnv(i,1) = ro*uo
     fgdnv(i,2) = ro*uo*uo+Ptoto
     fgdnv(i,3) = (etoto+Ptoto)*uo
     ! Transverse velocities
#if NDIM > 1
     do ivar = 4,ndim+2
        if(ustar>0)then
           fgdnv(i,ivar) = ro*uo*qleft (i,ivar)
        else
           fgdnv(i,ivar) = ro*uo*qright(i,ivar)
        endif
     end do
#endif
     ! Non-thermal energies
#if NENER>0
     do irad = 1,nener
        fgdnv(i,ndim+2+irad) = uo*erado(irad)
     end do
#endif
     ! Other passively advected quantities
#if NVAR > 2+NDIM+NENER
     do ivar = 3+ndim+nener,nvar
        if(ustar>0)then
           fgdnv(i,ivar) = ro*uo*qleft (i,ivar)
        else
           fgdnv(i,ivar) = ro*uo*qright(i,ivar)
        endif
     end do
#endif
     ! Thermal energy
     fgdnv(i,nvar+1) = uo*eo

  end do

end subroutine riemann_hllc
!###########################################################
!###########################################################
!###########################################################
!###########################################################
  


