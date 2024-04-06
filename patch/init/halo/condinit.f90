module halo_parameters
  use amr_parameters, only:dp

  ! Galactic merger IC
  real(dp), dimension(3)::halo_center = 0.0D0
  real(dp)::v_200 = 150
  real(dp)::concentration = 10
  real(dp)::baryon_fraction = 0.15
  real(dp)::lambda = 0.04
  real(dp)::halo_nmin = 1d-6
  real(dp)::halo_tmin = 1d6
  real(dp)::halo_eps = 0.01
  
end module halo_parameters

subroutine read_halo_params(g)
  use amr_commons, only: run_t, global_t
  use halo_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(global_t)::g
  !--------------------------------------------------
  ! Local variables  
  !--------------------------------------------------
  logical::nml_ok=.true.
  character(LEN=80)::infile
  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/halo_params/halo_center,v_200,concentration &
       & ,halo_eps,baryon_fraction,halo_nmin,halo_tmin,lambda

  CALL getarg(1,infile)
  open(1,file=infile)
  rewind(1)
  read(1,NML=halo_params,END=106)
  goto 107
106 write(*,*)' You need to set up namelist &HALO_PARAMS in the parameter file'
  stop
107 continue
  close(1)

  if(.not. nml_ok)then
     if(g%myid==1)write(*,*)'Too many errors in the namelist'
     if(g%myid==1)write(*,*)'Aborting...'
     stop
  end if

end subroutine read_halo_params
!===============================
!=== hydro ic for a NFW halo ===
!===============================
subroutine condinit(r,g,x,q,dx,nn)
  use amr_parameters, only: dp, ndim, nvector
  use hydro_parameters, only: nvar, nener
  use halo_parameters
  use amr_commons, only: run_t, global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  integer ::nn                             ! Number of cells
  real(dp)::dx                             ! Cell size
  real(dp),dimension(1:nvector,1:nprim)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x  ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user (aka code) units:
  ! x(i,1:ndim) are in [0,boxlen]**ndim.
  ! Q is the primitve variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4): d.u,d.v,d.w and Q(i,5): P.
  ! If nvar >= 6, remaining variables are treated as passive
  ! scalars or non-thermal energies in the hydro solver.
  ! Q(:,:) are in user (aka code) units.
  !================================================================
  integer::ivar,i
  real(dp)::xx,yy,zz,rc,rr,v,xc,yc,zc,eps,rho,rrmin,rmax,rrmax,tol,c,pmin
  real(dp)::fc,PI,IN,factor,Mr,v200,Hub,M200,dmin,fb,zhalo,r200,rs,hsmall
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,rhos,rhocrit
  real(dp)::j_max
  logical, save:: init_nml=.false.
  
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Read user-defined halo parameters from the namelist
  if (.not. init_nml) then
     call read_halo_params(g)
     init_nml = .true.
  end if

  pi = acos(-1d0)

  hsmall = 0.7d0
  Hub = hsmall*100.0d0/1.0d3 ! [km/s/kpc]
  rhocrit = 3.0*Hub**2/8.0d0/pi ! in the internal 2.3262e5 Msun/kpc^3 as G=1

  dmin = halo_nmin/scale_nH
  pmin = dmin*halo_tmin/scale_T2

  ! Halo parameters from namelist
  eps = halo_eps ! small like 10 pc
  xc = halo_center(1)+r%boxlen/2
  yc = halo_center(2)+r%boxlen/2
  zc = halo_center(3)+r%boxlen/2
  v200 = v_200 ! in [km/s]
  r200 = v200/hsmall ! in [kpc]
  M200 = (r200*hsmall/1.63d-2)**3/hsmall ! in [Msun]
  c = concentration ! concentration
  j_max = lambda*v200*r200  ! lambda -> specific ang. mom.
  fb = baryon_fraction ! baryon fraction
  rs = r200/c
  eps = eps/rs

  M200 = M200/2.3262d5 ! internal units
  tol = 1.0d-6
  rmax = 2.*r200/rs ! in units of rs
  rhos = rhocrit*200./3.0*c*c*c/(log(1d0+c)-c/(1d0+c))
  
  do i=1,nn ! we are now workin in [kpc], see units
     xx = x(i,1)-xc
     yy = x(i,2)-yc
     zz = x(i,3)-zc
     rc = sqrt(xx**2+yy**2) ! kpc
     rr = sqrt(xx**2+yy**2+zz**2)/rs ! units of rs
     
     rr = max(rr,eps)
     rc = max(rc,eps*rs)

     q(i,1) = fb*rhos/rr/(1.+rr)**2
     q(i,1) = max(q(i,1),dmin)

     Mr = 4.*pi*rhos*rs**3*(log(1+rr)-(rr/(1.+rr)))  ! M(<r)
     v = j_max*Mr/M200/(rr*rs)

     q(i,2) = -v*yy/rc
     q(i,3) = +v*xx/rc

     rrmin = log(rr)
     rrmax = log(rmax)

     q(i,5) = romberg(rrmin,rrmax,tol)
     q(i,5) = max(q(i,5),pmin)

  enddo

  q(1:nn,4) = 0.0

  if(r%metal) then
     q(1:nn,r%imetal) = r%z_ave*0.02d0
  endif

  if(r%entropy)then
     q(1:nn,r%ientropy) = q(1:nn,5)/q(1:nn,1)**r%gamma
  endif

contains
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function fffy(rint)
    implicit none
    ! Computes the integrand
    real(dp)::fffy
    real(dp),intent(in)::rint
    real(dp)::rrr,Mrr,rhorr
 
    rrr = exp(rint)
    rrr = max(rrr,eps)

    rhorr = fb*rhos/rrr/(1.+rrr)**2
    Mrr = 4.*pi*rhos*rs**3*(log(1+rrr)-(rrr/(1.+rrr))) ! M(<r)

    rrr = rrr*rs ! physical kpc
    fffy = rhorr*Mrr/rrr

    return
  end function fffy
  !cccccccccccccccccccccccccccccccccccccccccccccccccccc
  function romberg(a,b,tol)
    implicit none
    real(dp)::romberg
    !
    !     Romberg returns the integral from a to b of f(x)dx using Romberg 
    !     integration. The method converges provided that f(x) is continuous 
    !     in (a,b). The function f must be double precision and must be 
    !     declared external in the calling routine.  
    !     tol indicates the desired relative accuracy in the integral.
    !
    integer::maxiter=16,maxj=5
    real(dp),dimension(100)::g
    real(dp)::a,b,tol,fourj
    real(dp)::h,error,gmax,g0,g1
    integer::nint,i,j,k,jmax

    h=0.5d0*(b-a)
    gmax=h*(fffy(a)+fffy(b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if(.not.  (i>maxiter.or.(i>5.and.abs(error)<tol)))then
       ! Calculate next trapezoidal rule approximation to integral.   
       g0=0.0d0
       do k=1,nint
          g0=g0+fffy(a+(k+k-1)*h)
       end do
       g0=0.5d0*g(1)+h*g0
       h=0.5d0*h
       nint=nint+nint
       jmax=min(i,maxj)
       fourj=1.0d0
       
       do j=1,jmax
          ! Use Richardson extrapolation.
          fourj=4.0d0*fourj
          g1=g0+(g0-g(j))/(fourj-1.0d0)
          g(j)=g0
          g0=g1
       enddo
       if (abs(g0).gt.tol) then
          error=1.0d0-gmax/g0
       else
          error=gmax
       end if
       gmax=g0
       g(jmax+1)=g0
       go to 10
    end if
    romberg=g0

    if (i>maxiter.and.abs(error)>tol) &
         &    write(*,*) 'Romberg failed to converge; integral, error=', &
         &    romberg,error

    return
  end function romberg

end subroutine condinit

