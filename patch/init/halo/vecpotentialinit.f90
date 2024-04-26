!================================================================
!================================================================
!================================================================
!================================================================
subroutine vecpotentialinit(r,g,x,A,idim,nn)
  use amr_parameters, only: dp, ndim, nvector
  use amr_commons, only: run_t, global_t
  use halo_parameters
  use const
  implicit none
  type(run_t)::r
  type(global_t)::g
  integer::nn                             ! Number of cells
  integer::idim                           ! Direction of the component
  real(dp),dimension(1:nvector)::A        ! Vector potential component
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user (aka code) units:
  ! x(i,1:ndim) are in [0,boxlen]**ndim.
  ! A is the component of the vector potential corresponding
  ! to direction idim.
  ! A(:) is in user (aka code) units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,rc,rr,xc,yc,zc,eps,c,dnfw
  real(dp)::v200,r200,rs,hsmall
  
  ! Halo parameters from namelist
  hsmall = 0.7d0
  eps = halo_eps ! small like 10 pc
  xc = halo_center(1)+r%boxlen/2
  yc = halo_center(2)+r%boxlen/2
  zc = halo_center(3)+r%boxlen/2
  v200 = v_200 ! in [km/s]
  r200 = v200/hsmall ! in [kpc]
  c = concentration ! concentration
  rs = r200/c
  eps = eps/rs
  
  do i = 1,nn ! we are now workin in [kpc], see units

     xx = x(i,1)-xc
     yy = x(i,2)-yc
     zz = x(i,3)-zc
     rc = sqrt(xx**2+yy**2) ! kpc
     rr = sqrt(xx**2+yy**2+zz**2)/rs ! units of rs
     
     rr = max(rr,eps)
     rc = max(rc,eps*rs)

     dnfw = 1./rr/(1.+rr)**2

     select case (mag_type)
        
     case (0) ! constant done elsewhere
        A(i)=0

     case (1) ! toroidal
        if(idim==1)A(i)=0
        if(idim==2)A(i)=0        
        if(idim==3)A(i)=B_s*rc*dnfw**two3rd

     case (2) ! dipole
        if(idim==1)A(i)=-B_s*yy*dnfw**two3rd
        if(idim==2)A(i)=+B_s*xx*dnfw**two3rd
        if(idim==3)A(i)=0

     case (3) ! quadrupole
        if(idim==1)A(i)=-B_s*zz*dnfw**two3rd*yy/rc
        if(idim==2)A(i)=+B_s*zz*dnfw**two3rd*xx/rc
        if(idim==3)A(i)=0
        
     end select

  enddo
  
end subroutine vecpotentialinit

