!================================================================
!================================================================
!================================================================
!================================================================
subroutine vecpotentialinit(r,g,x,A,idim,nn)
  use amr_parameters, only: dp, ndim, nvector
  use amr_commons, only: run_t, global_t
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
#define LOOP 4
#define OT 5
#define PONO 6
  integer::i
#if INIT==LOOP
  real(dp)::R0, A0, xx, yy
#endif
#if INIT==OT
  real(dp)::B0, pi, xx, yy
#endif
#if INIT==PONO
  real(dp)::A0, twopi, xx, yy, zz, rr, tt
#else
  do i = 1,nn
     A(i)=0.0
  end do
#endif
  
#if INIT==LOOP
  R0 = 0.3
  A0 = 1d-3
  do i = 1,nn
     xx = x(i,1)-r%boxlen/2.0
     yy = x(i,2)-r%boxlen/2.0
     if(idim==1)A(i) = 0.0
     if(idim==2)A(i) = 0.0
     if(idim==3)A(i) = A0*max(R0-sqrt(xx**2+yy**2),0.0_dp)
  end do
#endif

#if INIT==OT
  pi = ACOS(-1.0d0)
  B0 = 1.0/sqrt(4.0*pi)
  do i = 1,nn
     xx = x(i,1)
     yy = x(i,2)
     A(i) = B0*(cos(4.0*pi*xx)/(4.0*pi)+cos(2.0*pi*yy)/(2.0*pi))
  end do
#endif

#if INIT==PONO
  A0 = 0.1
  twopi = 2.0d0*ACOS(-1.0d0)
  do i = 1,nn
     xx = x(i,1)-r%boxlen/2.0
     yy = x(i,2)-r%boxlen/2.0
     zz = x(i,3)
     rr = sqrt(xx**2+yy**2)
     if(rr > 0.0.and.rr < 2.0)then
        if(yy>0)then
           tt=acos(xx/rr)
        else
           tt=-acos(xx/rr)+twopi
        endif
        A(i) = A0*(2d0-rr)*rr*cos(tt)*sin(twopi*zz/r%boxlen)
        if(idim==1)A(i) = A(i)*cos(tt)
        if(idim==2)A(i) = A(i)*sin(tt)
        if(idim==3)A(i) = 0.0
     else
        A(i)=0.0
     endif
  end do
#endif

end subroutine vecpotentialinit
