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
  ! Linear advection test: uniform pressure and velocity, with a
  ! density step in the last spatial dimension.
  !================================================================
  integer::i
  real(kind=8)::zz,rho,p0,vadv

  p0 = 1.0d0
  vadv = 1.0d0

  do i=1,nn
     zz = x(i,ndim)/r%box_size(ndim)
     rho = 1.0d0
     if(zz >= 1.0d0/3.0d0 .and. zz < 2.0d0/3.0d0) then
        rho = 2.0d0
     end if
     q(i,1) = rho
     q(i,2) = 0.0d0
     q(i,3) = 0.0d0
     q(i,4) = 0.0d0
     q(i,1+ndim) = vadv
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
