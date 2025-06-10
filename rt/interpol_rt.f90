!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine interpol_rt(u1,u2,interpol_type)
  use amr_parameters, only: ndim, twotondim, twondim
  use rt_parameters, only: nrtvar
  implicit none
  integer::interpol_type
  real(kind=8),dimension(0:twondim  ,1:nrtvar)::u1
  real(kind=8),dimension(1:twotondim,1:nrtvar)::u2
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  ! The interpolation method is:
  ! interpol_type=0 straight injection
  ! interpol_type=1 linear interpolation with MinMod slope
  ! interpol_type=2 linear interpolation with Monotonized Central slope
  ! interpol_type=3 linear interpolation with unlimited Central slope
  !----------------------------------------------------------
  integer::j,ivar,idim,ind,ix,iy,iz
  real(kind=8),dimension(1:8,1:3)::xc
  real(kind=8),dimension(0:twondim)::a
  real(kind=8),dimension(1:ndim)::w

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! Loop over interpolation variables
  do ivar=1,nrtvar

     ! Load father variable
     do j=0,twondim
        a(j)=u1(j,ivar)
     end do

     ! Reset gradient
     w(1:ndim)=0.0D0

     ! Compute gradient with chosen limiter
     if(interpol_type==1)call compute_limiter_minmod(a,w)
     if(interpol_type==2)call compute_limiter_central(a,w)
     if(interpol_type==3)call compute_central(a,w)

     ! Interpolate over children cells
     do ind=1,twotondim
        u2(ind,ivar)=a(0)
        do idim=1,ndim
           u2(ind,ivar)=u2(ind,ivar)+w(idim)*xc(ind,idim)
        end do
     end do

  end do
  ! End loop over variables

end subroutine interpol_rt
