!###########################################################
!########################################################### 
!###########################################################
!###########################################################
subroutine upload_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the hydro variables.
  !----------------------------------------------------------------------
  integer::ioct,parent_cell,get_parent_cell
  integer::ind,ivar,igrid,icell
  integer(kind=8),dimension(0:ndim)::hash_key
  real(dp)::average

  if(ilevel==nlevelmax)return
  if(noct_tot(ilevel)==0)return
  if(noct_tot(ilevel+1)==0)return
  if(verbose)write(*,111)ilevel
 
  ! Set conservative variable to zero in refined cells
  do ioct=head(ilevel),tail(ilevel)
     do ivar=1,nvar
        do ind=1,twotondim
           if(grid(ioct)%refined(ind))grid(ioct)%uold(ind,ivar)=0.0
        end do
     end do
  end do

  cache_operation=operation_upload
  call open_cache

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ioct=head(ilevel+1),tail(ilevel+1)
     hash_key(1:ndim)=grid(ioct)%ckey(1:ndim)
     parent_cell=get_parent_cell(hash_key,.true.,.false.)
     igrid=(parent_cell-1)/twotondim+1
     icell=parent_cell-(igrid-1)*twotondim
     ! Average conservative variables
     do ivar=1,nvar
        average=0.0d0
        do ind=1,twotondim
           average=average+grid(ioct)%uold(ind,ivar)
        end do
        ! Scatter result to cell
        grid(igrid)%uold(icell,ivar)=average/dble(twotondim)
     end do
  end do

  call close_cache

111 format('   Entering upload_fine for level',i2)

end subroutine upload_fine
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
#ifdef TOTO
subroutine interpol_hydro(u1,u2,nn)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar)::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar)::u2
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  ! The interpolated variables are:
  ! interpol_var=0: rho, rho u and E
  ! interpol_var=1: rho, rho u and rho epsilon
  ! interpol_var=2: rho, u and rho epsilon
  ! The interpolation method is:
  ! interpol_type=0 straight injection
  ! interpol_type=1 linear interpolation with MinMod slope
  ! interpol_type=2 linear interpolation with Monotonized Central slope
  ! interpol_type=3 linear interpolation with unlimited Central slope
  ! interpol_type=4 in combination with interpol_var==2
  !                 type 3 for velocity and type 2 for density and
  !                 internal energy.
  !----------------------------------------------------------
  integer::i,j,ivar,irad,idim,ind,ix,iy,iz,ind2
  real(dp)::oneover_twotondim
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,0:twondim),save::a
  real(dp),dimension(1:nvector,1:ndim),save::w
  real(dp),dimension(1:nvector),save::ekin,mom
  real(dp),dimension(1:nvector),save::erad

  ! volume fraction of a fine cell realtive to a coarse cell
  oneover_twotondim=1.D0/dble(twotondim)

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! If necessary, convert father total energy into internal energy
  if(interpol_var==1 .or. interpol_var==2)then
     do j=0,twondim
        ekin(1:nn)=0.0d0
        do idim=1,ndim
           do i=1,nn
              ekin(i)=ekin(i)+0.5d0*u1(i,j,idim+1)**2/max(u1(i,j,1),smallr)
           end do
        end do
        erad(1:nn)=0.0d0
#if NENER>0
        do irad=1,nener
           do i=1,nn
              erad(i)=erad(i)+u1(i,j,ndim+2+irad)
           end do
        end do
#endif
        do i=1,nn
           u1(i,j,ndim+2)=u1(i,j,ndim+2)-ekin(i)-erad(i)
        end do

        ! and momenta to velocities
        if(interpol_var==2)then
           do idim=1,ndim
              do i=1,nn
                 u1(i,j,idim+1)=u1(i,j,idim+1)/max(u1(i,j,1),smallr)
              end do
           end do
        end if
     end do
  end if

  ! Loop over interpolation variables
  do ivar=1,nvar

     ! Load father variable
     do j=0,twondim
        do i=1,nn 
           a(i,j)=u1(i,j,ivar)
        end do
     end do

     ! Reset gradient
     w(1:nn,1:ndim)=0.0D0

     ! Compute gradient with chosen limiter
     if(interpol_type==1)call compute_limiter_minmod(a,w,nn)
     if(interpol_type==2)call compute_limiter_central(a,w,nn)
     if(interpol_type==3)call compute_central(a,w,nn)
     ! choose central limiter for velocities, mon-cen for 
     ! quantities that should not become negative.
     if(interpol_type==4)then
        if (interpol_var .ne. 2)then
           write(*,*)'interpol_type=4 is designed for interpol_var=2'
           call clean_stop
        end if
        if (ivar>1 .and. (ivar <= 1+ndim))then
           call compute_central(a,w,nn)
        else
           call compute_limiter_central(a,w,nn)
        end if
     end if

     ! Interpolate over children cells
     do ind=1,twotondim
        u2(1:nn,ind,ivar)=a(1:nn,0)
        do idim=1,ndim
           do i=1,nn
              u2(i,ind,ivar)=u2(i,ind,ivar)+w(i,idim)*xc(ind,idim)
           end do
        end do
     end do

  end do
  ! End loop over variables
  
  ! If necessary, convert children internal energy into total energy
  ! and velocities back to momenta
  if(interpol_var==1 .or. interpol_var==2)then
     if(interpol_var==2)then
        do ind=1,twotondim        
           do idim=1,ndim
              do i=1,nn
                 u2(i,ind,idim+1)=u2(i,ind,idim+1)*u2(i,ind,1)
              end do
           end do
        end do
        
        !correct total momentum keeping the slope fixed
        do idim=1,ndim
           mom(1:nn)=0.
           do ind=1,twotondim
              do i=1,nn
                 ! total momentum in children
                 mom(i)=mom(i)+u2(i,ind,idim+1)*oneover_twotondim
              end do
           end do
           do i=1,nn
              ! error in momentum
              mom(i)=mom(i)-u1(i,0,idim+1)*u1(i,0,1)
              ! correct children
              u2(i,1:twotondim,idim+1)=u2(i,1:twotondim,idim+1)-mom(i)              
           end do
        end do
     end if

     ! convert children internal energy into total energy
     do ind=1,twotondim
        ekin(1:nn)=0.0d0
        do idim=1,ndim
           do i=1,nn
              ekin(i)=ekin(i)+0.5d0*u2(i,ind,idim+1)**2/max(u2(i,ind,1),smallr)
           end do
        end do
        erad(1:nn)=0.0d0
#if NENER>0
        do irad=1,nener
           do i=1,nn
              erad(i)=erad(i)+u2(i,ind,ndim+2+irad)
           end do
        end do
#endif
        do i=1,nn
           u2(i,ind,ndim+2)=u2(i,ind,ndim+2)+ekin(i)+erad(i)
        end do
     end do
  end if
  
end subroutine interpol_hydro
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_limiter_minmod(a,w,nn)
  use amr_commons
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  !---------------
  ! MinMod slope
  !---------------
  integer::i,idim
  real(dp)::diff_left,diff_right,minmod

  do idim=1,ndim
     do i=1,nn
        diff_left=0.5*(a(i,2*idim)-a(i,0))
        diff_right=0.5*(a(i,0)-a(i,2*idim-1))
        if(diff_left*diff_right<=0.0)then
           minmod=0.0
        else
           minmod=MIN(ABS(diff_left),ABS(diff_right)) &
                &   *diff_left/ABS(diff_left)
        end if
        w(i,idim)=minmod
     end do
  end do

end subroutine compute_limiter_minmod
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_limiter_central(a,w,nn)
  use amr_commons
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  !---------------------------
  ! Monotonized Central slope
  !---------------------------
  integer::i,j,idim,ind,ix,iy,iz
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp)::xxc
  real(dp),dimension(1:nvector,1:twotondim),save::ac
  real(dp),dimension(1:nvector),save::corner,kernel,diff_corner,diff_kernel
  real(dp),dimension(1:nvector),save::max_limiter,min_limiter,limiter

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! Second order central slope
  do idim=1,ndim
     do i=1,nn
        w(i,idim)=0.25D0*(a(i,2*idim)-a(i,2*idim-1))
     end do
  end do

  ! Compute corner interpolated values
  do ind=1,twotondim
     do i=1,nn
        ac(i,ind)=a(i,0)
     end do
  end do
  do idim=1,ndim
     do ind=1,twotondim
        xxc = xc(ind,idim)
        do i=1,nn
           corner(i)=ac(i,ind)+2.D0*w(i,idim)*xxc
        end do
        do i=1,nn
           ac(i,ind)=corner(i)
        end do
     end do
  end do

  ! Compute max of corners
  do i=1,nn
     corner(i)=ac(i,1)
  end do
  do j=2,twotondim
     do i=1,nn
        corner(i)=MAX(corner(i),ac(i,j))
     end do
  end do

  ! Compute max of gradient kernel
  do i=1,nn
     kernel(i)=a(i,1)
  end do
  do j=2,twondim
     do i=1,nn
        kernel(i)=MAX(kernel(i),a(i,j))
     end do
  end do

  ! Compute differences
  do i=1,nn
     diff_kernel(i)=a(i,0)-kernel(i)
     diff_corner(i)=a(i,0)-corner(i)
  end do

  ! Compute max_limiter
  max_limiter=0.0D0
  do i=1,nn
     if(diff_kernel(i)*diff_corner(i) > 0.0D0)then
        max_limiter(i)=MIN(1.0_dp,diff_kernel(i)/diff_corner(i))
     end if
  end do

  ! Compute min of corners
  do i=1,nn
     corner(i)=ac(i,1)
  end do
  do j=2,twotondim
     do i=1,nn
        corner(i)=MIN(corner(i),ac(i,j))
     end do
  end do

  ! Compute min of gradient kernel
  do i=1,nn
     kernel(i)=a(i,1)
  end do
  do j=2,twondim
     do i=1,nn
        kernel(i)=MIN(kernel(i),a(i,j))
     end do
  end do

  ! Compute differences
  do i=1,nn
     diff_kernel(i)=a(i,0)-kernel(i)
     diff_corner(i)=a(i,0)-corner(i)
  end do

  ! Compute max_limiter
  min_limiter=0.0D0
  do i=1,nn
     if(diff_kernel(i)*diff_corner(i) > 0.0D0)then
        min_limiter(i)=MIN(1.0_dp,diff_kernel(i)/diff_corner(i))
     end if
  end do

  ! Compute limiter
  do i=1,nn
     limiter(i)=MIN(min_limiter(i),max_limiter(i))
  end do

  ! Correct gradient with limiter
  do idim=1,ndim
     do i=1,nn
        w(i,idim)=w(i,idim)*limiter(i)
     end do
  end do

end subroutine compute_limiter_central
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_central(a,w,nn)
  use amr_commons
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  !---------------------------
  ! Unlimited Central slope
  !---------------------------
  integer::i,idim

  ! Second order central slope
  do idim=1,ndim
     do i=1,nn
        w(i,idim)=0.25D0*(a(i,2*idim)-a(i,2*idim-1))
     end do
  end do

end subroutine compute_central
#endif
