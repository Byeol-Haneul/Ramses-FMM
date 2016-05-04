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
  integer::ind,ivar,igrid,icell,idim
  integer(kind=8),dimension(0:ndim)::hash_key
  real(dp)::average,ekin,erad

  if(ilevel==nlevelmax)return
  if(noct_tot(ilevel)==0)return
  if(noct_tot(ilevel+1)==0)return
  if(verbose)write(*,111)ilevel
 
  ! Set conservative variable to zero in refined cells
  do ioct=head(ilevel),tail(ilevel)
     do ivar=1,nvar
        do ind=1,twotondim
           if(grid(ioct)%refined(ind))then
              grid(ioct)%uold(ind,ivar)=0.0
           endif
        end do
     end do
  end do

  call open_cache(operation_upload,domain_decompos_amr)

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ioct=head(ilevel+1),tail(ilevel+1)

     ! Get cell and grid index
     hash_key(1:ndim)=grid(ioct)%ckey(1:ndim)
     parent_cell=get_parent_cell(hash_key,grid_dict,.true.,.false.)
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

     ! Average internal energy instead of total energy
     if(interpol_var==1 .or. interpol_var==2)then
        average=0.0d0
        do ind=1,twotondim
           ekin=0.0d0
           do idim=1,ndim
              ekin=ekin+0.5d0*grid(ioct)%uold(ind,idim+1)**2/max(grid(ioct)%uold(ind,1),smallr)
           end do
           erad=0.0d0
#if NENER>0
           do irad=1,nener
              erad=erad+grid(ioct)%uold(ind,ndim+2+irad)
           end do
#endif
           average=average+grid(ioct)%uold(ind,ndim+2)-ekin-erad
        end do
        ! Scatter result to cell
        ekin=0.0d0
        do idim=1,ndim
           ekin=ekin+0.5d0*grid(igrid)%uold(icell,idim+1)**2/max(grid(igrid)%uold(icell,1),smallr)
        end do
        erad=0.0d0
#if NENER>0
        do irad=1,nener
           erad=erad+grid(igrid)%uold(icell,ndim+2+irad)
        end do
#endif
        grid(igrid)%uold(icell,ndim+2)=average/dble(twotondim)+ekin+erad
     endif
  end do

  call close_cache(grid_dict)

111 format('   Entering upload_fine for level',i2)

end subroutine upload_fine
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine interpol_hydro(u1,u2)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::nn
  real(dp),dimension(0:twondim  ,1:nvar)::u1
  real(dp),dimension(1:twotondim,1:nvar)::u2
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
  real(dp),dimension(1:8,1:3)::xc
  real(dp),dimension(0:twondim)::a
  real(dp),dimension(1:ndim)::w
  real(dp)::ekin,mom
  real(dp)::erad

  ! Volume fraction of a fine cell realtive to a coarse cell
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
        ekin=0.0d0
        do idim=1,ndim
           ekin=ekin+0.5d0*u1(j,idim+1)**2/max(u1(j,1),smallr)
        end do
        erad=0.0d0
#if NENER>0
        do irad=1,nener
           erad=erad+u1(j,ndim+2+irad)
        end do
#endif
        u1(j,ndim+2)=u1(j,ndim+2)-ekin-erad

        ! and momenta to velocities
        if(interpol_var==2)then
           do idim=1,ndim
              u1(j,idim+1)=u1(j,idim+1)/max(u1(j,1),smallr)
           end do
        end if
     end do
  end if

  ! Loop over interpolation variables
  do ivar=1,nvar

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

     ! choose central limiter for velocities, mon-cen for 
     ! quantities that should not become negative.
     if(interpol_type==4)then
        if (interpol_var .ne. 2)then
           write(*,*)'interpol_type=4 is designed for interpol_var=2'
           call clean_stop
        end if
        if (ivar>1 .and. (ivar <= 1+ndim))then
           call compute_central(a,w)
        else
           call compute_limiter_central(a,w)
        end if
     end if

     ! Interpolate over children cells
     do ind=1,twotondim
        u2(ind,ivar)=a(0)
        do idim=1,ndim
           u2(ind,ivar)=u2(ind,ivar)+w(idim)*xc(ind,idim)
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
              u2(ind,idim+1)=u2(ind,idim+1)*u2(ind,1)
           end do
        end do
        
        ! correct total momentum keeping the slope fixed
        do idim=1,ndim
           mom=0.
           do ind=1,twotondim
              ! total momentum in children
              mom=mom+u2(ind,idim+1)*oneover_twotondim
           end do
           ! error in momentum
           mom=mom-u1(0,idim+1)*u1(0,1)
           ! correct children
           u2(1:twotondim,idim+1)=u2(1:twotondim,idim+1)-mom
        end do
     end if

     ! convert children internal energy into total energy
     do ind=1,twotondim
        ekin=0.0d0
        do idim=1,ndim
           ekin=ekin+0.5d0*u2(ind,idim+1)**2/max(u2(ind,1),smallr)
        end do
        erad=0.0d0
#if NENER>0
        do irad=1,nener
           erad=erad+u2(ind,ndim+2+irad)
        end do
#endif
        u2(ind,ndim+2)=u2(ind,ndim+2)+ekin+erad
     end do
  end if
  
end subroutine interpol_hydro
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_limiter_minmod(a,w)
  use amr_commons
  use hydro_commons
  implicit none
  real(dp),dimension(0:twondim)::a
  real(dp),dimension(1:ndim)::w
  !---------------
  ! MinMod slope
  !---------------
  integer::idim
  real(dp)::diff_left,diff_right,minmod

  do idim=1,ndim
     diff_left=0.5*(a(2*idim)-a(0))
     diff_right=0.5*(a(0)-a(2*idim-1))
     if(diff_left*diff_right<=0.0)then
        minmod=0.0
     else
        minmod=MIN(ABS(diff_left),ABS(diff_right)) &
             &   *diff_left/ABS(diff_left)
     end if
     w(idim)=minmod
  end do

end subroutine compute_limiter_minmod
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_limiter_central(a,w)
  use amr_commons
  use hydro_commons
  implicit none
  real(dp),dimension(0:twondim)::a
  real(dp),dimension(1:ndim)::w
  !---------------------------
  ! Monotonized Central slope
  !---------------------------
  integer::i,j,idim,ind,ix,iy,iz
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp)::xxc
  real(dp),dimension(1:twotondim)::ac
  real(dp)::corner,kernel,diff_corner,diff_kernel
  real(dp)::max_limiter,min_limiter,limiter

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
     w(idim)=0.25D0*(a(2*idim)-a(2*idim-1))
  end do

  ! Compute corner interpolated values
  do ind=1,twotondim
     ac(ind)=a(0)
  end do
  do idim=1,ndim
     do ind=1,twotondim
        xxc = xc(ind,idim)
        corner=ac(ind)+2.D0*w(idim)*xxc
        ac(ind)=corner
     end do
  end do

  ! Compute max of corners
  corner=ac(1)
  do j=2,twotondim
     corner=MAX(corner,ac(j))
  end do

  ! Compute max of gradient kernel
  kernel=a(1)
  do j=2,twondim
     kernel=MAX(kernel,a(j))
  end do

  ! Compute differences
  diff_kernel=a(0)-kernel
  diff_corner=a(0)-corner

  ! Compute max_limiter
  max_limiter=0.0D0
  if(diff_kernel*diff_corner > 0.0D0)then
     max_limiter=MIN(1.0_dp,diff_kernel/diff_corner)
  end if

  ! Compute min of corners
  corner=ac(1)
  do j=2,twotondim
     corner=MIN(corner,ac(j))
  end do

  ! Compute min of gradient kernel
  kernel=a(1)
  do j=2,twondim
     kernel=MIN(kernel,a(j))
  end do

  ! Compute differences
  diff_kernel=a(0)-kernel
  diff_corner=a(0)-corner

  ! Compute min_limiter
  min_limiter=0.0D0
  if(diff_kernel*diff_corner > 0.0D0)then
     min_limiter=MIN(1.0_dp,diff_kernel/diff_corner)
  end if

  ! Compute limiter
  limiter=MIN(min_limiter,max_limiter)

  ! Correct gradient with limiter
  do idim=1,ndim
     w(idim)=w(idim)*limiter
  end do

end subroutine compute_limiter_central
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_central(a,w)
  use amr_commons
  use hydro_commons
  implicit none
  real(dp),dimension(0:twondim)::a
  real(dp),dimension(1:ndim)::w
  !---------------------------
  ! Unlimited Central slope
  !---------------------------
  integer::idim

  ! Second order central slope
  do idim=1,ndim
     w(idim)=0.25D0*(a(2*idim)-a(2*idim-1))
  end do

end subroutine compute_central
