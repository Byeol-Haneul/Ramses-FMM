#ifdef GRAV
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Mask restriction (bottom-up)
! ------------------------------------------------------------------------

recursive subroutine r_restrict_mask(pst,input_array,input_size,output_array,output_size)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer,dimension(1:output_size)::next_output_array
  integer::ilevel
  logical::allmasked
  
  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_RESTRICT_MASK,pst%iUpper+1,input_size,output_size,input_array)
     call r_restrict_mask(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size,next_output_array)
     output_array(1)=output_array(1)*next_output_array(1)
  else
     ilevel=input_array(1)
     call restrict_mask(pst%s,ilevel,allmasked)
     if(allmasked)then
        output_array(1)=1
     else
        output_array(1)=0
     endif
  endif

end subroutine r_restrict_mask

subroutine restrict_mask(s,ifinelevel,allmasked)
  use amr_parameters, only: dp,nvector,nhilbert,ndim,twotondim
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils_p
  use cache_commons
  use hilbert
  use hash
  implicit none
  type(ramses_t)::s
  integer, intent(in) :: ifinelevel
  logical, intent(inout) :: allmasked

  integer(kind=8),dimension(0:ndim) :: hash_key
  integer :: ichild,ind,igrid,icell
  real(dp) :: ngpmask, mask_max
  real(dp) :: dtwotondim = (twotondim)
  type(oct),pointer::gridp

  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Initialize volume fraction to zero at coarse level
  do igrid=m%head_mg(ifinelevel-1),m%tail_mg(ifinelevel-1)
     do ind=1,twotondim
        m%grid(igrid)%f(ind,3)=0d0
     end do
  end do
  
  hash_key(0)=ifinelevel
  
  call open_cache(s,operation_restrict_mask,domain_decompos_mg)
  
  ! Loop over grids
  do ichild=m%head_mg(ifinelevel),m%tail_mg(ifinelevel)

     ! Loop over cells
     do ind=1,twotondim        

        hash_key(1:ndim)=m%grid(ichild)%ckey(1:ndim)

        ! Get parent cell using write-only cache
        call get_parent_cell_p(s,hash_key,m%mg_dict,gridp,icell,.true.,.false.)

        ! Convert mask value to volume fraction
        ngpmask=(1d0+m%grid(ichild)%f(ind,3))/2d0/dtwotondim
        gridp%f(icell,3)=gridp%f(icell,3)+ngpmask

     end do
  end do

  call close_cache(s,m%mg_dict)

  ! Convert volume fraction back to to mask value for coarse level
  do igrid=m%head_mg(ifinelevel-1),m%tail_mg(ifinelevel-1)
     do ind=1,twotondim
        m%grid(igrid)%f(ind,3)=2d0*m%grid(igrid)%f(ind,3)-1d0
     end do
  end do
  
  ! Check mask state at coarse level
  mask_max=-1.0
  do igrid=m%head_mg(ifinelevel-1),m%tail_mg(ifinelevel-1)
     do ind=1,twotondim
        mask_max=MAX(mask_max,m%grid(igrid)%f(ind,3))
     end do
  end do
  allmasked=(mask_max<=0d0)

  end associate

end subroutine restrict_mask

subroutine init_flush_restrict_mask(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  
#ifdef GRAV
  do ind=1,twotondim
     grid%f(ind,3)=0.0d0
  end do
#endif
  
end subroutine init_flush_restrict_mask

subroutine pack_flush_restrict_mask(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp(ind)=grid%f(ind,3)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_restrict_mask

subroutine unpack_flush_restrict_mask(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

  msg=transfer(msg_array,msg)
  
#ifdef GRAV
  do ind=1,twotondim
     grid%f(ind,3)=grid%f(ind,3)+msg%realdp(ind)
  end do
#endif

end subroutine unpack_flush_restrict_mask

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Residual computation
! ------------------------------------------------------------------------

recursive subroutine r_cmp_residual_mg(pst,input_array,input_size,output_array,output_size)
  use ramses_commons, only: pst_t
  use mdl_parameters
  use hash
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel,ifine

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_CMP_RESIDUAL_MG,pst%iUpper+1,input_size,output_size,input_array)
     call r_cmp_residual_mg(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     ilevel=input_array(1)
     ifine=input_array(2)
     if(ifine==ilevel)then
        call cmp_residual_mg(pst%s,pst%s%m%grid_dict,ifine)
     else
        call cmp_residual_mg(pst%s,pst%s%m%mg_dict,ifine)
     endif
  endif

end subroutine r_cmp_residual_mg

subroutine cmp_residual_mg(s,hash_dict, ilevel)
  use amr_parameters, only: dp,nvector,nhilbert,ndim,twondim,twotondim
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils_p
  use cache_commons
  use hilbert
  use hash
  implicit none
  type(ramses_t)::s
  integer, intent(in) :: ilevel
  type(hash_table) :: hash_dict

  ! Computes the residual for MG levels, and stores it into grid(igrid)%f(ind,1)
    
  integer, dimension(1:3,1:2,1:8) :: iii, jjj
  real(dp),dimension(1:twotondim,0:twondim)::phi_nbor,dis_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer(kind=8),dimension(0:ndim) :: hash_nbor
  real(dp) :: dx, oneoverdx2, phi_c, dis_c, nb_sum
  integer  :: igrid, ind, inbor, idim, igridn, id, ig
  real(dp) :: dtwondim = (twondim)
  type(oct),pointer::gridp

  associate(r=>s%r,g=>s%g,m=>s%m)
  
  ! Set constants
  dx = r%boxlen/2**ilevel
  oneoverdx2 = 1.0d0/(dx*dx)
  
  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  
  call open_cache(s,operation_mg,domain_decompos_mg)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=m%head_mg(ilevel),m%tail_mg(ilevel)

     ! Get central oct potential
     do ind=1,twotondim
        phi_nbor(ind,0)=m%grid(igrid)%phi(ind)
        dis_nbor(ind,0)=m%grid(igrid)%f(ind,3)
     end do

     ! Get neighboring octs potential
     do inbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=m%grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=m%ckey_max(ilevel)-1
           if(hash_nbor(idim)==m%ckey_max(ilevel))hash_nbor(idim)=0
        enddo

        ! Get neighbouring grid using read-only cache
        call get_grid_p(s,hash_nbor,hash_dict,gridp,.false.,.true.)

        ! If grid exists, then copy into array
        if(associated(gridp))then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=gridp%phi(ind)
              dis_nbor(ind,inbor)=gridp%f(ind,3)
           end do

        ! Otherwise set to zero and outside
        else
           do ind=1,twotondim
              phi_nbor(ind,inbor)=0.0
              dis_nbor(ind,inbor)=-1.0
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim

        ! Compute residual using 6 neighbors potential
        phi_c=m%grid(igrid)%phi(ind)
        dis_c=m%grid(igrid)%f(ind,3)

        nb_sum=0.0

        if(.not. btest(m%grid(igrid)%flag2(ind),0))then ! No scan needed

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 nb_sum=nb_sum+phi_nbor(id,ig)
              end do
           end do

        else ! Scan is required

           ! If cell is outside, set residual to zero
           if(dis_c<=0.0)then
              m%grid(igrid)%f(ind,1)=0.0
              cycle
           else

              ! Loop over neighbours
              do inbor=1,2
                 do idim=1,ndim
                    id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                    if(dis_nbor(id,ig)<=0.0)then
                       nb_sum=nb_sum+phi_c/dis_c*dis_nbor(id,ig)
                    else
                       nb_sum=nb_sum+phi_nbor(id,ig)
                    endif
                 end do
              end do

           endif

        endif

        ! Store residual in f(ind,1)
        m%grid(igrid)%f(ind,1)=-oneoverdx2*( nb_sum - dtwondim*phi_c )+m%grid(igrid)%f(ind,2)

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache(s,hash_dict)

  end associate
  
end subroutine cmp_residual_mg

subroutine pack_fetch_mg(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_twin_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_twin_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp_phi(ind)=grid%phi(ind)
     msg%realdp_dis(ind)=grid%f(ind,3)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_mg

subroutine unpack_fetch_mg(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_twin_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_twin_realdp)::msg

  msg=transfer(msg_array,msg)
  
#ifdef GRAV
  do ind=1,twotondim
     grid%phi(ind)=msg%realdp_phi(ind)
     grid%f(ind,3)=msg%realdp_dis(ind)
  end do
#endif

end subroutine unpack_fetch_mg

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Gauss-Seidel Red-Black sweeps
! ------------------------------------------------------------------------

recursive subroutine r_gauss_seidel_mg(pst,input_array,input_size,output_array,output_size)
  use ramses_commons, only: pst_t
  use mdl_parameters
  use hash
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel,ifine,isafe,iredstep
  logical::safe,redstep

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_GAUSS_SEIDEL_MG,pst%iUpper+1,input_size,output_size,input_array)
     call r_gauss_seidel_mg(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     ilevel=input_array(1)
     ifine=input_array(2)
     isafe=input_array(3)
     iredstep=input_array(4)
     safe=(isafe==1)
     redstep=(iredstep==1)
     if(ifine==ilevel)then
        call gauss_seidel_mg(pst%s,pst%s%m%grid_dict,ifine,safe,redstep)
     else
        call gauss_seidel_mg(pst%s,pst%s%m%mg_dict,ifine,safe,redstep)
     endif
  endif

end subroutine r_gauss_seidel_mg

subroutine gauss_seidel_mg(s,hash_dict,ilevel,safe,redstep)
  use amr_parameters, only: dp,nvector,nhilbert,ndim,twondim,twotondim
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils_p
  use cache_commons
  use hilbert
  use hash
  implicit none
  type(ramses_t)::s
  integer, intent(in) :: ilevel
  logical, intent(in) :: safe
  logical, intent(in) :: redstep
  type(hash_table) :: hash_dict

  ! Perform a Gauss-Seidel update of grid(igrid)%phi(ind).
  ! The domain mask is also needed.
  
  integer, dimension(1:3,1:2,1:8) :: iii, jjj
  real(dp),dimension(1:twotondim,0:twondim)::phi_nbor,dis_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer(kind=8),dimension(0:ndim) :: hash_nbor
  real(dp) :: phi_c, dis_c, dx2, nb_sum, weight
  integer  :: igrid, ind, inbor, idim, igridn, id, ig, ind0
  real(dp) :: dtwondim = (twondim)
  type(oct),pointer::gridp

  integer, dimension(1:4) :: ired, iblack
  
  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Set constants
  dx2  = (r%boxlen/2**ilevel)**2
  
  ired  (1:4)=(/1,4,6,7/)
  iblack(1:4)=(/2,3,5,8/)
  
  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  
  call open_cache(s,operation_mg,domain_decompos_mg)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=m%head_mg(ilevel),m%tail_mg(ilevel)
     
     ! Loop over cells
     do ind=1,twotondim

        ! Get central oct potential and distance
        phi_nbor(ind,0)=m%grid(igrid)%phi(ind)
        dis_nbor(ind,0)=m%grid(igrid)%f(ind,3)

     end do

     ! Get neighboring octs potential
     do inbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=m%grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=m%ckey_max(ilevel)-1
           if(hash_nbor(idim)==m%ckey_max(ilevel))hash_nbor(idim)=0
        enddo

        ! Get neighbouring grid using a read-only cache
        call get_grid_p(s,hash_nbor,hash_dict,gridp,.false.,.true.)

        ! If grid exists, then copy into array
        if(associated(gridp))then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=gridp%phi(ind)
              dis_nbor(ind,inbor)=gridp%f(ind,3)
           end do

        ! Otherwise set to zero and outside
        else
           
           do ind=1,twotondim
              phi_nbor(ind,inbor)=0.0
              dis_nbor(ind,inbor)=-1.0
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells, with red/black ordering
     do ind0=1,twotondim/2      ! Only half of the cells for a red or black sweep

        if(redstep) then
           ind = ired  (ind0)
        else
           ind = iblack(ind0)
        end if

        ! Compute residual using 6 neighbors potential
        phi_c=m%grid(igrid)%phi(ind)
        dis_c=m%grid(igrid)%f(ind,3)

        nb_sum=0.0

        if(.not. btest(m%grid(igrid)%flag2(ind),0))then ! No scan needed

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 nb_sum=nb_sum+phi_nbor(id,ig)
              end do
           end do

           ! Update the potential, solving for potential on icell_amr
           m%grid(igrid)%phi(ind)=(nb_sum-dx2*m%grid(igrid)%f(ind,2))/dtwondim

        else ! Scan is required

           ! If cell is outside, don't update phi
           if(dis_c<=0.0)cycle
           if(safe .and. dis_c<1.0)cycle

           weight=0.0d0   ! Central weight for "Solve G-S"

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 if(dis_nbor(id,ig)<=0.0)then
                    weight=weight+dis_nbor(id,ig)/dis_c
                 else
                    nb_sum=nb_sum+phi_nbor(id,ig)
                 endif
              end do
           end do

           ! Update the potential
           m%grid(igrid)%phi(ind)=(nb_sum-dx2*m%grid(igrid)%f(ind,2))/(dtwondim - weight)

        endif

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache(s,hash_dict)

  end associate

end subroutine gauss_seidel_mg

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Reset correction
! ------------------------------------------------------------------------

recursive subroutine r_reset_correction(pst,input_array,input_size,output_array,output_size)
  use amr_parameters, only: twotondim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::igrid,ilevel

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_RESET_CORRECTION,pst%iUpper+1,input_size,output_size,input_array)
     call r_reset_correction(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     ilevel=input_array(1)
     do igrid=pst%s%m%head_mg(ilevel),pst%s%m%tail_mg(ilevel)
        pst%s%m%grid(igrid)%phi(1:twotondim)=0.0d0
     end do
  endif

end subroutine r_reset_correction

! ------------------------------------------------------------------------
! Residual restriction
! ------------------------------------------------------------------------

recursive subroutine r_restrict_residual(pst,input_array,input_size,output_array,output_size)
  use amr_parameters, only: twotondim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_RESTRICT_RESIDUAL,pst%iUpper+1,input_size,output_size,input_array)
     call r_restrict_residual(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     ilevel=input_array(1)
     call restrict_residual(pst%s,ilevel)
  endif

end subroutine r_restrict_residual

subroutine restrict_residual(s,ifinelevel)
  use amr_parameters, only: dp,nvector,nhilbert,ndim,twondim,twotondim
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils_p
  use cache_commons
  use hilbert
  use hash
  implicit none
  type(ramses_t)::s
  integer, intent(in) :: ifinelevel

  ! Restrict the residual of the fine level (stored in grid(ichild)%f(ind,1))
  ! into the rhs of the coarse level (stored in grid(igrid)%f(icell,2))
  ! For interior coarse cell only (we need the mask stored in grid(igrid)%f(icell,3))
  
  integer :: ichild, ind
  integer :: igrid, icell
  real(dp) :: dtwotondim = (twotondim)
  integer(kind=8),dimension(0:ndim) :: hash_key
  type(oct),pointer::gridp

  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Set rhs to zero in coarse cells
  do igrid=m%head_mg(ifinelevel-1),m%tail_mg(ifinelevel-1)
     do ind=1,twotondim
        m%grid(igrid)%f(ind,2)=0.0
     end do
  end do

  hash_key(0)=ifinelevel

  call open_cache(s,operation_restrict_res,domain_decompos_mg)
  
  ! Loop over grids
  do ichild=m%head_mg(ifinelevel),m%tail_mg(ifinelevel)
     
     ! Loop over cells
     do ind=1,twotondim        
        
        ! Is fine cell masked?
        if(m%grid(ichild)%f(ind,3)<=0d0)cycle

        hash_key(1:ndim)=m%grid(ichild)%ckey(1:ndim)
        
        ! Get parent cell using read-write cache
        call get_parent_cell_p(s,hash_key,m%mg_dict,gridp,icell,.true.,.true.)
        
        ! Is coarse cell masked?
        if(gridp%f(icell,3)<=0d0)cycle
        
        ! Stack fine cell residual in coarse cell rhs
        gridp%f(icell,2)=gridp%f(icell,2)+m%grid(ichild)%f(ind,1)/dtwotondim
        
     end do
  end do
  
  call close_cache(s,m%mg_dict)

  end associate
  
end subroutine restrict_residual

subroutine pack_fetch_restrict_res(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp(ind)=grid%f(ind,3)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_restrict_res

subroutine unpack_fetch_restrict_res(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

  msg=transfer(msg_array,msg)
  
#ifdef GRAV
  do ind=1,twotondim
     grid%f(ind,3)=msg%realdp(ind)
  end do
#endif

end subroutine unpack_fetch_restrict_res

subroutine init_flush_restrict_res(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  
#ifdef GRAV
  do ind=1,twotondim
     grid%f(ind,2)=0.0d0
  end do
#endif
  
end subroutine init_flush_restrict_res

subroutine pack_flush_restrict_res(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp(ind)=grid%f(ind,2)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_restrict_res

subroutine unpack_flush_restrict_res(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

  msg=transfer(msg_array,msg)
  
#ifdef GRAV
  do ind=1,twotondim
     grid%f(ind,2)=grid%f(ind,2)+msg%realdp(ind)
  end do
#endif

end subroutine unpack_flush_restrict_res

! ------------------------------------------------------------------------
! Interpolation and correction
! ------------------------------------------------------------------------

recursive subroutine r_interpolate_and_correct(pst,input_array,input_size,output_array,output_size)
  use amr_parameters, only: twotondim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_INTERPOLATE_AND_CORRECT,pst%iUpper+1,input_size,output_size,input_array)
     call r_interpolate_and_correct(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     ilevel=input_array(1)
     call interpolate_and_correct(pst%s,ilevel)
  endif

end subroutine r_interpolate_and_correct

subroutine interpolate_and_correct(s,ifinelevel)
  use amr_parameters, only: dp,nvector,nhilbert,ndim,twondim,twotondim,threetondim
  use amr_commons, only: nbor,oct
  use ramses_commons, only: ramses_t
  use nbors_utils_p
  use cache_commons
  use hilbert
  use hash
  implicit none
  type(ramses_t)::s
  integer, intent(in) :: ifinelevel
  
  ! Interpolate the solution of the coarse level (stored in grid(igrid)%phi(icell))
  ! and corrct the solutionn of the fine level (stored in grid(ichild)%phi(ind))
  
  integer(kind=8),dimension(0:ndim) :: hash_key
  integer,dimension(1:threetondim) :: igrid_nbor,ind_nbor
  type(nbor),dimension(1:threetondim) :: grid_nbor
  integer  :: ichild, ind
  real(dp) :: aa, bb, cc, dd, coeff
  real(dp), dimension(1:8)     :: bbb
  integer,  dimension(1:8,1:8) :: ccc
  integer::ind_average,ind_father
  integer::igrid_nbr,ind_nbr
  real(dp),dimension(1:twotondim)::corr
  type(oct),pointer::gridp
  
  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Local constants
  aa = 1.0D0/4.0D0**ndim
  bb = 3.0D0*aa
  cc = 9.0D0*aa
  dd = 27.D0*aa 
  bbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)
  
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)
  
  call open_cache(s,operation_phi,domain_decompos_mg)

  hash_key(0)=ifinelevel

  ! Loop over grids
  do ichild=m%head_mg(ifinelevel),m%tail_mg(ifinelevel)

     ! For fine level, correction is interpolated from coarser level solution
     hash_key(1:ndim)=m%grid(ichild)%ckey(1:ndim)
     
     ! Get 3**ndim neighbouring parent cell using a read-only cache
     call get_threetondim_nbor_parent_cell_p(s,hash_key,m%mg_dict,grid_nbor,ind_nbor,.false.,.true.)
     
     ! Loop over cells
     do ind=1,twotondim

        ! Set correction to zero
        corr(ind)=0d0

        ! Fine cell is masked as "outside": no correction
        if(m%grid(ichild)%f(ind,3)<=0.0)cycle

        ! Loop over relevant parent cells
        do ind_average=1,twotondim
           ind_father=ccc(ind_average,ind)
           coeff=bbb(ind_average)
           gridp=>grid_nbor(ind_father)%p
           ind_nbr=ind_nbor(ind_father)
           if (associated(gridp)) then
              corr(ind)=corr(ind)+coeff*gridp%phi(ind_nbr)
           endif
        end do

     end do
     ! End loop over cells
     
     do ind=1,threetondim
        call unlock_cache_p(s,grid_nbor(ind)%p)
     end do
     
     ! Add correction to fine level solution
     do ind=1,twotondim
        m%grid(ichild)%phi(ind)=m%grid(ichild)%phi(ind)+corr(ind)
     end do

  end do
  ! End loop over grids

  call close_cache(s,m%mg_dict)

  end associate
  
end subroutine interpolate_and_correct

subroutine pack_fetch_phi(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp(ind)=grid%phi(ind)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_phi

subroutine unpack_fetch_phi(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

  msg=transfer(msg_array,msg)
  
#ifdef GRAV
  do ind=1,twotondim
     grid%phi(ind)=msg%realdp(ind)
  end do
#endif

end subroutine unpack_fetch_phi

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Flag settings used to speed-up the sweeps
! ------------------------------------------------------------------------

recursive subroutine r_set_scan_flag(pst,input_array,input_size,output_array,output_size)
  use ramses_commons, only: pst_t
  use mdl_parameters
  use hash
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel,ifine

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_SET_SCAN_FLAG,pst%iUpper+1,input_size,output_size,input_array)
     call r_set_scan_flag(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     ilevel=input_array(1)
     ifine=input_array(2)
     if(ifine==ilevel)then
        call set_scan_flag(pst%s,pst%s%m%grid_dict,ifine)
     else
        call set_scan_flag(pst%s,pst%s%m%mg_dict,ifine)
     endif
  endif

end subroutine r_set_scan_flag

subroutine set_scan_flag(s,hash_dict,ilevel)
  use amr_parameters, only: dp,nvector,nhilbert,ndim,twondim,twotondim,threetondim
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils_p
  use cache_commons
  use hilbert
  use hash
  implicit none
  type(ramses_t)::s
  integer, intent(in) :: ilevel
  type(hash_table) :: hash_dict
  !
  integer :: ind, igrid, igridn, inbor, idim, id, ig
  integer, dimension(1:3,1:2,1:8)::iii, jjj
  real(dp),dimension(1:twotondim,0:twondim)::dis_nbor
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  real(dp)::dis_c
  type(oct),pointer::gridp
  
  associate(r=>s%r,g=>s%g,m=>s%m)

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  
  call open_cache(s,operation_scan,domain_decompos_mg)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=m%head_mg(ilevel),m%tail_mg(ilevel)

     ! Get central oct potential
     do ind=1,twotondim
        dis_nbor(ind,0)=m%grid(igrid)%f(ind,3)
     end do

     ! Get neighboring octs potential
     do inbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=m%grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=m%ckey_max(ilevel)-1
           if(hash_nbor(idim)==m%ckey_max(ilevel))hash_nbor(idim)=0
        enddo

        ! Get neighbouring grid using read-only cache
        call get_grid_p(s,hash_nbor,hash_dict,gridp,.false.,.true.)

        ! If grid exists, then copy into array
        if(associated(gridp))then
           do ind=1,twotondim
              dis_nbor(ind,inbor)=gridp%f(ind,3)
           end do

        ! Otherwise set to "outside"
        else
           do ind=1,twotondim
              dis_nbor(ind,inbor)=-1.0
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim

        ! Compute residual using 6 neighbors potential
        dis_c=m%grid(igrid)%f(ind,3)

        ! If cell is entirely inside, set flag tentatively to 0 (no scan needed)
        if(dis_c==1.0)then
           m%grid(igrid)%flag2(ind)=0

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 ! If one neighbour is outside, then scan needed
                 if(dis_nbor(id,ig)<=0.0)then
                    m%grid(igrid)%flag2(ind)=max(m%grid(igrid)%flag2(ind),1)
                 endif
              end do
           end do

         ! If cell is even partially outside, then scan needed
        else
           m%grid(igrid)%flag2(ind)=1           
        endif

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache(s,hash_dict)

  end associate
  
end subroutine set_scan_flag

subroutine pack_fetch_scan(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp(ind)=grid%f(ind,3)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_scan

subroutine unpack_fetch_scan(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

  msg=transfer(msg_array,msg)
  
#ifdef GRAV
  do ind=1,twotondim
     grid%f(ind,3)=msg%realdp(ind)
  end do
#endif

end subroutine unpack_fetch_scan

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Compute norm of residual 
! ------------------------------------------------------------------------

recursive subroutine r_cmp_residual_norm2(pst,input_array,input_size,output_array,output_size)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer,dimension(1:output_size)::next_output_array
  integer::ilevel
  real(kind=8)::norm2,next_norm2

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_CMP_RESIDUAL_NORM2,pst%iUpper+1,input_size,output_size,input_array)
     call r_cmp_residual_norm2(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size,next_output_array)
     norm2=transfer(output_array,norm2)
     next_norm2=transfer(next_output_array,next_norm2)
     norm2=norm2+next_norm2
     output_array=transfer(norm2,output_array)
  else
     ilevel=input_array(1)
     call cmp_residual_norm2(pst%s%r,pst%s%m,ilevel,norm2)
     output_array=transfer(norm2,output_array)
  endif

end subroutine r_cmp_residual_norm2

subroutine cmp_residual_norm2(r,m,ilevel, norm2)
  use amr_parameters, only: dp,ndim,twotondim
  use amr_commons, only: run_t,mesh_t
  implicit none
  type(run_t)::r
  type(mesh_t)::m
  integer,  intent(in)  :: ilevel
  real(kind=8), intent(out) :: norm2
  
  real(kind=8) :: dx2
  integer  :: ind, igrid
  
  ! Set constants
  dx2  = (r%boxlen/2**ilevel)**2
  norm2 = 0.0d0
  
  ! Loop over grids
  do igrid=m%head_mg(ilevel),m%tail_mg(ilevel)
     
     ! Loop over cells
     do ind=1,twotondim
        if(m%grid(igrid)%f(ind,3)<=0.0)cycle      ! Do not count masked cells
        norm2 = norm2 + m%grid(igrid)%f(ind,1)**2
     end do
     ! End loop over cells
     
  end do
  ! End loop over grids
  
  norm2 = dx2*norm2
  
end subroutine cmp_residual_norm2

#endif

