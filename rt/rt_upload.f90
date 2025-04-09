module rt_upload_module
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_rt_upload_fine(pst,ilevel)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine is the master procedure to upload HYDRO variables
  ! from level ilevel+1 to ilevel (averaging down or restriction).
  !--------------------------------------------------------------------
  if(ilevel==pst%s%r%nlevelmax)return
  if(pst%s%m%noct_tot(ilevel)==0)return
  if(pst%s%m%noct_tot(ilevel+1)==0)return
  if(pst%s%r%verbose)write(*,111)ilevel
111 format(' Entering upload_fine for level',i2)

  call r_rt_upload_fine(pst,ilevel,1)

end subroutine m_rt_upload_fine
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_rt_upload_fine(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RT_UPLOAD_FINE,pst%iUpper+1,input_size,0,ilevel)
     call r_rt_upload_fine(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call rt_upload_fine(pst%s,ilevel)
  endif

end subroutine r_rt_upload_fine
!###########################################################
!########################################################### 
!###########################################################
!###########################################################
subroutine rt_upload_fine(s,ilevel)
  use mdl_module
  use rt_parameters, only: nrtvar
  use amr_parameters, only: dp,ndim,twotondim
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use rt_flag_module, only: pack_fetch_rt,unpack_fetch_rt
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !----------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the RT variables.
  !----------------------------------------------------------------------
  integer::ioct,ind,ivar,icell
  integer(kind=8),dimension(0:ndim)::hash_key
  integer,dimension(1:6,1:4)::hh
  real(dp)::average
  type(oct),pointer::gridp
  type(msg_realdp)::dummy_realdp

  hh(1,1:4)=(/1,3,5,7/)
  hh(2,1:4)=(/2,4,6,8/)
  hh(3,1:4)=(/1,2,5,6/)
  hh(4,1:4)=(/3,4,7,8/)
  hh(5,1:4)=(/1,2,3,4/)
  hh(6,1:4)=(/5,6,7,8/)

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ! Set conservative variable to zero in refined cells
  do ioct=m%head(ilevel),m%tail(ilevel)
     do ivar=1,nrtvar
        do ind=1,twotondim
           if(m%grid(ioct)%refined(ind))then
#ifdef RT
              m%grid(ioct)%rtuold(ind,ivar)=0.0
#endif
           endif
        end do
     end do
  end do

  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_realdp)/32,&
                     pack=pack_fetch_rt, unpack=unpack_fetch_rt,&
                     init=init_flush_upload_rt, flush=pack_flush_upload_rt, combine=unpack_flush_upload_rt)

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ioct=m%head(ilevel+1),m%tail(ilevel+1)

     ! Get parent cell and grid index
     hash_key(1:ndim)=m%grid(ioct)%ckey(1:ndim)
     call get_parent_cell(s,hash_key,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.false.)

     ! Average conservative variables
     do ivar=1,nrtvar
        average=0.0d0
        do ind=1,twotondim
#ifdef RT
           average=average+m%grid(ioct)%rtuold(ind,ivar)
#endif
        end do
        ! Scatter result to parent cell
#ifdef RT
        gridp%rtuold(icell,ivar)=average/dble(twotondim)
        ! Rescale according to light speed difference between levels
        if (mod(ivar,ndim+1).eq.1) &
            gridp%rtuold(icell,ivar) = gridp%rtuold(icell,ivar) * g%rt_c(ilevel+1)/g%rt_c(ilevel)
#endif
     end do

  end do

  call close_cache(s,m%grid_dict)

  end associate

end subroutine rt_upload_fine
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine init_flush_upload_rt(grid,hash_key)
  use amr_parameters, only: ndim,twotondim
  use rt_parameters, only: nrtvar
  use amr_commons, only: oct
  type(oct)::grid
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  do ivar=1,nrtvar
     do ind=1,twotondim
#ifdef RT  
        grid%rtuold(ind,ivar)=0.0d0
#endif
     end do
  end do
end subroutine init_flush_upload_rt
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine pack_flush_upload_rt(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use rt_parameters, only: nrtvar
  use amr_commons, only: oct
  use cache_commons, only: msg_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_realdp)::msg
  do ivar=1,nrtvar
     do ind=1,twotondim
#ifdef RT
        msg%realdp_rt(ind,ivar)=grid%rtuold(ind,ivar)
#endif
     end do
  end do
  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_upload_rt
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine unpack_flush_upload_rt(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use rt_parameters, only: nrtvar
  use amr_commons, only: oct
  use cache_commons, only: msg_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  type(msg_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)
  do ivar=1,nrtvar
     do ind=1,twotondim
        if(grid%refined(ind))then
#ifdef RT
           grid%rtuold(ind,ivar)=grid%rtuold(ind,ivar)+msg%realdp_rt(ind,ivar)
#endif
        endif
     end do
  end do
end subroutine unpack_flush_upload_rt
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
end module rt_upload_module
