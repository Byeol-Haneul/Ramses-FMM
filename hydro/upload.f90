!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_upload_fine(pst,ilevel)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine is the master procedure to upload HYDRO variables
  ! from level ilevel+1 to ilevel (averaging down or restriction).
  !--------------------------------------------------------------------
  integer::dummy

  if(ilevel==pst%s%r%nlevelmax)return
  if(pst%s%m%noct_tot(ilevel)==0)return
  if(pst%s%m%noct_tot(ilevel+1)==0)return
  if(pst%s%r%verbose)write(*,111)ilevel
111 format(' Entering upload_fine for level',i2)

  call r_upload_fine(pst,ilevel,1,dummy,0)

end subroutine m_upload_fine
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_upload_fine(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_UPLOAD_FINE,pst%iUpper+1,input_size,output_size,input_array)
     call r_upload_fine(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     ilevel=input_array(1)
     call upload_fine(pst%s,ilevel)
  endif

end subroutine r_upload_fine
!###########################################################
!########################################################### 
!###########################################################
!###########################################################
subroutine upload_fine(s,ilevel)
  use mdl_module
  use hydro_parameters, only: nvar,nener
  use amr_parameters, only: dp,ndim,twotondim
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils_p
  use cache_commons
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !----------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the hydro variables.
  !----------------------------------------------------------------------
#if NENER>0
  integer::irad
#endif
  integer::ioct,ind,ivar,icell,idim
  integer(kind=8),dimension(0:ndim)::hash_key
  real(dp)::average,ekin,erad
  type(oct),pointer::gridp

#ifdef HYDRO

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ! Set conservative variable to zero in refined cells
  do ioct=m%head(ilevel),m%tail(ilevel)
     do ivar=1,nvar
        do ind=1,twotondim
           if(m%grid(ioct)%refined(ind))then
              m%grid(ioct)%uold(ind,ivar)=0.0
           endif
        end do
     end do
  end do

  call open_cache(s,operation_upload,domain_decompos_amr)

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ioct=m%head(ilevel+1),m%tail(ilevel+1)

     ! Get cell and grid index
     hash_key(1:ndim)=m%grid(ioct)%ckey(1:ndim)
     call get_parent_cell_p(s,hash_key,m%grid_dict,gridp,icell,.true.,.false.)

     ! Average conservative variables
     do ivar=1,nvar
        average=0.0d0
        do ind=1,twotondim
           average=average+m%grid(ioct)%uold(ind,ivar)
        end do
        ! Scatter result to cell
        gridp%uold(icell,ivar)=average/dble(twotondim)
     end do

     ! Average internal energy instead of total energy
     if(r%interpol_var==1 .or. r%interpol_var==2)then
        average=0.0d0
        do ind=1,twotondim
           ekin=0.0d0
           do idim=1,ndim
              ekin=ekin+0.5d0*m%grid(ioct)%uold(ind,idim+1)**2/max(m%grid(ioct)%uold(ind,1),r%smallr)
           end do
           erad=0.0d0
#if NENER>0
           do irad=1,nener
              erad=erad+m%grid(ioct)%uold(ind,ndim+2+irad)
           end do
#endif
           average=average+m%grid(ioct)%uold(ind,ndim+2)-ekin-erad
        end do
        ! Scatter result to cell
        ekin=0.0d0
        do idim=1,ndim
           ekin=ekin+0.5d0*gridp%uold(icell,idim+1)**2/max(gridp%uold(icell,1),r%smallr)
        end do
        erad=0.0d0
#if NENER>0
        do irad=1,nener
           erad=erad+gridp%uold(icell,ndim+2+irad)
        end do
#endif
        gridp%uold(icell,ndim+2)=average/dble(twotondim)+ekin+erad
     endif
  end do

  call close_cache(s,m%grid_dict)

  ! Set conservative variable to zero in refined cells
  do ioct=m%head(ilevel),m%tail(ilevel)
     do ivar=1,nvar
        do ind=1,twotondim
           if(m%grid(ioct)%uold(ind,1)==0.0)then
              write(*,*)'Sorry zero density cell'
              write(*,*)m%grid(ioct)%uold(ind,1:nvar)
              write(*,*)m%grid(ioct)%refined(ind)
              call mdl_abort(mdl)
           endif
        end do
     end do
  end do

  end associate

#endif

end subroutine upload_fine
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine init_flush_upload(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  
#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        grid%uold(ind,ivar)=0.0d0
     end do
  end do
#endif
  
end subroutine init_flush_upload
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine pack_flush_upload(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_realdp)::msg

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        msg%realdp(ind,ivar)=grid%uold(ind,ivar)
     end do
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_upload
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine unpack_flush_upload(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_realdp)::msg

  msg=transfer(msg_array,msg)
  
#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        if(grid%refined(ind))then
           grid%uold(ind,ivar)=grid%uold(ind,ivar)+msg%realdp(ind,ivar)
        endif
     end do
  end do
#endif

end subroutine unpack_flush_upload
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
