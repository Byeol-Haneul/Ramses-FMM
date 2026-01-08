module upload_module
contains
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

  call r_upload_fine(pst,ilevel,1)

end subroutine m_upload_fine
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_upload_fine(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_UPLOAD_FINE,pst%iUpper+1,input_size,0,ilevel)
     call r_upload_fine(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
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
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
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
  integer,dimension(1:6,1:4)::hh
  real(kind=8)::average,ekin,erad,emag
  type(oct),pointer::gridp
  type(msg_upload_hydro_mflux_mhd)::dummy_upload

#ifdef HYDRO

  hh(1,1:4)=(/1,3,5,7/)
  hh(2,1:4)=(/2,4,6,8/)
  hh(3,1:4)=(/1,2,5,6/)
  hh(4,1:4)=(/3,4,7,8/)
  hh(5,1:4)=(/1,2,3,4/)
  hh(6,1:4)=(/5,6,7,8/)

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
#ifdef HYDRO
     do ind=1,twotondim
        if(m%grid(ioct)%refined(ind))then
           m%grid(ioct)%mflux(ind,1:6)=0.0d0
        endif
     end do
#endif
#ifdef MHD
     do ivar=1,6
        do ind=1,twotondim
           if(m%grid(ioct)%refined(ind))then
              m%grid(ioct)%bold(ind,ivar)=0.0
           endif
        end do
     end do
#endif
  end do

  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_upload)/32,&
                     pack=pack_fetch_upload, unpack=unpack_fetch_upload,&
                     init=init_flush_upload, flush=pack_flush_upload, combine=unpack_flush_upload)

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ioct=m%head(ilevel+1),m%tail(ilevel+1)

     ! Get parent cell and grid index
     hash_key(1:ndim)=m%grid(ioct)%ckey(1:ndim)
     call get_parent_cell(s,hash_key,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.false.)

     ! Average conservative variables
     do ivar=1,nvar
        average=0.0d0
        do ind=1,twotondim
           average=average+m%grid(ioct)%uold(ind,ivar)
        end do
        ! Scatter result to parent cell
        gridp%uold(icell,ivar)=average/dble(twotondim)
     end do

     ! Average face-centered time-integrated mass flux (book-keeping, like face-centered B)
     ! mflux(:,1:ndim)   : left faces
     ! mflux(:,4:3+ndim) : right faces
     do idim=1,ndim
        ! Left face in parent cell
        average=0.0d0
        do ind=1,twotondim/2
           average=average+m%grid(ioct)%mflux(hh(2*idim-1,ind),idim)
        end do
        gridp%mflux(icell,idim)=average/dble(twotondim/2)
        ! Right face in parent cell
        average=0.0d0
        do ind=1,twotondim/2
           average=average+m%grid(ioct)%mflux(hh(2*idim,ind),idim+3)
        end do
        gridp%mflux(icell,idim+3)=average/dble(twotondim/2)
     end do

     ! Average cell-centered magnetic field
#ifdef MHD
#if NDIM<3
     ! Average cell-centered Bz
     average=0.0d0
     do ind=1,twotondim
        average=average+m%grid(ioct)%bold(ind,3)
     end do
     ! Scatter result to parent cell
     gridp%bold(icell,3)=average/dble(twotondim)
     gridp%bold(icell,6)=average/dble(twotondim)
#if NDIM==1
     ! Average cell-centered Bx and By
     do idim=1,2
        average=0.0d0
        do ind=1,twotondim
           average=average+m%grid(ioct)%bold(ind,idim)
        end do
        ! Scatter result to parentcell
        gridp%bold(icell,idim)=average/dble(twotondim)
        gridp%bold(icell,idim+3)=average/dble(twotondim)
     end do
#endif
#endif
#endif

     ! Average face-centered magnetic field
#ifdef MHD
#if NDIM>1
     ! Loop over dimensions
     do idim=1,ndim
        ! Left magnetic field in parent cell
        average=0.0d0
        do ind=1,twotondim/2
           average=average+m%grid(ioct)%bold(hh(2*idim-1,ind),idim)
        end do
        gridp%bold(icell,idim)=average/dble(twotondim/2)
        ! Right magnetic field in parent cell
        average=0.0d0
        do ind=1,twotondim/2
           average=average+m%grid(ioct)%bold(hh(2*idim,ind),idim+3)
        end do
        gridp%bold(icell,idim+3)=average/dble(twotondim/2)
     end do
#endif
#endif

     ! Average internal energy instead of total energy
     if(r%interpol_var==1)then
        average=0.0d0
        do ind=1,twotondim
           ekin=0.0d0
           do idim=1,3
              ekin=ekin+0.5d0*m%grid(ioct)%uold(ind,idim+1)**2/max(dble(m%grid(ioct)%uold(ind,1)),r%smallr)
           end do
           emag=0.0d0
#ifdef MHD
           do idim=1,3
              emag=emag+0.125d0*(m%grid(ioct)%bold(ind,idim)+m%grid(ioct)%bold(ind,idim+3))
           end do
#endif
           erad=0.0d0
#if NENER>0
           do irad=1,nener
              erad=erad+m%grid(ioct)%uold(ind,5+irad)
           end do
#endif
           average=average+m%grid(ioct)%uold(ind,5)-ekin-erad-emag
        end do
        ! Scatter result to parent cell
        ekin=0.0d0
        do idim=1,3
           ekin=ekin+0.5d0*gridp%uold(icell,idim+1)**2/max(dble(gridp%uold(icell,1)),r%smallr)
        end do
        emag=0.0d0
#ifdef MHD
        do idim=1,3
           emag=emag+0.125d0*(gridp%bold(icell,idim)+gridp%bold(icell,idim+3))
        end do
#endif
        erad=0.0d0
#if NENER>0
        do irad=1,nener
           erad=erad+gridp%uold(icell,5+irad)
        end do
#endif
        gridp%uold(icell,5)=average/dble(twotondim)+ekin+erad+emag
     endif
  end do

  call close_cache(s,m%grid_dict)

  end associate

#endif

end subroutine upload_fine

!##########################################################################
! Fetch pack/unpack for upload_fine (explicitly includes mflux)
!##########################################################################
subroutine pack_fetch_upload(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_upload_hydro_mflux_mhd
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer::ind,ivar
  type(msg_upload_hydro_mflux_mhd)::msg

  do ind=1,twotondim
     msg%int4(ind)=merge(1,0,grid%refined(ind))
  end do

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        msg%realdp_hydro(ind,ivar)=grid%uold(ind,ivar)
     end do
  end do
  msg%realdp_mflux=grid%mflux
#endif

#ifdef MHD
  msg%realdp_mhd=grid%bold
#endif

  msg_array=transfer(msg,msg_array)
end subroutine pack_fetch_upload

subroutine unpack_fetch_upload(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_upload_hydro_mflux_mhd
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key
  integer::ind,ivar
  type(msg_upload_hydro_mflux_mhd)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     grid%refined(ind)= (msg%int4(ind)==1)
  end do

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        grid%uold(ind,ivar)=msg%realdp_hydro(ind,ivar)
     end do
  end do
  grid%mflux=msg%realdp_mflux
#endif

#ifdef MHD
  grid%bold=msg%realdp_mhd
#endif

end subroutine unpack_fetch_upload
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine init_flush_upload(grid,hash_key)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  type(oct)::grid
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  
  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        grid%uold(ind,ivar)=0.0d0
     end do
  end do
  grid%mflux=0.0d0
#endif
#ifdef MHD
  grid%bold=0.0d0
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
  use cache_commons, only: msg_upload_hydro_mflux_mhd
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_upload_hydro_mflux_mhd)::msg

  do ind=1,twotondim
     msg%int4(ind)=merge(1,0,grid%refined(ind))
  end do

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        msg%realdp_hydro(ind,ivar)=grid%uold(ind,ivar)
     end do
  end do
#endif

#ifdef HYDRO
  msg%realdp_mflux=grid%mflux
#endif

#ifdef MHD
  msg%realdp_mhd=grid%bold
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_upload
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine unpack_flush_upload(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_upload_hydro_mflux_mhd
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  type(msg_upload_hydro_mflux_mhd)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)
  
#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        if(grid%refined(ind))then
           grid%uold(ind,ivar)=grid%uold(ind,ivar)+msg%realdp_hydro(ind,ivar)
        endif
     end do
  end do
#endif

#ifdef HYDRO
  do ind=1,twotondim
     if(grid%refined(ind))then
        grid%mflux(ind,1:6)=grid%mflux(ind,1:6)+msg%realdp_mflux(ind,1:6)
     endif
  end do
#endif

#ifdef MHD
  do ivar=1,6
     do ind=1,twotondim
        if(grid%refined(ind))then
           grid%bold(ind,ivar)=grid%bold(ind,ivar)+msg%realdp_mhd(ind,ivar)
        endif
     end do
  end do
#endif

end subroutine unpack_flush_upload
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
end module upload_module
