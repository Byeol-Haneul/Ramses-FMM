!################################################################
!################################################################
!################################################################
!################################################################

module flag_utils

contains

subroutine m_flag_fine(pst,ilevel,icount)
  use ramses_commons, only: pst_t
  use smooth_module, only: r_smooth_fine
  implicit none
  type(pst_t)::pst
  integer::ilevel,icount
  !---------------------------------------------------------------
  ! This master routine builds the refinement map at level ilevel.
  !---------------------------------------------------------------
  integer::iexpand
  integer::nflag_tot

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,mdl=>pst%s%mdl)
  
  if(ilevel==r%nlevelmax)return
  if(ilevel<r%levelmin)return
  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,111)ilevel
111 format('   Entering flag_fine for level ',I2)

  ! Step 1: initialize refinement map to minimal refinement rules
  call r_init_flag(pst,ilevel,1,nflag_tot,1)
  if(r%verbose)write(*,*) '  ==> end step 1',nflag_tot
  
  ! Step 2: make one cubic buffer around flagged cells,
  ! in order to enforce numerical rule.
  call r_smooth_fine(pst,ilevel,1,nflag_tot,1)
  if(r%verbose)write(*,*) '  ==> end step 2',nflag_tot

  ! Step 3: if cell satisfies user-defined physical citeria,
  ! then flag cell for refinement.
  call r_user_flag(pst,ilevel,1,nflag_tot,1)
  if(r%verbose)write(*,*) '  ==> end step 3',nflag_tot

  ! Step 4: make nexpand cubic buffers around flagged cells.
  do iexpand=1,r%nexpand(ilevel)
     call r_smooth_fine(pst,ilevel,1,nflag_tot,1)
  end do
  if(r%verbose)write(*,*) '  ==> end step 4',nflag_tot

  ! In case of adaptive time step ONLY, check for refinement rules
  ! and unflag cells that will not be refined.
  if(ilevel>r%levelmin)then
     if(icount<r%nsubcycle(ilevel-1))then
        call r_ensure_ref_rules(pst,ilevel,1)
     end if
  end if

  end associate
  
end subroutine m_flag_fine
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_init_flag(pst,ilevel,input_size,noct,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel,noct

  integer::next_noct
  integer::nflag
  integer::rID
  
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INIT_FLAG,pst%iUpper+1,input_size,output_size,ilevel)
     call r_init_flag(pst%pLower,ilevel,input_size,noct,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_noct)
     noct=noct+next_noct
  else
     call init_flag(pst%s,ilevel,nflag)
     noct=nflag
  endif

end subroutine r_init_flag
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flag(s,ilevel,nflag)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use cache_commons
  use nbors_utils_p
  implicit none
  type(ramses_t)::s
  integer::ilevel,nflag
  !-------------------------------------------
  ! This routine initialize the refinement map
  ! to a minimal state in order to satisfy the
  ! refinement rules.
  !-------------------------------------------
  integer::igrid,ichild,icell,ind
  logical::ok
  integer(kind=8),dimension(0:ndim)::hash_key
  type(oct),pointer::gridp

  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Initialize flag1 to 0 for ilevel grids
  g%nflag=0
  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim
        m%grid(igrid)%flag1(ind)=0
     end do
  end do
  !---------------------------------------------------------
  ! Set flag1 to 1 if cell is refined and  contains a 
  ! flagged son or a refined son.
  ! This ensures that refinement rules are satisfied.
  !---------------------------------------------------------
  call open_cache(s,operation_initflag,domain_decompos_amr)

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ichild=m%head(ilevel+1),m%tail(ilevel+1)
     hash_key(1:ndim)=m%grid(ichild)%ckey(1:ndim)
     call get_parent_cell_p(s,hash_key,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.false.)
     ok=.false.
     ! Loop over cells
     do ind=1,twotondim
        ok=(ok.or.(m%grid(ichild)%refined(ind)))
        ok=(ok.or.(m%grid(ichild)%flag1(ind)==1))
     end do
     if(ok)then
        gridp%flag1(icell)=1
        g%nflag=g%nflag+1
     endif
  end do

  call close_cache(s,m%grid_dict)

  nflag=g%nflag

  end associate
  
end subroutine init_flag
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_flag(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_int4)::msg

  do ind=1,twotondim
     msg%int4(ind)=grid%flag1(ind)
  end do

  msg_array=transfer(msg,msg_array)
  
end subroutine pack_fetch_flag
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine unpack_fetch_flag(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_int4)::msg

  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     grid%flag1(ind)=msg%int4(ind)
  end do
  
end subroutine unpack_fetch_flag
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine init_flush_initflag(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  
  do ind=1,twotondim
     grid%flag1(ind)=0
  end do

end subroutine init_flush_initflag
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine pack_flush_initflag(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_int4)::msg

  do ind=1,twotondim
     msg%int4(ind)=grid%flag1(ind)
  end do

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_initflag
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine unpack_flush_initflag(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_int4)::msg

  msg=transfer(msg_array,msg)
  
  do ind=1,twotondim
     grid%flag1(ind)=MAX(grid%flag1(ind),msg%int4(ind))
  end do

end subroutine unpack_flush_initflag
!###############################################################
!###############################################################
!###############################################################
!###############################################################
recursive subroutine r_user_flag(pst,ilevel,input_size,noct,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel,noct

  integer::next_noct
  integer::nflag
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_USER_FLAG,pst%iUpper+1,input_size,output_size,ilevel)
     call r_user_flag(pst%pLower,ilevel,input_size,noct,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_noct)
     noct=noct+next_noct
  else
     call user_flag(pst%s,ilevel,nflag)
     noct=nflag
  endif

end subroutine r_user_flag
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine user_flag(s,ilevel,nflag)
  use ramses_commons, only: ramses_t
  use hydro_flag_module, only: hydro_flag
  implicit none
  type(ramses_t)::s
  integer::ilevel,nflag
  ! -------------------------------------------------------------------
  ! This routine flag for refinement cells that satisfies
  ! some user-defined physical criteria at the level ilevel. 
  ! -------------------------------------------------------------------

  ! Refinement rules for the gravity solver
  if(s%r%poisson)call poisson_flag(s,ilevel)

  ! Refinement rules for the hydro solver
  if(s%r%hydro)call hydro_flag(s,ilevel)

  nflag=s%g%nflag

end subroutine user_flag
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_ensure_ref_rules(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::ilevel
  integer,VALUE::input_size

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_ENSURE_REF_RULES,pst%iUpper+1,input_size,0,ilevel)
     call r_ensure_ref_rules(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call ensure_ref_rules(pst%s,ilevel)
  endif

end subroutine r_ensure_ref_rules
!############################################################
!############################################################
!############################################################
!############################################################
subroutine ensure_ref_rules(s,ilevel)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use cache_commons
  use nbors_utils_p
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !-----------------------------------------------------------------
  ! This routine determines if all grids at level ilevel are 
  ! surrounded by 26 neighboring grids, in order to enforce the 
  ! strict refinement rule. 
  ! Used in case of adaptive time steps only.
  !-----------------------------------------------------------------
  integer::idim,ind,igrid,ichild
  integer::i1,j1,k1
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer(kind=8),dimension(0:ndim)::hash_nbor
  logical::ok
  type(oct),pointer::gridp

  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Integer constants
  i1min=0; i1max=0; j1min=0; j1max=0; k1min=0; k1max=0
#if NDIM>0
  i1max=2
#endif
#if NDIM>1
  j1max=2
#endif
#if NDIM>2
  k1max=2
#endif

  call open_cache(s,operation_smooth,domain_decompos_amr)

  hash_nbor(0)=ilevel
  do igrid=m%head(ilevel),m%tail(ilevel)
     
     ok=.true.

     ! Loop over 3x3x3 neighboring father cells
     do k1=k1min,k1max
        do j1=j1min,j1max
           do i1=i1min,i1max

              ! Compute neighboring grid Cartesian index
#if NDIM>0
              hash_nbor(1)=m%grid(igrid)%ckey(1)+i1-1
#endif
#if NDIM>1
              hash_nbor(2)=m%grid(igrid)%ckey(2)+j1-1
#endif
#if NDIM>2
              hash_nbor(3)=m%grid(igrid)%ckey(3)+k1-1
#endif
              ! Periodic boundary conditons
              do idim=1,ndim
                 if(hash_nbor(idim)<0)hash_nbor(idim)=m%ckey_max(ilevel)-1
                 if(hash_nbor(idim)==m%ckey_max(ilevel))hash_nbor(idim)=0
              enddo

              ! Get neighboring grid index
              call get_grid_p(s,hash_nbor,m%grid_dict,gridp,flush_cache=.false.,fetch_cache=.true.)
              ok=ok.and.(associated(gridp))

           end do
        end do
     end do

     if(.not. ok)then
        do ind=1,twotondim
           m%grid(igrid)%flag1(ind)=0
        end do
     end if

  end do

  call close_cache(s,m%grid_dict)

  end associate
  
end subroutine ensure_ref_rules
!############################################################
!############################################################
!############################################################
!############################################################
end module flag_utils
