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
  use marshal, only: pack_fetch_flag, unpack_fetch_flag
  use cache_commons
  use cache
  use nbors_utils
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
  type(msg_int4)::dummy_int4

  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Initialize flag1 to 0 for ilevel grids
  g%nflag=0
  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim
        m%grid(igrid)%flag1(ind)=0
     end do
  end do
  !---------------------------------------------------------
  ! Set flag1 to 1 if cell is refined and contains a
  ! flagged son or a refined son.
  ! This ensures that refinement rules are satisfied.
  !---------------------------------------------------------
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_int4)/32,&
                pack=pack_fetch_flag,unpack=unpack_fetch_flag,&
                init=init_flush_initflag, flush=pack_flush_initflag, combine=unpack_flush_initflag)

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ichild=m%head(ilevel+1),m%tail(ilevel+1)
     hash_key(1:ndim)=m%grid(ichild)%ckey(1:ndim)
     call get_parent_cell(s,hash_key,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.false.)
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
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine init_flush_initflag(grid,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  type(oct)::grid
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  
  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
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
subroutine unpack_flush_initflag(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_int4)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
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
  use rt_flag_module, only: rt_flag
  implicit none
  type(ramses_t)::s
  integer::ilevel,nflag
  ! -------------------------------------------------------------------
  ! This routine flag for refinement cells that satisfies
  ! some user-defined physical criteria at the level ilevel. 
  ! -------------------------------------------------------------------
  integer::level_lock

  ! Unlock levels progressively
  if(s%r%cosmo.and.s%r%aexp_lock_refine>0d0)then
     if(ilevel.GT.s%g%nlevelmax_part+3)then
        level_lock=s%r%nlevelmax+int(log(s%g%aexp/s%r%aexp_lock_refine/2.0d0)/log(2d0))
        if(ilevel.GE.level_lock)then
           nflag=s%g%nflag
           return
        endif
     endif
  endif

  ! Refinement rules for the gravity solver
  if(s%r%poisson)call poisson_flag(s,ilevel)

  ! Refinement rules for the hydro solver
  if(s%r%hydro)call hydro_flag(s,ilevel)

  ! Refinement rules for the radiative transfer solver
  if(s%r%rt)call rt_flag(s,ilevel)
  ! Refinement rules around sink particles
  if(s%r%sink.and.s%r%sink_refine)call sink_flag(s,s%sink,ilevel)

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
  use marshal, only: pack_fetch_flag, unpack_fetch_flag
  use boundaries, only: init_bound_flag
  use cache_commons
  use cache
  use nbors_utils
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
  type(msg_int4)::dummy_int4

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

  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_int4)/32,&
                pack=pack_fetch_flag,unpack=unpack_fetch_flag,&
                bound=init_bound_flag)

  hash_nbor(0)=ilevel
  do igrid=m%head(ilevel),m%tail(ilevel)
     
     ok=.true.

     ! Loop over 3x3x3 neighboring father cells
     do k1=k1min,k1max
        do j1=j1min,j1max
           do i1=i1min,i1max

              ! Compute neighboring grid Cartesian index
#if NDIM>0
              hash_nbor(1)=m%grid(igrid)%ckey(1)+i1-3*(i1/2)
#endif
#if NDIM>1
              hash_nbor(2)=m%grid(igrid)%ckey(2)+j1-3*(j1/2)
#endif
#if NDIM>2
              hash_nbor(3)=m%grid(igrid)%ckey(3)+k1-3*(k1/2)
#endif
              ! Periodic boundary conditions
              do idim=1,ndim
                 if(r%periodic(idim))then
                    if(hash_nbor(idim)<m%box_ckey_min(idim,ilevel))hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
                    if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel))hash_nbor(idim)=m%box_ckey_min(idim,ilevel)
                 endif
              enddo

              ! Get neighboring grid index
              call get_grid(s,hash_nbor,m%grid_dict,gridp,flush_cache=.false.,fetch_cache=.true.)
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
subroutine sink_flag(s,p,ilevel)
   use amr_parameters, only: ndim,twotondim,dp
   use amr_commons, only: nbor,oct
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   use nbors_utils
   use cache_commons
   use cache
   use marshal, only: pack_fetch_flag, unpack_fetch_flag
   !use flag_utils, only: init_flush_initflag, pack_flush_initflag, unpack_flush_initflag
   use hilbert
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ilevel
   !==================================================================
   ! This routine flag for refinement cells that are close enough
   ! from sink particles.
   !==================================================================
   real(kind=8),dimension(:,:),allocatable::x_nei
   integer,dimension(1:ndim)::ckey,ckey_ref,ckey_nbor
   integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
   integer::i,j,k,ipart,icellp,icelln,ind,idim,i_nei,n_nei,nrad
   integer,dimension(1:ndim)::ix
   real(kind=8)::dx_loc,vol_loc,x,y,z,rrad,rr
   real(kind=8),dimension(1:3)::xcen,xnei
   type(oct),pointer::gridp,gridn
   type(msg_int4)::dummy_int4
 
#ifdef HYDRO
#if NDIM==3
   associate(r=>s%r,g=>s%g,m=>s%m)
 
   ! Mesh spacing in that level
   dx_loc = r%boxlen / 2**ilevel 
   vol_loc = dx_loc**ndim
 
   ! Compute number of cells within sink sphere
   nrad = r%sink_b_spline_order
   rrad = dble(nrad)
   n_nei = 0
   do k = -nrad, nrad
      z = dble(k) / dble(nrad) * rrad
      do j = -nrad, nrad
         y = dble(j) / dble(nrad) * rrad
         do i = -nrad, nrad
            x = dble(i) / dble(nrad) * rrad
            rr = sqrt(dble(x*x+y*y+z*z))
            if(rr .le. rrad) n_nei = n_nei + 1
         enddo
      enddo
   enddo
   allocate(x_nei(1:ndim,1:n_nei))
   i_nei = 0
   do k = -nrad, nrad
      z = dble(k) / dble(nrad) * rrad
      do j = -nrad, nrad
         y = dble(j) / dble(nrad) * rrad
         do i = -nrad, nrad
            x = dble(i) / dble(nrad) * rrad
            rr = sqrt(dble(x*x+y*y+z*z))
            if(rr .le. rrad)then
               i_nei = i_nei + 1
               x_nei(1, i_nei) = x
               x_nei(2, i_nei) = y
               x_nei(3, i_nei) = z
            endif
         enddo
      enddo
   enddo
 
   ! Open cache for array uold (fetch) and unew (flush)
   call open_cache(s, table=m%grid_dict, data_size=storage_size(m%grid(1))/32, &
                 hilbert=m%domain, pack_size=storage_size(dummy_int4)/32, &
                 pack=pack_fetch_flag, unpack=unpack_fetch_flag, &
                 init=init_flush_initflag, flush=pack_flush_initflag, combine=unpack_flush_initflag)
 
   ! Loop over sink particles at current level and finer levels
   hash_nbor(0) = ilevel+1
   do ipart = p%headp(ilevel), p%npart
 
      ! Sink sphere center in units of current level Cartesian coordinates
      xcen(1:ndim) = p%xp(ipart,1:ndim) / dx_loc
 
      ! Collect sink sphere sampling points
      do i_nei = 1, n_nei
 
         ! Compute neighboring cell coordinates
         xnei(1:ndim) = xcen(1:ndim) + x_nei(1:ndim, i_nei)
         ! Periodic boundary conditions
         do idim=1,ndim
            if(xnei(idim)<                0.0d0)xnei(idim)=xnei(idim)+m%ckey_max(ilevel+1)
            if(xnei(idim)>=m%ckey_max(ilevel+1))xnei(idim)=xnei(idim)-m%ckey_max(ilevel+1)
         end do
         ! Get neighboring cell at current level
         hash_nbor(1:ndim) = int(xnei(1:ndim))
         call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.false.)
         ! If missing, then cycle. This should never happens if sink_refine=.true.
         if(.not.associated(gridn))cycle
 
         ! Set refinement map flag1 to 1
         gridn%flag1(icelln)=1
 
      end do
      ! End loop over sink sphere sampling points
 
   end do
   ! End loop over particles
 
   call close_cache(s,m%grid_dict)
 
   deallocate(x_nei)
 
   end associate
#endif
#endif
end subroutine sink_flag
!############################################################
!############################################################
!############################################################
!############################################################
end module flag_utils
