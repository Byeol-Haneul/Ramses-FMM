module fmm_init_refine_grid_module
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_fmm_init_refine_grid(pst)
  use ramses_commons, only: pst_t
  use flag_utils, only: m_flag_fine
  use init_refine_basegrid_module, only: r_collect_noct,r_noct_tot,r_noct_min,r_noct_max,r_noct_used_max
#ifdef GRAV
  use rho_fine_module, only: m_rho_fine
#endif
  use input_hydro_grafic_module, only: r_input_refmap_grafic
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine is the master procedure to set the  grid
  ! and fmm_initialize all cell-based variables within it.
  !--------------------------------------------------------------------
  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)

  allocate(m%head_mg(1:r%nlevelmax))
  allocate(m%tail_mg(1:r%nlevelmax))
  allocate(m%noct_mg(1:r%nlevelmax))

  if(r%verbose)write(*,*)'Entering fmm_init_refine_grid'

  write(*,*)'Building fmm_initial  grid at level ',r%levelmin

  ! Call recursive slave routine
  do ilevel = 1, r%levelmin - g%level_fmm_to_amr ! TODO: exception for this levelmin vs. level_fmm_to amr. 
    call r_fmm_init_refine_grid(pst, ilevel, 1)
    ! Get total, min and max grid count (only in master).
    call r_noct_tot(pst,ilevel,1,m%noct_tot(ilevel),2)
    call r_noct_min(pst,ilevel,1,m%noct_min(ilevel),1)
    call r_noct_max(pst,ilevel,1,m%noct_max(ilevel),1)
    call r_noct_used_max(pst,ilevel,1,m%noct_used_max,1)
  end do

  end associate

end subroutine m_fmm_init_refine_grid
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_fmm_init_refine_grid(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_FMM_INIT_REFINE_GRID,pst%iUpper+1,input_size,0,ilevel)
     call r_fmm_init_refine_grid(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call fmm_init_refine_grid(pst%s,ilevel)
  endif

end subroutine r_fmm_init_refine_grid
!################################################################
!################################################################
!################################################################
!################################################################
subroutine fmm_init_refine_grid(s,ilevel)
  use mdl_module, only: mdl_abort
  use amr_parameters, only: nhilbert,ndim,twotondim
  use ramses_commons, only: ramses_t
  use hilbert
  use hash, only: hash_setp, hash_is_clean, hash_stats
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !-------------------------------------------------------
  ! This routine builds a fully refined Cartesian grid
  ! at level ilevel. Always starts at levelmin.
  !-------------------------------------------------------
  logical::clean
  integer::i,igrid,ioct,ilev,istart,i1,j1,k1
  integer(kind=8)::ikey
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
  integer(kind=8),dimension(0:ndim)::hash_key
  integer(kind=8),dimension(1:nhilbert,1:s%r%nlevelmax)::key_ref
  integer(kind=8),dimension(1:nhilbert)::coarse_key
  integer,dimension(1:s%r%nlevelmax)::n_same,npatch

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ! Compute starting grid index at that level
  if(ilevel == 1)then
     m%ifree=m%noct_used+1 ! Jun-Young: start at index of the first free variable
     istart=m%ifree 
  else
     istart=m%tail_mg(ilevel-1)+1
  endif

  hk=0
  ! New grid in current level
  igrid=istart-1

  ! Loop over the Cartesian grid in Hilbert order
  do ikey=m%domain(ilevel)%b(1,g%myid-1), m%domain(ilevel)%b(1,g%myid)-1
     ! Compute Cartesian index from Hilbert index
     hk(1)=ikey
     ix=hilbert_reverse(hk,ilevel-1)
     if(ix(1).ge.m%box_ckey_min(1,ilevel).and.ix(1).lt.m%box_ckey_max(1,ilevel))then ! Jun-Young Lee if in the domain
#if NDIM>1
     if(ix(2).ge.m%box_ckey_min(2,ilevel).and.ix(2).lt.m%box_ckey_max(2,ilevel))then
#endif
#if NDIM>2
     if(ix(3).ge.m%box_ckey_min(3,ilevel).and.ix(3).lt.m%box_ckey_max(3,ilevel))then
#endif
        ! Insert new grid in main array
        igrid=igrid+1
        if(igrid.GT.r%ngridmax)then
           write(*,*)'No more free memory'
           write(*,*)'Increase ngridmax'
           call mdl_abort(mdl)
        end if
        if(igrid==istart)m%head_mg(ilevel)=istart
        m%tail_mg(ilevel)=igrid
        m%noct_mg(ilevel)=m%noct_mg(ilevel)+1
        m%noct(ilevel)=m%noct(ilevel)+1
        m%noct_used=m%noct_used+1
        m%grid(igrid)%lev=ilevel
        m%grid(igrid)%ckey(1:ndim)=int(ix(1:ndim),kind=4)
        m%grid(igrid)%hkey(1:nhilbert)=hk(1:nhilbert)
        m%grid(igrid)%refined(1:twotondim)=.false.
        ! Insert new grid in hash table
        hash_key(0)=ilevel
        hash_key(1:ndim)=ix(1:ndim)
        call hash_setp(m%mg_dict,hash_key,m%grid(igrid)) !Jun-Young grid_fmm
        !write(*,'(A,I0,A,I0,A,I0)') "[ID:", g%myid, "]: inserted grid at ilevel:", ilevel, " at igrid:", igrid
     endif
#if NDIM>1
     endif
#endif
#if NDIM>2
     endif
#endif
  end do

  !-----------
  ! Super-octs
  !-----------
  do i=1,ilevel
     npatch(i)=twotondim**i
  end do
  ilev=ilevel
  n_same=0
  key_ref=0
  key_ref(1,1:r%nlevelmax)=-1
  do ioct=m%head_mg(ilev),m%tail_mg(ilev)
     m%grid(ioct)%superoct=1
     coarse_key(1:nhilbert)=m%grid(ioct)%hkey(1:nhilbert)
     do i=1,MIN(ilev-1,r%nsuperoct)
        coarse_key(1:nhilbert)=coarsen_key(coarse_key(1:nhilbert),ilev-1) ! ilev-1 used to speed up only
        if(eq_keys(coarse_key(1:nhilbert),key_ref(1:nhilbert,i)))then
           n_same(i)=n_same(i)+1
        else
           n_same(i)=1
           key_ref(1:nhilbert,i)=coarse_key(1:nhilbert)
        endif
        if(n_same(i).EQ.npatch(i))then
           m%grid(ioct-npatch(i)+1:ioct)%superoct=npatch(i)
        endif
     end do
  end do

  !---------------------
  ! Clean and dirty octs
  !---------------------
  m%head_clean(ilev)=1
  m%head_dirty(ilev)=1
  m%noct_clean(ilev)=0
  m%noct_dirty(ilev)=0
  hash_key(0)=ilev
  do ioct=m%head_mg(ilev),m%tail_mg(ilev)
     clean=.true.
#if NDIM>2
     do k1=-1,1
     hash_key(3)=m%grid(ioct)%ckey(3)+k1
#endif
#if NDIM>1
     do j1=-1,1
     hash_key(2)=m%grid(ioct)%ckey(2)+j1
#endif
     do i1=-1,1
        hash_key(1)=m%grid(ioct)%ckey(1)+i1
        clean=clean.and.hash_is_clean(m%mg_dict,hash_key)
     end do
#if NDIM>1
     end do
#endif
#if NDIM>2
     end do
#endif
     if(clean)then
        m%indx_clean(m%head_clean(ilev)+m%noct_clean(ilev))=ioct
        m%noct_clean(ilev)=m%noct_clean(ilev)+1
     else
        m%indx_dirty(m%head_dirty(ilev)+m%noct_dirty(ilev))=ioct
        m%noct_dirty(ilev)=m%noct_dirty(ilev)+1
     endif
  end do

  end associate

end subroutine fmm_init_refine_grid
!################################################################
!################################################################
!################################################################
!################################################################
end module fmm_init_refine_grid_module
