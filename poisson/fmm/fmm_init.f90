module init_fmm_module
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
#ifdef GRAV
recursive subroutine r_init_fmm(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INIT_FMM,pst%iUpper+1,input_size,0,ilevel)
     call r_init_fmm(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call init_fmm(pst%s,ilevel)
  endif

end subroutine r_init_fmm
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_fmm(s,ilevel)
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
  integer::i,igrid,ioct,ilev,istart,i1,j1,k1,idom
  integer(kind=8)::ikey
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
  integer(kind=8),dimension(0:ndim)::hash_key
  integer(kind=8),dimension(1:nhilbert,1:s%r%nlevelmax)::key_ref
  integer(kind=8),dimension(1:nhilbert)::coarse_key
  integer,dimension(1:s%r%nlevelmax)::n_same,npatch

  associate(r=>s%r,g=>s%g,m=>s%m,m_mg=>s%m_mg,mdl=>s%mdl)

  call m_mg%domain(ilevel)%copy(m%domain(ilevel))
  do ilev=ilevel-1,1,-1
     call m_mg%domain(ilev)%copy(m_mg%domain(ilev+1))
     do idom=0,m_mg%domain(ilev)%ncpu
        m_mg%domain(ilev)%b(1:nhilbert,idom)=coarsen_key(m_mg%domain(ilev+1)%b(1:nhilbert,idom),ilev)
     end do
  end do
  
  do ilev=r%bound_levelmin, r%levelmin-r%level_fmm_to_amr
    ! Compute starting grid index at that level
    if(ilev == r%bound_levelmin)then
      m_mg%ifree=m_mg%noct_used+1 ! Jun-Young: start at index of the first free variable
      istart=m_mg%ifree 
    else
      istart=m_mg%tail(ilev-1)+1
    endif
    hk=0
    ! New grid in current level
    igrid=istart-1
    m_mg%head(ilev)=istart
    m_mg%tail(ilev)=igrid
    m_mg%noct(ilev)=0

    ! Loop over the Cartesian grid in Hilbert order
    do ikey=m_mg%domain(ilev)%b(1,g%myid-1), m_mg%domain(ilev)%b(1,g%myid)-1
      ! Compute Cartesian index from Hilbert index
      hk(1)=ikey
      ix=hilbert_reverse(hk,ilev-1)
      if(ix(1).ge.m_mg%box_ckey_min(1,ilev).and.ix(1).lt.m_mg%box_ckey_max(1,ilev))then ! Jun-Young Lee if in the domain
#if NDIM>1
      if(ix(2).ge.m_mg%box_ckey_min(2,ilev).and.ix(2).lt.m_mg%box_ckey_max(2,ilev))then
#endif
#if NDIM>2
      if(ix(3).ge.m_mg%box_ckey_min(3,ilev).and.ix(3).lt.m_mg%box_ckey_max(3,ilev))then
#endif
          ! Insert new grid in main array
          igrid=igrid+1
          if(igrid.GT.m_mg%ngridmax)then
            write(*,*)'No more free memory'
            write(*,*)'Increase ngridmax'
            call mdl_abort(mdl)
          end if
          if(igrid==istart)m_mg%head(ilev)=istart
          m_mg%tail(ilev)=igrid
          m_mg%noct(ilev)=m_mg%noct(ilev)+1
          !m_mg%noct(ilev)=m_mg%noct(ilev)+1
          !m_mg%noct_used=m_mg%noct_used+1
          m_mg%grid(igrid)%lev=ilev
          m_mg%grid(igrid)%ckey(1:ndim)=int(ix(1:ndim),kind=4)
          m_mg%grid(igrid)%hkey(1:nhilbert)=hk(1:nhilbert)
          m_mg%grid(igrid)%refined(1:twotondim)=.true.
          ! Insert new grid in hash table
          hash_key(0)=ilev
          hash_key(1:ndim)=ix(1:ndim)
          call hash_setp(m_mg%grid_dict,hash_key,igrid)
      endif
#if NDIM>1
      endif
#endif
#if NDIM>2
      endif
#endif
    end do
  end do
  end associate
end subroutine init_fmm
#endif
!################################################################
!################################################################
!################################################################
!################################################################
end module init_fmm_module
