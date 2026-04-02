module init_fmm_module
integer, parameter :: FMM_BUILD_STANDARD = 0
integer, parameter :: FMM_BUILD_MERGED = 1
integer, parameter :: FMM_MULTIPOLE_STANDARD = 0
integer, parameter :: FMM_MULTIPOLE_MERGED = 1
integer, parameter :: FMM_TREE_SOURCE_STANDARD = 0
integer, parameter :: FMM_TREE_SOURCE_MERGED = 1

type :: double_level_t
    integer::ilevel,ifine
    integer::mode = FMM_BUILD_STANDARD
end type double_level_t

type :: fmm_level_t
    integer::ilev,flev
    integer::mode = FMM_MULTIPOLE_STANDARD
end type fmm_level_t

type :: downward_level_t
    integer::ilev,jlev,flev ! to call trees built on target(i)/source(j) level amr leaf cells
    integer::mode = FMM_TREE_SOURCE_STANDARD
end type downward_level_t
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
     call init_fmm(pst%s%r,pst%s%m,pst%s%m_fmm_list(ilevel),ilevel)
  endif

end subroutine r_init_fmm

subroutine init_fmm(r,m,m_fmm,ilevel)
  use amr_parameters, only: nhilbert
  use amr_commons, only: run_t,mesh_t
  use hilbert
  implicit none
  type(run_t)::r
  type(mesh_t)::m,m_fmm
  integer::ilevel
  integer::ilev,idom

  ! Compute multigrid Hilbert key tick marks
  call m_fmm%domain(ilevel)%copy(m%domain(ilevel))
  do ilev=ilevel-1,1,-1
     call m_fmm%domain(ilev)%copy(m_fmm%domain(ilev+1))
     do idom=0,m_fmm%domain(ilev)%ncpu
        m_fmm%domain(ilev)%b(1:nhilbert,idom)=coarsen_key(m_fmm%domain(ilev+1)%b(1:nhilbert,idom),ilev)
     end do
  end do
end subroutine init_fmm

recursive subroutine r_build_fmm(pst,input,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(double_level_t)::input

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BUILD_FMM,pst%iUpper+1,input_size,0,input)
     call r_build_fmm(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
      if (input%mode == FMM_BUILD_MERGED) then
          call build_fmm_merged(pst%s,input%ilevel)
      else if(input%ifine==input%ilevel)then
          call build_fmm(pst%s,pst%s%m,pst%s%m_fmm_list(input%ilevel), input)
      else
          call build_fmm(pst%s,pst%s%m_fmm_list(input%ilevel),pst%s%m_fmm_list(input%ilevel),input)
      end if
  endif

end subroutine r_build_fmm

subroutine build_fmm(s, m, m_fmm, input)
  use mdl_module
  use amr_parameters, only: nhilbert, ndim, twotondim, multipole_size, taylor_size
  use ramses_commons, only: ramses_t
  use amr_commons, only: mesh_t
  use cache_commons
  use cache
  use nbors_utils
  use hilbert
  use hash
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  type(double_level_t),intent(in)::input

  type(mesh_t)::m,m_fmm
  integer::ifinelevel,icoarselevel,igrid,idim,ichild,grid_cpu,child_cpu,ind,ifather
  integer::i1,j1,k1,i1min,i1max,j1min,j1max,k1min,k1max
  integer(kind=8),dimension(0:ndim)::hash_key,hash_father,hash_child
  integer(kind=4),dimension(1:ndim)::cart_key
  integer(kind=8),dimension(1:nhilbert)::hk, hk_child
  integer(kind=8),dimension(1:ndim)::ix
  logical::in_domain,check_refinement,flag_all_refined
  type(msg_small_realdp)::dummy_small_realdp

  associate(r=>s%r,g=>s%g,mdl=>s%mdl)
  ifinelevel=input%ifine
  icoarselevel=ifinelevel-1
  
  m_fmm%ifree=m_fmm%noct_used+1
  m_fmm%head(icoarselevel)=m_fmm%ifree
  hash_father(0)=icoarselevel

  check_refinement = (input%ilevel == ifinelevel)

  if (g%ncpu > 1 .and. check_refinement) then
     i1min = -1; i1max = 1
     j1min = -1*(ndim/2); j1max = 1*(ndim/2)
     k1min = -1*(ndim/3); k1max = 1*(ndim/3)
  else
     i1min = 0; i1max = 0
     j1min = 0; j1max = 0
     k1min = 0; k1max = 0
  end if

  call open_cache(mdl, m_fmm, pack_size=storage_size(dummy_small_realdp)/32, &
       flush=pack_flush_build_fmm, combine=unpack_flush_build_fmm)

  ! Loop over fine grids
  do igrid=m%head(ifinelevel),m%tail(ifinelevel)

    hash_key(1:ndim)=m%grid(igrid)%ckey(1:ndim)

    ! prob not worry
    if (check_refinement) then
      flag_all_refined = all(m%grid(igrid)%refined)
    else
      flag_all_refined = .false.
    end if

    do k1=k1min,k1max
       do j1=j1min,j1max
          do i1=i1min,i1max
             if ((i1 == 0) .and. (j1 == 0) .and. (k1 == 0) .and. flag_all_refined) cycle

             hash_child(0) = ifinelevel
#if NDIM>0
             hash_child(1) = hash_key(1) + i1
#endif
#if NDIM>1
             hash_child(2) = hash_key(2) + j1
#endif
#if NDIM>2
             hash_child(3) = hash_key(3) + k1
#endif

             in_domain = .true.
             do idim = 1, ndim
                if (r%periodic(idim)) then
                   if (hash_child(idim) < m%box_ckey_min(idim,ifinelevel)) then
                      hash_child(idim) = m%box_ckey_max(idim,ifinelevel) - 1
                   end if
                   if (hash_child(idim) >= m%box_ckey_max(idim,ifinelevel)) then
                      hash_child(idim) = m%box_ckey_min(idim,ifinelevel)
                   end if
                end if
                in_domain = in_domain .and. hash_child(idim) .ge. m%box_ckey_min(idim,ifinelevel) &
                     &                  .and. hash_child(idim) .lt. m%box_ckey_max(idim,ifinelevel)
             end do
             if (.not. in_domain) cycle

             hash_father(1:ndim)=hash_child(1:ndim)/2

             in_domain = .true.
             do idim = 1, ndim
                in_domain = in_domain .and. hash_father(idim) .ge. m_fmm%box_ckey_min(idim,icoarselevel) &
                     &                  .and. hash_father(idim) .lt.  m_fmm%box_ckey_max(idim,icoarselevel)
             end do

             if(.not. in_domain) cycle

        ! Access hash table
        ifather=hash_getp(m_fmm%grid_dict,hash_father)

        ! If grid does not exist, create it in memory
        if(ifather<=0)then

          ! Compute Cartesian keys of new oct
          cart_key(1:ndim)=int(hash_father(1:ndim),kind=4)

          ! Compute Hilbert keys of new octs
          ix(1:ndim)=cart_key(1:ndim)
          hk(1:nhilbert)=hilbert_key(ix,icoarselevel-1)

          ! Use the unique domain owner for coarse FMM grids. After the
          ! Hilbert boundaries are coarsened, the lower/upper-bound test can
          ! become ambiguous across ranks and duplicate parent ownership.
          grid_cpu = m_fmm%domain(icoarselevel)%get_rank(hk)

          if ((i1 /= 0 .or. j1 /= 0 .or. k1 /= 0) .and. check_refinement) then
             cart_key(1:ndim)=int(hash_child(1:ndim),kind=4)
             ix(1:ndim)=cart_key(1:ndim)
             hk_child(1:nhilbert)=hilbert_key(ix,ifinelevel-1)
             child_cpu = m%domain(ifinelevel)%get_rank(hk_child)
             if (child_cpu == mdl_self(mdl) .or. grid_cpu /= mdl_self(mdl)) cycle
             cart_key(1:ndim)=int(hash_father(1:ndim),kind=4)
          end if

          if(grid_cpu == mdl_self(mdl))then

              ! Set grid index to a virtual grid in local main memory
              ichild=m_fmm%ifree
              ! Go to next main memory free line
              m_fmm%ifree=m_fmm%ifree+1
              if(m_fmm%ifree.GT.m_fmm%ngridmax)then
                write(*,*)'No more free memory'
                write(*,*)'in multigrid'
                write(*,*)'Increase ngridmax'
                call mdl_abort(mdl)
              end if
              ! Insert new grid in hash table
              call hash_setp(m_fmm%grid_dict,hash_father,ichild)

          else

              ! Otherwise, stage the parent grid in the cache for its owner.
              ! If next cache line is occupied, free it.
              if(m_fmm%occupied(m_fmm%free_cache))call destage(mdl,m_fmm%ngridmax+m_fmm%free_cache)
              ! Set grid index to a virtual grid in local cache memory
              ichild=m_fmm%ngridmax+m_fmm%free_cache
              m_fmm%occupied(m_fmm%free_cache)=.true.
              m_fmm%parent_cpu(m_fmm%free_cache)=grid_cpu
              m_fmm%dirty(m_fmm%free_cache)=.true.
              m_fmm%ghost_parent_grid(m_fmm%free_cache)=0
              m_fmm%ghost_parent_cell(m_fmm%free_cache)=0
              ! Go to next free cache line
              m_fmm%free_cache=m_fmm%free_cache+1
              m_fmm%ncache=m_fmm%ncache+1
              if(m_fmm%free_cache.GT.m_fmm%ncachemax)m_fmm%free_cache=1
              if(m_fmm%ncache.GT.m_fmm%ncachemax)m_fmm%ncache=m_fmm%ncachemax
              ! Insert new grid in hash table
              call hash_setp(m_fmm%grid_dict,hash_father,ichild)

          endif

          ! Set oct properties
          m_fmm%grid(ichild)%lev=icoarselevel
          m_fmm%grid(ichild)%ckey(1:ndim)=cart_key(1:ndim)
          m_fmm%grid(ichild)%hkey(1:nhilbert)=hk(1:nhilbert)
          m_fmm%grid(ichild)%refined(1:twotondim)=.true.
          m_fmm%grid(ichild)%superoct=1

          ! Set flag arrays
          m_fmm%flag1(1:twotondim,ichild)=0
          m_fmm%flag2(1:twotondim,ichild)=0

          ! Intitialize gravity variables
          do ind=1,twotondim
#ifdef FMM
              m_fmm%multipole(ind,1:multipole_size,ichild)=0
              m_fmm%taylor_coeff(ind,1:taylor_size,ichild)=0
#endif
          enddo

        end if

          end do
       end do
    end do
  end do
  ! End loop over grids

  call close_cache(mdl)

  ! Finalize level bounds after close_cache because COMBINER_CREATE may
  ! materialize owner-local parent grids during the cache drain.
  m_fmm%tail(icoarselevel)=m_fmm%ifree-1
  m_fmm%noct(icoarselevel)=m_fmm%tail(icoarselevel)-m_fmm%head(icoarselevel)+1
  m_fmm%noct_used=m_fmm%tail(icoarselevel)

  if(r%verbose) print *, "      <LEV>: ", icoarselevel, "| created: ", m_fmm%noct(icoarselevel)
  end associate

end subroutine build_fmm

subroutine build_fmm_merged(s,active_levelmin)
  use mdl_module
  use amr_parameters, only: nhilbert, ndim, twotondim, multipole_size, taylor_size
  use ramses_commons, only: ramses_t
  use amr_commons, only: mesh_t
  use cache_commons
  use cache
  use hilbert
  use hash
  implicit none

  type(ramses_t)::s
  integer, intent(in) :: active_levelmin
  type(mesh_t), pointer :: m_merged
  integer::flev, jlev, ioct, idim, ind, igrid_new, first_ifree, igrid
  integer(kind=8),dimension(0:ndim)::hash_key
  integer(kind=4),dimension(1:ndim)::cart_key
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
  integer :: grid_cpu
  logical::in_domain
  type(msg_small_realdp)::dummy_small_realdp

  associate(r=>s%r, mdl=>s%mdl)
  if (.not. associated(s%m_fmm_merged)) return

  m_merged => s%m_fmm_merged

  call reset_entire_hash(m_merged%grid_dict,.false.)
  m_merged%head=1
  m_merged%tail=0
  m_merged%noct=0
  m_merged%ifree=1
  m_merged%noct_tot=0
  m_merged%noct_used=0

  call init_fmm(r,s%m,m_merged,r%nlevelmax)
  call open_cache(mdl, m_merged, pack_size=storage_size(dummy_small_realdp)/32, &
       flush=pack_flush_build_fmm, combine=unpack_flush_build_fmm)

  do flev=r%bound_levelmin,r%nlevelmax
     first_ifree = m_merged%ifree
     hash_key(0) = flev

     do jlev=max(active_levelmin, flev+1),r%nlevelmax
        if (s%m_fmm_list(jlev)%tail(flev) < s%m_fmm_list(jlev)%head(flev)) cycle

        do ioct=s%m_fmm_list(jlev)%head(flev),s%m_fmm_list(jlev)%tail(flev)
           hash_key(1:ndim)=s%m_fmm_list(jlev)%grid(ioct)%ckey(1:ndim)

           in_domain = .true.
           do idim = 1, ndim
              in_domain = in_domain .and. hash_key(idim) .ge. m_merged%box_ckey_min(idim,flev) &
                   &                .and. hash_key(idim) .lt. m_merged%box_ckey_max(idim,flev)
           end do
           if(.not. in_domain) cycle

           igrid = hash_getp(m_merged%grid_dict,hash_key)
           if(igrid>0) cycle

           cart_key(1:ndim)=int(hash_key(1:ndim),kind=4)
           ix(1:ndim)=cart_key(1:ndim)
           hk(1:nhilbert)=hilbert_key(ix,flev-1)

           grid_cpu = m_merged%domain(flev)%get_rank(hk)

           if(grid_cpu == mdl_self(mdl)) then
              igrid_new=m_merged%ifree
              m_merged%ifree=m_merged%ifree+1
              if(m_merged%ifree.GT.m_merged%ngridmax)then
                 write(*,*)'No more free memory'
                 write(*,*)'in merged fmm tree'
                 write(*,*)'Increase ngridmax'
                 call mdl_abort(mdl)
              end if
              call hash_setp(m_merged%grid_dict,hash_key,igrid_new)
           else
              if(m_merged%occupied(m_merged%free_cache))call destage(mdl,m_merged%ngridmax+m_merged%free_cache)
              igrid_new=m_merged%ngridmax+m_merged%free_cache
              m_merged%occupied(m_merged%free_cache)=.true.
              m_merged%parent_cpu(m_merged%free_cache)=grid_cpu
              m_merged%dirty(m_merged%free_cache)=.true.
              m_merged%ghost_parent_grid(m_merged%free_cache)=0
              m_merged%ghost_parent_cell(m_merged%free_cache)=0
              m_merged%free_cache=m_merged%free_cache+1
              m_merged%ncache=m_merged%ncache+1
              if(m_merged%free_cache.GT.m_merged%ncachemax)m_merged%free_cache=1
              if(m_merged%ncache.GT.m_merged%ncachemax)m_merged%ncache=m_merged%ncachemax
              call hash_setp(m_merged%grid_dict,hash_key,igrid_new)
           end if

           m_merged%grid(igrid_new)%lev=flev
           m_merged%grid(igrid_new)%ckey(1:ndim)=cart_key(1:ndim)
           m_merged%grid(igrid_new)%hkey(1:nhilbert)=hk(1:nhilbert)
           m_merged%grid(igrid_new)%refined(1:twotondim)=.true.
           m_merged%grid(igrid_new)%superoct=1

           m_merged%flag1(1:twotondim,igrid_new)=0
           m_merged%flag2(1:twotondim,igrid_new)=0
           do ind=1,twotondim
#ifdef FMM
              m_merged%multipole(ind,1:multipole_size,igrid_new)=0.0D0
              m_merged%taylor_coeff(ind,1:taylor_size,igrid_new)=0.0D0
#endif
           end do
        end do
     end do

     if (m_merged%ifree>first_ifree) then
        m_merged%head(flev)=first_ifree
        m_merged%tail(flev)=m_merged%ifree-1
        m_merged%noct(flev)=m_merged%tail(flev)-m_merged%head(flev)+1
     else
        m_merged%head(flev)=1
        m_merged%tail(flev)=0
        m_merged%noct(flev)=0
     end if
  end do

  call close_cache(mdl)

  m_merged%noct_used=m_merged%ifree-1
  if(r%verbose) print *, "      <MERGED TREE> created: ", m_merged%noct_used
  end associate
end subroutine build_fmm_merged

subroutine pack_flush_build_fmm(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_small_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp(ind)=0.0d0
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_build_fmm

subroutine unpack_flush_build_fmm(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_small_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::idim,ind
  type(msg_small_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     mesh%grid(igrid)%refined(ind)=.true.
  end do

  do ind=1,twotondim
#ifdef FMM
     mesh%multipole(ind,:,igrid)=0
     mesh%taylor_coeff(ind,:,igrid)=0
#endif
  end do
end subroutine unpack_flush_build_fmm
#endif
!################################################################
!################################################################
!################################################################
!################################################################
end module init_fmm_module
