module init_fmm_module
type :: double_level_t
    integer::ilevel,ifine
end type double_level_t

type :: fmm_level_t
    integer::ilev,flev
end type fmm_level_t

type :: downward_level_t
    integer::ilev,jlev,flev ! to call trees built on target(i)/source(j) level amr leaf cells
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
      if(input%ifine==input%ilevel)then
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
  integer::ifinelevel,icoarselevel,igrid,idim,ichild,grid_cpu,ind,ifather
  integer(kind=8),dimension(0:ndim)::hash_key,hash_father
  integer(kind=4),dimension(1:ndim)::cart_key
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
  logical::in_rank,in_domain,check_refinement,flag_all_refined
  type(msg_small_realdp)::dummy_small_realdp

  associate(r=>s%r,g=>s%g,mdl=>s%mdl)
  ifinelevel=input%ifine
  icoarselevel=ifinelevel-1
  
  m_fmm%ifree=m_fmm%noct_used+1
  m_fmm%head(icoarselevel)=m_fmm%ifree
  hash_father(0)=icoarselevel

  check_refinement = (input%ilevel == ifinelevel)

  call open_cache(mdl, m_fmm, pack_size=storage_size(dummy_small_realdp)/32, &
       flush=pack_flush_build_fmm, combine=unpack_flush_build_fmm)

  ! Loop over fine grids
  do igrid=m%head(ifinelevel),m%tail(ifinelevel)

    hash_key(1:ndim)=m%grid(igrid)%ckey(1:ndim)
    hash_father(1:ndim)=hash_key(1:ndim)/2

    ! prob not worry
    in_domain = .true.
    do idim = 1, ndim
        in_domain = in_domain .and. hash_father(idim) .ge. m_fmm%box_ckey_min(idim,icoarselevel) &
            &                .and. hash_father(idim) .lt.  m_fmm%box_ckey_max(idim,icoarselevel)
    end do

    ! OPTIMIZATION NEEDED
    if (check_refinement) then
      flag_all_refined = all(m%grid(igrid)%refined)
    else
      flag_all_refined = .false.
    end if

    if(in_domain .and. (.not. flag_all_refined))then

        ! Access hash table
        ifather=hash_getp(m_fmm%grid_dict,hash_father)

        ! If grid does not exist, create it in memory
        if(ifather<=0)then

          ! Compute Cartesian keys of new oct
          cart_key(1:ndim)=int(hash_father(1:ndim),kind=4)

          ! Compute Hilbert keys of new octs
          ix(1:ndim)=cart_key(1:ndim)
          hk(1:nhilbert)=hilbert_key(ix,icoarselevel-1)

          ! Check if grid sits inside processor boundaries
          in_rank = ge_keys(hk,m_fmm%domain(icoarselevel)%b(1:nhilbert,mdl_self(mdl)-1)).and. &
                &    gt_keys(m_fmm%domain(icoarselevel)%b(1:nhilbert,mdl_self(mdl)),hk)

          if(in_rank)then

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

              ! Otherwise, determine parent processor and use the cache
              grid_cpu = m_fmm%domain(icoarselevel)%get_rank(hk)
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

    end if
  end do
  ! End loop over grids

  call close_cache(mdl)

  ! Multigrid oct statistics
  m_fmm%tail(icoarselevel)=m_fmm%ifree-1
  m_fmm%noct(icoarselevel)=m_fmm%tail(icoarselevel)-m_fmm%head(icoarselevel)+1
  m_fmm%noct_used=m_fmm%tail(icoarselevel)

  print *, "      <LEV>: ", icoarselevel, "| created: ", m_fmm%noct(icoarselevel)

  end associate

end subroutine build_fmm

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

#ifdef GRAV
  do idim=1,ndim
     do ind=1,twotondim
        mesh%f(ind,idim,igrid)=0.0d0
     end do
  end do
  do ind=1,twotondim
     mesh%phi(ind,igrid)=0.0d0
  end do
#endif

end subroutine unpack_flush_build_fmm
#endif
!################################################################
!################################################################
!################################################################
!################################################################
end module init_fmm_module
