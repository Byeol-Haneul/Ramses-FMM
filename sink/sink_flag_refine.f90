module sink_flag_module
contains
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine sink_flag(s,p,ilevel)
  use amr_parameters, only: ndim,twotondim,dp
  use amr_commons, only: nbor,oct
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use nbors_utils
  use cache_commons
  use cache
  use marshal, only: pack_fetch_flag, unpack_fetch_flag
  use flag_utils, only: init_flush_initflag, pack_flush_initflag, unpack_flush_initflag
  use hilbert
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  !==================================================================
  ! This routine flag for refinement cells that are close enough
  ! from sink particles.
  !==================================================================
  real(dp),dimension(:,:),allocatable::x_nei
  integer,dimension(1:ndim)::ckey,ckey_ref,ckey_nbor
  integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
  integer::i,j,k,ipart,icellp,icelln,ind,idim,i_nei,n_nei,nrad
  integer,dimension(1:ndim)::ix
  real(dp)::dx_loc,vol_loc,x,y,z,rrad,rr
  real(dp),dimension(1:3)::xcen,xnei
  type(oct),pointer::gridp,gridn
  type(msg_int4)::dummy_int4

#ifdef HYDRO
#if NDIM==3
  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Mesh spacing in that level
  dx_loc = r%boxlen / 2**ilevel 
  vol_loc = dx_loc**ndim

  ! Compute number of cells within sink sphere
  rrad = r%sink_accretion_radius/dx_loc
  nrad = ceiling(rrad)
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

  ! Loop over sink particles at current level
  hash_nbor(0) = ilevel+1
  do ipart = p%headp(ilevel), p%tailp(ilevel)

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
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module sink_flag_module
