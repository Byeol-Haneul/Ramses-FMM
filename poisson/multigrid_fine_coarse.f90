module multigrid_fine_coarse

  type :: level_count_t
     integer::ilevel,icount
  end type level_count_t

  type :: double_level_t
     integer::ilevel,ifine
  end type double_level_t

  type :: gs_step_t
     integer::ilevel,ifine
     logical::safe,redstep
  end type gs_step_t

#ifdef _CUDA
  use gpu_runner, only: gpu_restrict_mask, gpu_cmp_residual, gpu_gauss_seidel, gpu_restrict_residual, gpu_interpolate_correct
#endif

contains

#ifdef GRAV

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Mask restriction (bottom-up)
! ------------------------------------------------------------------------

recursive subroutine r_restrict_mask(pst,input,input_size,masked,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::masked,next_masked
  type(double_level_t)::input

  logical::allmasked
  integer::rID
  
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RESTRICT_MASK,pst%iUpper+1,input_size,output_size,input)
     call r_restrict_mask(pst%pLower,input,input_size,masked,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_masked)
     masked=masked*next_masked
  else
#ifdef _CUDA
     call gpu_restrict_mask(pst%s, input%ilevel, input%ifine, allmasked)
#else
     if(input%ifine==input%ilevel)then
        call restrict_mask(pst%s,pst%s%m,input%ifine,allmasked)
     else
        call restrict_mask(pst%s,pst%s%m_mg,input%ifine,allmasked)
     endif
#endif
     if(allmasked)then
        masked=1
     else
        masked=0
     endif
  endif

end subroutine r_restrict_mask

subroutine restrict_mask(s,m,ifinelevel,allmasked)
  use amr_parameters, only: nvector, nhilbert, ndim, twotondim
  use ramses_commons, only: ramses_t
  use amr_commons, only: mesh_t
  use nbors_utils
  use cache_commons
  use cache
  use hilbert
  use hash
  implicit none
  type(ramses_t)::s
  type(mesh_t)::m
  integer, intent(in) :: ifinelevel
  logical, intent(inout) :: allmasked

  integer(kind=8),dimension(0:ndim) :: hash_key
  integer :: ichild,ind,igrid,icell
  real(kind=8) :: ngpmask, mask_max
  real(kind=8) :: dtwotondim = (twotondim)
  type(msg_small_realdp)::dummy_small_realdp

  associate(r=>s%r,g=>s%g,m_mg=>s%m_mg,mdl=>s%mdl)

  ! Initialize volume fraction to zero at coarse level
  do igrid=m_mg%head(ifinelevel-1),m_mg%tail(ifinelevel-1)
     do ind=1,twotondim
        m_mg%f(ind,3,igrid)=0d0
     end do
  end do

  hash_key(0)=ifinelevel

  call open_cache(mdl, m_mg, pack_size=storage_size(dummy_small_realdp)/32, &
       init=init_flush_restrict_mask, flush=pack_flush_restrict_mask, combine=unpack_flush_restrict_mask)

  ! Loop over fine level grids
  do ichild=m%head(ifinelevel),m%tail(ifinelevel)

     ! Loop over cells
     do ind=1,twotondim

        hash_key(1:ndim)=m%grid(ichild)%ckey(1:ndim)

        ! Get parent cell using write-only cache
        call get_parent_cell(s,hash_key,igrid,icell,flush_cache=.true.,fetch_cache=.false.)

        ! Convert mask value to volume fraction
        ngpmask=(1d0+m%f(ind,3,ichild))/2d0/dtwotondim
        m_mg%f(icell,3,igrid)=m_mg%f(icell,3,igrid)+ngpmask

     end do
  end do

  call close_cache(mdl)

  ! Convert volume fraction back to to mask value for coarse level
  do igrid=m_mg%head(ifinelevel-1),m_mg%tail(ifinelevel-1)
     do ind=1,twotondim
        m_mg%f(ind,3,igrid)=2d0*m_mg%f(ind,3,igrid)-1d0
     end do
  end do

  ! Check mask state at coarse level
  mask_max=-1.0
  do igrid=m_mg%head(ifinelevel-1),m_mg%tail(ifinelevel-1)
     do ind=1,twotondim
        mask_max=MAX(mask_max,real(m_mg%f(ind,3,igrid),kind=8))
     end do
  end do
  allmasked=(mask_max<=0d0)

  end associate

end subroutine restrict_mask

subroutine init_flush_restrict_mask(mesh,igrid,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  type(mesh_t)::mesh
  integer::igrid
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)

#ifdef GRAV
  do ind=1,twotondim
     mesh%f(ind,3,igrid)=0.0d0
  end do
#endif

end subroutine init_flush_restrict_mask

subroutine pack_flush_restrict_mask(mesh,igrid,msg_size,msg_array)
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
     msg%realdp(ind)=mesh%f(ind,3,igrid)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_restrict_mask

subroutine unpack_flush_restrict_mask(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_small_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_small_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef GRAV
  do ind=1,twotondim
     mesh%f(ind,3,igrid)=mesh%f(ind,3,igrid)+msg%realdp(ind)
  end do
#endif

end subroutine unpack_flush_restrict_mask

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Residual computation
! ------------------------------------------------------------------------

recursive subroutine r_cmp_residual_mg(pst,input,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  use hash
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(double_level_t)::input

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CMP_RESIDUAL_MG,pst%iUpper+1,input_size,0,input)
     call r_cmp_residual_mg(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_cmp_residual(pst, input%ilevel, input%ifine)
#else
     if(input%ifine==input%ilevel)then
        call cmp_residual_mg(pst%s,pst%s%m,input%ifine)
     else
        call cmp_residual_mg(pst%s,pst%s%m_mg,input%ifine)
     endif
#endif
  endif

end subroutine r_cmp_residual_mg

subroutine cmp_residual_mg(s,m,ilevel)
  use amr_parameters, only: nvector, nhilbert, ndim, twondim, twotondim
  use ramses_commons, only: ramses_t
  use amr_commons, only: mesh_t
  use nbors_utils
  use cache_commons
  use cache
  use hilbert
  use hash
  use boundaries, only: init_bound_mg
  implicit none
  type(ramses_t)::s
  type(mesh_t)::m
  integer, intent(in) :: ilevel

  ! Computes the residual for MG levels, and stores it into grid(igrid)%f(ind,1)
    
  integer, dimension(1:3,1:2,1:8) :: iii, jjj
  real(kind=8),dimension(1:twotondim,0:twondim)::phi_nbor,dis_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer(kind=8),dimension(0:ndim) :: hash_nbor
  real(kind=8) :: dx, oneoverdx2, phi_c, dis_c, nb_sum
  integer  :: igrid, ind, inbor, idim, igridn, id, ig
  real(kind=8) :: dtwondim = (twondim)
  type(msg_twin_realdp)::dummy_twin_realdp

  associate(r=>s%r,g=>s%g,mdl=>s%mdl)

  ! Set constants
  dx = r%boxlen/2**ilevel
  oneoverdx2 = 1.0d0/(dx*dx)

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  call open_cache(mdl, m, pack_size=storage_size(dummy_twin_realdp)/32, &
       pack=pack_fetch_mg, unpack=unpack_fetch_mg, bound=init_bound_mg)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Get central oct potential and distance
     do ind=1,twotondim
        phi_nbor(ind,0)=m%phi(ind,igrid)
        dis_nbor(ind,0)=m%f(ind,3,igrid)
     end do

     ! Get neighboring octs potential
     do inbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=m%grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditions
        do idim=1,ndim
           if(r%periodic(idim))then
              if(hash_nbor(idim)< m%box_ckey_min(idim,ilevel))hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
              if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel))hash_nbor(idim)=m%box_ckey_min(idim,ilevel)
           endif
        enddo

        ! Get neighbouring grid using read-only cache
        call get_grid(s,hash_nbor,igridn,flush_cache=.false.,fetch_cache=.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=m%phi(ind,igridn)
              dis_nbor(ind,inbor)=m%f(ind,3,igridn)
           end do

        ! Otherwise set to zero and outside
        else
           do ind=1,twotondim
              phi_nbor(ind,inbor)=0.0
              dis_nbor(ind,inbor)=-1.0
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim

        ! Compute residual using 6 neighbors potential
        phi_c=m%phi(ind,igrid)
        dis_c=m%f(ind,3,igrid)

        nb_sum=0.0

        if(.not. btest(m%flag2(ind,igrid),0))then ! No scan needed

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 nb_sum=nb_sum+phi_nbor(id,ig)
              end do
           end do

        else ! Scan is required

           ! If cell is outside, set residual to zero
           if(dis_c<=0.0)then
              m%f(ind,1,igrid)=0.0
              cycle
           else

              ! Loop over neighbours
              do inbor=1,2
                 do idim=1,ndim
                    id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                    if(dis_nbor(id,ig)<=0.0)then
                       nb_sum=nb_sum+phi_c/dis_c*dis_nbor(id,ig)
                    else
                       nb_sum=nb_sum+phi_nbor(id,ig)
                    endif
                 end do
              end do

           endif

        endif

        ! Store residual in f(ind,1)
        m%f(ind,1,igrid)=-oneoverdx2*( nb_sum - dtwondim*phi_c )+m%f(ind,2,igrid)

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache(mdl)

  end associate

end subroutine cmp_residual_mg

subroutine pack_fetch_mg(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_twin_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_twin_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp_phi(ind)=mesh%phi(ind,igrid)
     msg%realdp_dis(ind)=mesh%f(ind,3,igrid)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_mg

subroutine unpack_fetch_mg(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_twin_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_twin_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef GRAV
  do ind=1,twotondim
     mesh%phi(ind,igrid)=msg%realdp_phi(ind)
     mesh%f(ind,3,igrid)=msg%realdp_dis(ind)
  end do
#endif

end subroutine unpack_fetch_mg

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Gauss-Seidel Red-Black sweeps
! ------------------------------------------------------------------------

recursive subroutine r_gauss_seidel_mg(pst,input,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  use hash
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(gs_step_t)::input

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_GAUSS_SEIDEL_MG,pst%iUpper+1,input_size,0,input)
     call r_gauss_seidel_mg(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_gauss_seidel(pst%s, input%ilevel, input%ifine, input%safe, input%redstep)
#else
     if(input%ifine==input%ilevel)then
        call gauss_seidel_mg(pst%s,pst%s%m,input%ifine,input%safe,input%redstep)
     else
        call gauss_seidel_mg(pst%s,pst%s%m_mg,input%ifine,input%safe,input%redstep)
     endif
#endif
  endif

end subroutine r_gauss_seidel_mg

subroutine gauss_seidel_mg(s,m,ilevel,safe,redstep)
  use amr_parameters, only: nvector, nhilbert, ndim, twondim, twotondim
  use ramses_commons, only: ramses_t
  use amr_commons, only: mesh_t
  use nbors_utils
  use cache_commons
  use cache
  use hilbert
  use hash
  use boundaries, only: init_bound_mg
  implicit none
  type(ramses_t)::s
  type(mesh_t)::m
  integer, intent(in) :: ilevel
  logical, intent(in) :: safe
  logical, intent(in) :: redstep

  ! Perform a Gauss-Seidel update of grid(igrid)%phi(ind).
  ! The domain mask is also needed.

  integer, dimension(1:3,1:2,1:8) :: iii, jjj
  real(kind=8),dimension(1:twotondim,0:twondim)::phi_nbor,dis_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer(kind=8),dimension(0:ndim) :: hash_nbor
  real(kind=8) :: phi_c, dis_c, dx2, nb_sum, weight
  integer  :: igrid, ind, inbor, idim, igridn, id, ig, ind0
  real(kind=8) :: dtwondim = (twondim)
  type(msg_twin_realdp)::dummy_twin_realdp

  integer, dimension(1:4) :: ired, iblack

  associate(r=>s%r,g=>s%g,mdl=>s%mdl)

  ! Set constants
  dx2  = (r%boxlen/2**ilevel)**2

  ired  (1:4)=(/1,4,6,7/)
  iblack(1:4)=(/2,3,5,8/)

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  call open_cache(mdl, m, pack_size=storage_size(dummy_twin_realdp)/32, &
       pack=pack_fetch_mg, unpack=unpack_fetch_mg, bound=init_bound_mg)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Get central oct potential and distance
     do ind=1,twotondim
        phi_nbor(ind,0)=m%phi(ind,igrid)
        dis_nbor(ind,0)=m%f(ind,3,igrid)
     end do

     ! Get neighboring octs potential
     do inbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=m%grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditions
        do idim=1,ndim
           if(r%periodic(idim))then
              if(hash_nbor(idim)< m%box_ckey_min(idim,ilevel))hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
              if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel))hash_nbor(idim)=m%box_ckey_min(idim,ilevel)
           endif
        enddo

        ! Get neighbouring grid using read-only cache
        call get_grid(s,hash_nbor,igridn,flush_cache=.false.,fetch_cache=.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=m%phi(ind,igridn)
              dis_nbor(ind,inbor)=m%f(ind,3,igridn)
           end do

        ! Otherwise set to zero and outside
        else
           do ind=1,twotondim
              phi_nbor(ind,inbor)=0.0
              dis_nbor(ind,inbor)=-1.0
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells, with red/black ordering
     do ind0=1,twotondim/2      ! Only half of the cells for a red or black sweep

        if(redstep) then
           ind = ired  (ind0)
        else
           ind = iblack(ind0)
        end if

        ! Compute residual using 6 neighbors potential
        phi_c=m%phi(ind,igrid)
        dis_c=m%f(ind,3,igrid)

        nb_sum=0.0

        if(.not. btest(m%flag2(ind,igrid),0))then ! No scan needed

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 nb_sum=nb_sum+phi_nbor(id,ig)
              end do
           end do

           ! Update the potential, solving for potential on icell_amr
           m%phi(ind,igrid)=(nb_sum-dx2*m%f(ind,2,igrid))/dtwondim

        else ! Scan is required

           ! If cell is outside, don't update phi
           if(dis_c<=0.0)cycle
           if(safe .and. dis_c<1.0)cycle

           weight=0.0d0   ! Central weight for "Solve G-S"

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 if(dis_nbor(id,ig)<=0.0)then
                    weight=weight+dis_nbor(id,ig)/dis_c
                 else
                    nb_sum=nb_sum+phi_nbor(id,ig)
                 endif
              end do
           end do

           ! Update the potential
           m%phi(ind,igrid)=(nb_sum-dx2*m%f(ind,2,igrid))/(dtwondim - weight)

        endif

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache(mdl)

  end associate

end subroutine gauss_seidel_mg

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Reset correction
! ------------------------------------------------------------------------

recursive subroutine r_reset_correction(pst,ilevel,input_size)
  use mdl_module
  use amr_parameters, only: twotondim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::igrid
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RESET_CORRECTION,pst%iUpper+1,input_size,0,ilevel)
     call r_reset_correction(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_reset_corr(pst%s, ilevel)
#else
     do igrid=pst%s%m_mg%head(ilevel),pst%s%m_mg%tail(ilevel)
        pst%s%m_mg%phi(1:twotondim,igrid)=0.0d0
     end do
#endif
  endif

end subroutine r_reset_correction

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Residual restriction
! ------------------------------------------------------------------------

recursive subroutine r_restrict_residual(pst,input,input_size)
  use mdl_module
  use amr_parameters, only: twotondim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(double_level_t)::input

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RESTRICT_RESIDUAL,pst%iUpper+1,input_size,0,input)
     call r_restrict_residual(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_restrict_residual(pst%s, input%ilevel, input%ifine)
#else
     if(input%ifine==input%ilevel)then
        call restrict_residual(pst%s,pst%s%m,input%ifine)
     else
        call restrict_residual(pst%s,pst%s%m_mg,input%ifine)
     endif
#endif
  endif

end subroutine r_restrict_residual

subroutine restrict_residual(s,m,ifinelevel)
  use amr_parameters, only: nvector, nhilbert, ndim, twondim, twotondim
  use ramses_commons, only: ramses_t
  use amr_commons, only: mesh_t
  use nbors_utils
  use cache_commons
  use cache
  use hilbert
  use hash
  implicit none
  type(ramses_t)::s
  type(mesh_t)::m
  integer, intent(in) :: ifinelevel

  ! Restrict the residual of the fine level (stored in m%f(ind,1,igrid))
  ! into the rhs of the coarse level (stored in m%f(icell,2,igrid))
  ! For interior coarse cell only (we need the mask stored in m%f(icell,3,igrid))

  integer :: ichild, ind
  integer :: igrid, icell
  real(kind=8) :: dtwotondim = (twotondim)
  integer(kind=8),dimension(0:ndim) :: hash_key
  type(msg_small_realdp)::dummy_small_realdp

  associate(r=>s%r,g=>s%g,m_mg=>s%m_mg,mdl=>s%mdl)

  ! Set rhs to zero in coarse cells
  do igrid=m_mg%head(ifinelevel-1),m_mg%tail(ifinelevel-1)
     do ind=1,twotondim
        m_mg%f(ind,2,igrid)=0.0
     end do
  end do

  hash_key(0)=ifinelevel

  call open_cache(mdl, m_mg, pack_size=storage_size(dummy_small_realdp)/32, &
       pack=pack_fetch_restrict_res, unpack=unpack_fetch_restrict_res, &
       init=init_flush_restrict_res, flush=pack_flush_restrict_res, combine=unpack_flush_restrict_res)

  ! Loop over grids
  do ichild=m%head(ifinelevel),m%tail(ifinelevel)

     ! Loop over cells
     do ind=1,twotondim        

        ! Is fine cell masked?
        if(m%f(ind,3,ichild)<=0d0)cycle

        hash_key(1:ndim)=m%grid(ichild)%ckey(1:ndim)

        ! Get parent cell using read-write cache
        call get_parent_cell(s,hash_key,igrid,icell,flush_cache=.true.,fetch_cache=.true.)

        ! Is coarse cell masked?
        if(m_mg%f(icell,3,igrid)<=0d0)cycle

        ! Stack fine cell residual in coarse cell rhs
        m_mg%f(icell,2,igrid)=m_mg%f(icell,2,igrid)+m%f(ind,1,ichild)/dtwotondim

     end do
  end do

  call close_cache(mdl)

  end associate

end subroutine restrict_residual

subroutine pack_fetch_restrict_res(mesh,igrid,msg_size,msg_array)
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
     msg%realdp(ind)=mesh%f(ind,3,igrid)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_restrict_res

subroutine unpack_fetch_restrict_res(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_small_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_small_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef GRAV
  do ind=1,twotondim
     mesh%f(ind,3,igrid)=msg%realdp(ind)
  end do
#endif

end subroutine unpack_fetch_restrict_res

subroutine init_flush_restrict_res(mesh,igrid,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  type(mesh_t)::mesh
  integer::igrid
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)

#ifdef GRAV
  do ind=1,twotondim
     mesh%f(ind,2,igrid)=0.0d0
  end do
#endif

end subroutine init_flush_restrict_res

subroutine pack_flush_restrict_res(mesh,igrid,msg_size,msg_array)
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
     msg%realdp(ind)=mesh%f(ind,2,igrid)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_restrict_res

subroutine unpack_flush_restrict_res(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_small_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_small_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef GRAV
  do ind=1,twotondim
     mesh%f(ind,2,igrid)=mesh%f(ind,2,igrid)+msg%realdp(ind)
  end do
#endif

end subroutine unpack_flush_restrict_res

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Interpolation and correction
! ------------------------------------------------------------------------

recursive subroutine r_interpolate_and_correct(pst,input,input_size)
  use mdl_module
  use amr_parameters, only: twotondim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(double_level_t)::input

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INTERPOLATE_AND_CORRECT,pst%iUpper+1,input_size,0,input)
     call r_interpolate_and_correct(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_interpolate_correct(pst%s, input%ilevel, input%ifine)
#else
     if(input%ifine==input%ilevel)then
        call interpolate_and_correct(pst%s,pst%s%m,input%ifine)
     else
        call interpolate_and_correct(pst%s,pst%s%m_mg,input%ifine)
     end if
#endif
  endif

end subroutine r_interpolate_and_correct

subroutine interpolate_and_correct(s,m,ifinelevel)
  use amr_parameters, only: nvector, nhilbert, ndim, twondim, twotondim, threetondim
  use ramses_commons, only: ramses_t
  use amr_commons, only: mesh_t
  use nbors_utils
  use cache_commons
  use cache
  use hilbert
  use hash
  use boundaries, only: init_bound_mg
  implicit none
  type(ramses_t)::s
  type(mesh_t)::m
  integer, intent(in) :: ifinelevel

  ! Interpolate the solution of the coarse level (stored in m%phi(icell,igrid))
  ! and correct the solution of the fine level (stored in m%phi(ind,igrid))

  integer(kind=8),dimension(0:ndim) :: hash_key
  integer,dimension(1:threetondim) :: igrid_nbor,ind_nbor
  integer  :: ichild, ind
  real(kind=8) :: aa, bb, cc, dd, coeff
  real(kind=8), dimension(1:8)     :: bbb
  integer,  dimension(1:8,1:8) :: ccc
  integer::ind_average,ind_father
  integer::igrid_nbr,ind_nbr
  real(kind=8),dimension(1:twotondim)::corr
  type(msg_small_realdp)::dummy_small_realdp

  associate(r=>s%r,g=>s%g,m_mg=>s%m_mg,mdl=>s%mdl)

  ! Local constants
  aa = 1.0D0/4.0D0**ndim
  bb = 3.0D0*aa
  cc = 9.0D0*aa
  dd = 27.D0*aa 
  bbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  call open_cache(mdl, m_mg, pack_size=storage_size(dummy_small_realdp)/32, &
       pack=pack_fetch_phi, unpack=unpack_fetch_phi, bound=init_bound_mg)

  hash_key(0)=ifinelevel

  ! Loop over fine level grids
  do ichild=m%head(ifinelevel),m%tail(ifinelevel)

     ! For fine level, correction is interpolated from coarser level solution
     hash_key(1:ndim)=m%grid(ichild)%ckey(1:ndim)

     ! Get 3**ndim neighbouring parent cell using a read-only cache
     call get_threetondim_nbor_parent_cell(s,hash_key,igrid_nbor,ind_nbor,flush_cache=.false.,fetch_cache=.true.)

     ! Loop over cells
     do ind=1,twotondim

        ! Set correction to zero
        corr(ind)=0d0

        ! Fine cell is masked as "outside": no correction
        if(m%f(ind,3,ichild)<=0.0)cycle

        ! Loop over relevant parent cells
        do ind_average=1,twotondim
           ind_father=ccc(ind_average,ind)
           coeff=bbb(ind_average)
           igrid_nbr=igrid_nbor(ind_father)
           ind_nbr=ind_nbor(ind_father)
           if (igrid_nbr>0) then
              corr(ind)=corr(ind)+coeff*m_mg%phi(ind_nbr,igrid_nbr)
           endif
        end do

     end do
     ! End loop over cells

     do ind=1,threetondim
        call unlock_cache(m_mg,igrid_nbor(ind))
     end do

     ! Add correction to fine level solution
     do ind=1,twotondim
        m%phi(ind,ichild)=m%phi(ind,ichild)+corr(ind)
     end do

  end do
  ! End loop over grids

  call close_cache(mdl)

  end associate

end subroutine interpolate_and_correct

subroutine pack_fetch_phi(mesh,igrid,msg_size,msg_array)
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
     msg%realdp(ind)=mesh%phi(ind,igrid)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_phi

subroutine unpack_fetch_phi(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_small_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_small_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef GRAV
  do ind=1,twotondim
     mesh%phi(ind,igrid)=msg%realdp(ind)
  end do
#endif

end subroutine unpack_fetch_phi

subroutine pack_fetch_rho(mesh,igrid,msg_size,msg_array)
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
     msg%realdp(ind)=mesh%rho(ind,igrid)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_rho

subroutine unpack_fetch_rho(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_small_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_small_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef GRAV
  do ind=1,twotondim
     mesh%rho(ind,igrid)=msg%realdp(ind)
  end do
#endif

end subroutine unpack_fetch_rho

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Flag settings used to speed-up the sweeps
! ------------------------------------------------------------------------

recursive subroutine r_set_scan_flag(pst,input,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  use hash
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(double_level_t)::input

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SET_SCAN_FLAG,pst%iUpper+1,input_size,0,input)
     call r_set_scan_flag(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     ! Do nothing
#else
     if(input%ifine==input%ilevel)then
        call set_scan_flag(pst%s,pst%s%m,input%ifine)
     else
        call set_scan_flag(pst%s,pst%s%m_mg,input%ifine)
     endif
#endif
  endif

end subroutine r_set_scan_flag

subroutine set_scan_flag(s,m,ilevel)
  use amr_parameters, only: nvector, nhilbert, ndim, twondim, twotondim, threetondim
  use ramses_commons, only: ramses_t
  use amr_commons, only: mesh_t
  use nbors_utils
  use cache_commons
  use cache
  use hilbert
  use hash
  use boundaries, only: init_bound_mg
  implicit none
  type(ramses_t)::s
  type(mesh_t)::m
  integer, intent(in) :: ilevel
  !
  integer :: ind, igrid, igridn, inbor, idim, id, ig
  integer, dimension(1:3,1:2,1:8)::iii, jjj
  real(kind=8),dimension(1:twotondim,0:twondim)::dis_nbor
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  real(kind=8)::dis_c
  type(msg_small_realdp)::dummy_small_realdp

  associate(r=>s%r,g=>s%g,mdl=>s%mdl)

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  call open_cache(mdl, m, pack_size=storage_size(dummy_small_realdp)/32, &
       pack=pack_fetch_scan, unpack=unpack_fetch_scan, bound=init_bound_mg)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Get central oct potential
     do ind=1,twotondim
        dis_nbor(ind,0)=m%f(ind,3,igrid)
     end do

     ! Get neighboring octs potential
     do inbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=m%grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditions
        do idim=1,ndim
           if(r%periodic(idim))then
              if(hash_nbor(idim)< m%box_ckey_min(idim,ilevel))hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
              if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel))hash_nbor(idim)=m%box_ckey_min(idim,ilevel)
           endif
        enddo

        ! Get neighbouring grid using read-only cache
        call get_grid(s,hash_nbor,igridn,flush_cache=.false.,fetch_cache=.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              dis_nbor(ind,inbor)=m%f(ind,3,igridn)
           end do

        ! Otherwise set to "outside"
        else
           do ind=1,twotondim
              dis_nbor(ind,inbor)=-1.0
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim

        ! Compute residual using 6 neighbors potential
        dis_c=m%f(ind,3,igrid)

        ! If cell is entirely inside, set flag tentatively to 0 (no scan needed)
        if(dis_c==1.0)then
           m%flag2(ind,igrid)=0

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 ! If one neighbour is outside, then scan needed
                 if(dis_nbor(id,ig)<=0.0)then
                    m%flag2(ind,igrid)=max(m%flag2(ind,igrid),1)
                 endif
              end do
           end do

         ! If cell is even partially outside, then scan needed
        else
           m%flag2(ind,igrid)=1           
        endif

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache(mdl)

  end associate

end subroutine set_scan_flag

subroutine pack_fetch_scan(mesh,igrid,msg_size,msg_array)
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
     msg%realdp(ind)=mesh%f(ind,3,igrid)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_scan

subroutine unpack_fetch_scan(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_small_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_small_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef GRAV
  do ind=1,twotondim
     mesh%f(ind,3,igrid)=msg%realdp(ind)
  end do
#endif

end subroutine unpack_fetch_scan

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Compute norm of residual 
! ------------------------------------------------------------------------

recursive subroutine r_cmp_residual_norm2(pst,ilevel,input_size,norm2,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  real(kind=8)::norm2,next_norm2

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CMP_RESIDUAL_NORM2,pst%iUpper+1,input_size,output_size,ilevel)
     call r_cmp_residual_norm2(pst%pLower,ilevel,input_size,norm2,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_norm2)
     norm2=norm2+next_norm2
  else
#ifdef _CUDA
     call gpu_residual_norm(pst%s, ilevel, norm2)
#else
     call cmp_residual_norm2(pst%s%r,pst%s%m,ilevel,norm2)
#endif
  endif

end subroutine r_cmp_residual_norm2

subroutine cmp_residual_norm2(r,m,ilevel, norm2)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: run_t, mesh_t
  implicit none
  type(run_t)::r
  type(mesh_t)::m
  integer,  intent(in)  :: ilevel
  real(kind=8), intent(out) :: norm2

  real(kind=8) :: dx2
  integer  :: ind, igrid

  ! Set constants
  dx2  = (r%boxlen/2**ilevel)**2
  norm2 = 0.0d0

  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Loop over cells
     do ind=1,twotondim
        if(m%f(ind,3,igrid)<=0.0)cycle      ! Do not count masked cells
        norm2 = norm2 + m%f(ind,1,igrid)**2
     end do
     ! End loop over cells

  end do
  ! End loop over grids

  norm2 = dx2*norm2

end subroutine cmp_residual_norm2

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

#endif

end module multigrid_fine_coarse
