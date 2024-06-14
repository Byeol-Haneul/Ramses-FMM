module rho_fine_module
contains
!###############################################
!###############################################
!###############################################
!###############################################
#ifdef GRAV
subroutine m_rho_fine(pst,ilevel,rtype)
  use amr_parameters, only: dp,ndim
  use ramses_commons, only: pst_t
  use amr_commons, only: multipole_t
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !------------------------------------------------------------------
  ! This master routine computes the mass density field to be used
  ! as source term in the Poisson solver.
  ! The density field is computed for all levels greater than ilevel.
  ! On output, particles are sorted according to their grid level of
  ! refinement, and inside their level, they are sorted according to
  ! their grid Hilbert order.
  !------------------------------------------------------------------
  type(multipole_t)::multipole_tot
  integer::i,input_size,rtype ! rtype 1 all 2 dm 3 star 4 gas
  integer,allocatable,dimension(:)::input_array
  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)

  if(.not. r%poisson)return
  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,'(" Entering rho_fine for level ",I2)')ilevel

  !---------------------------
  ! Reset multipole to zero
  !---------------------------
  if(ilevel==r%levelmin)then
     multipole_tot%q=0d0
     input_size=storage_size(multipole_tot)/32
     call r_broadcast_multipole(pst,multipole_tot,input_size)
  endif
  
  !-------------------------------------------------------
  ! Initialize rho to analytical and baryon density field
  !-------------------------------------------------------
  ! Loop over all finer levels from fine to coarse
  do i=r%nlevelmax,ilevel,-1

     ! Compute gas multipole expansion
     if(r%hydro)then

        ! Set multipoles in all leaf cells
        if(m%noct_tot(i)>0)then
           if(r%verbose)write(*,'(" Compute leaf multipoles for level ",I2)')i
           call r_multipole_leaf_cells(pst,i,1)
        endif

        ! Average down multipoles in all split cells
        if(i<r%nlevelmax)then
           if(m%noct_tot(i+1)>0)then
              if(r%verbose)write(*,'(" Compute split multipoles for level ",I2)')i
              call r_multipole_split_cells(pst,i,1)
           endif
        endif

     endif

     ! Reset array rho to zero
     if(m%noct_tot(i)>0)then
        call r_reset_rho(pst,i,1)
     endif

     ! Gas mass deposition using pseudo-particles
     if(r%hydro.AND.m%noct_tot(i)>0.AND.(rtype==1 .or. rtype==4))then
        if(r%verbose)write(*,'(" Compute rho from multipoles for level ",I2)')i
        call r_cic_multipole(pst,i,1)
     endif

  end do
  ! End loop over finer levels

  !-------------------------------------------------------
  ! Compute particle contribution to density field
  !-------------------------------------------------------
  if(r%pic.AND.(rtype.ne.4))then
     do i=ilevel,r%nlevelmax
        if(m%noct_tot(i)>0)then
           if(r%verbose)write(*,'(" Compute rho from particles for level ",I2)')i
           allocate(input_array(1:2))
           input_array(1)=i
           input_array(2)=rtype
           call r_cic_part(pst,input_array,2)
           deallocate(input_array)
        endif
        if(m%noct_tot(i)>0.AND.i<r%nlevelmax)then
           if(r%verbose)write(*,'(" Split particles for level ",I2)')i
           allocate(input_array(1:2))
           input_array(1)=i
           input_array(2)=rtype
           call r_split_part(pst,input_array,2)
           deallocate(input_array)
        endif
     end do
  endif

  !---------------------------------------------------------------------
  ! Collect multipole contribution from all CPU and broadcast rho_tot
  !---------------------------------------------------------------------
  if(ilevel==r%levelmin)then

     ! Collect local multipole from all CPU
     call r_collect_multipole(pst,ilevel,1,multipole_tot,storage_size(multipole_tot)/32)

     ! Broadcast total multipole to all CPU
     call r_broadcast_multipole(pst,multipole_tot,storage_size(multipole_tot)/32)

     if(r%verbose)write(*,*)'rho_average=',g%rho_tot
  endif  

  end associate

end subroutine m_rho_fine
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_multipole_leaf_cells(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_MULTIPOLE_LEAF_CELLS,pst%iUpper+1,input_size,0,ilevel)
     call r_multipole_leaf_cells(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call multipole_leaf_cells(pst%s%r,pst%s%g,pst%s%m,ilevel)
  endif

end subroutine r_multipole_leaf_cells
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine multipole_leaf_cells(r,g,m,ilevel)
  use amr_parameters, only: ndim,dp,twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  use cache_commons
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !-------------------------------------------------------------------
  ! This routine compute the monopole and dipole of the gas mass and
  ! the analytical profile (if any) within each cell.
  ! For pure particle runs, this is not necessary and the
  ! routine is not even called.
  !-------------------------------------------------------------------
  integer::igrid,ind,idim,ivar,nstride,icell
  real(dp),dimension(1:ndim)::xx
  real(kind=8)::dx_loc,vol_loc,mmm,dd
  logical::leaf_cell

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

#ifdef HYDRO
  ! Initialize multipole fields to zero
  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim
        do idim=1,ndim+1
           m%grid(igrid)%unew(ind,idim)=0.0D0
        end do
     end do
  end do
#endif

  !-------------------------------------------------------
  ! Compute contribution of leaf cells to mass multipoles
  !-------------------------------------------------------
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim

        leaf_cell=m%grid(igrid)%refined(ind).EQV..FALSE.

        ! For leaf cells only
        if(leaf_cell)then

           ! Cell coordinates
           do idim=1,ndim
              nstride=2**(idim-1)
              xx(idim)=(2*m%grid(igrid)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx_loc-m%skip(idim)
           end do
#ifdef HYDRO
           ! Add gas mass
           mmm=max(m%grid(igrid)%uold(ind,1),r%smallr)*vol_loc
           m%grid(igrid)%unew(ind,1)=m%grid(igrid)%unew(ind,1)+mmm
           do idim=1,ndim
              m%grid(igrid)%unew(ind,idim+1)=m%grid(igrid)%unew(ind,idim+1)+mmm*xx(idim)
           end do
#endif
           ! Add analytical density profile
           if(r%gravity_type < 0)then
              call rho_ana(xx,dd,dx_loc,r%gravity_params)
              mmm=max(dd,r%smallr)*vol_loc
#ifdef HYDRO
              m%grid(igrid)%unew(ind,1)=m%grid(igrid)%unew(ind,1)+mmm
              do idim=1,ndim
                 m%grid(igrid)%unew(ind,idim+1)=m%grid(igrid)%unew(ind,idim+1)+mmm*xx(idim)
              end do
#endif
           end if
        endif
     end do
     ! End loop over cells
  end do
  ! End loop over grids

end subroutine multipole_leaf_cells
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_multipole_split_cells(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_MULTIPOLE_SPLIT_CELLS,pst%iUpper+1,input_size,0,ilevel)
     call r_multipole_split_cells(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call multipole_split_cells(pst%s,ilevel)
  endif

end subroutine r_multipole_split_cells
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine multipole_split_cells(s,ilevel)
  use amr_parameters, only: ndim,dp,twotondim
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use hydro_flag_module, only: pack_fetch_hydro,unpack_fetch_hydro
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !-------------------------------------------------------------------
  ! This routine compute the monopole and dipole of the gas mass and
  ! the analytical profile (if any) within each cell.
  ! For pure particle runs, this is not necessary and the
  ! routine is not even called.
  !-------------------------------------------------------------------
  integer::ind,idim,ivar,ioct,icell
  real(kind=8)::average
  integer(kind=8),dimension(0:ndim)::hash_key
  logical::leaf_cell
  type(oct),pointer::gridp
  type(msg_realdp)::dummy_realdp

  associate(r=>s%r,g=>s%g,m=>s%m)
  
  !-------------------------------------------------------
  ! Perform octree restriction from level ilevel+1
  !-------------------------------------------------------
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_realdp)/32,&
                     pack=pack_fetch_hydro,unpack=unpack_fetch_hydro,&
                     init=init_flush_multipole, flush=pack_flush_multipole, combine=unpack_flush_multipole)

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ioct=m%head(ilevel+1),m%tail(ilevel+1)
     hash_key(1:ndim)=m%grid(ioct)%ckey(1:ndim)
     ! Get parent cell using a write-only cache
     call get_parent_cell(s,hash_key,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.false.)
#ifdef HYDRO
     ! Average conservative variables
     do ivar=1,ndim+1
        average=0.0d0
        do ind=1,twotondim
           average=average+m%grid(ioct)%unew(ind,ivar)
        end do
        ! Scatter result to cell
        gridp%unew(icell,ivar)=average
     end do
#endif
  end do

  call close_cache(s,m%grid_dict)

  end associate

end subroutine multipole_split_cells
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_multipole(grid,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  type(oct)::grid
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  
  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
#ifdef HYDRO
  do ivar=1,ndim+1
     do ind=1,twotondim
        grid%unew(ind,ivar)=0.0
     end do
  end do
#endif
  
end subroutine init_flush_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_multipole(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_realdp)::msg

#ifdef HYDRO
  do ivar=1,ndim+1
     do ind=1,twotondim
        msg%realdp(ind,ivar)=grid%unew(ind,ivar)
     end do
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_multipole(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  type(msg_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)
  
#ifdef HYDRO
  do ivar=1,ndim+1
     do ind=1,twotondim
        if(grid%refined(ind))then
           grid%unew(ind,ivar)=grid%unew(ind,ivar)+msg%realdp(ind,ivar)
        endif
     end do
  end do
#endif

end subroutine unpack_flush_multipole
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_reset_rho(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RESET_RHO,pst%iUpper+1,input_size,0,ilevel)
     call r_reset_rho(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call reset_rho(pst%s%r,pst%s%g,pst%s%m,ilevel)
  endif

end subroutine r_reset_rho
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine reset_rho(r,g,m,ilevel)
  use amr_parameters, only: twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !-------------------------------------------------------------------
  ! This routine compute array rho (source term for Poisson equation)
  ! by first reseting array rho to zero, then 
  ! by depositing the gas multipole mass in each cells using CIC.
  ! For pure particle runs, the gas mass deposition is not done
  ! and the routine only set rho to zero.
  !-------------------------------------------------------------------
  integer::igrid,ind

#ifdef GRAV
  ! Initialize density field to zero
  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim
        m%grid(igrid)%rho(ind)=0.0D0
        m%grid(igrid)%nref(ind)=0.0D0
     end do
  end do
#endif  

end subroutine reset_rho
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_cic_multipole(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CIC_MULTIPOLE,pst%iUpper+1,input_size,0,ilevel)
     call r_cic_multipole(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call cic_multipole(pst%s,ilevel)
  endif

end subroutine r_cic_multipole
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cic_multipole(s,ilevel)
  use mdl_module
  use amr_parameters, only: ndim,twotondim,dp
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use multigrid_fine_coarse, only:pack_fetch_phi,unpack_fetch_phi
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !
  ! Local variables
  real(dp),dimension(1:ndim)::x,dd,dg
  integer,dimension(1:ndim)::ig,id
  real(dp),dimension(1:twotondim)::vol
  integer,dimension(1:ndim,1:twotondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::inbor,igrid,ind,idim
  integer::icell
  real(kind=8)::dx_loc,vol_loc,mmm,mask
  type(oct),pointer::gridp
  type(msg_twin_realdp)::dummy_twin_realdp

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)
    
  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Use hash table directly for cells (not for grids)
  hash_nbor(0)=ilevel+1

  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_twin_realdp)/32,&
                pack=pack_fetch_phi,unpack=unpack_fetch_phi,&
                init=init_flush_rho, flush=pack_flush_rho, combine=unpack_flush_rho)

  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Loop over cells
     do ind=1,twotondim

#ifdef HYDRO        
        ! Compute pseudo particle mass
        mmm=m%grid(igrid)%unew(ind,1)

        ! Compute pseudo particle (centre of mass) position
        if(mmm==0)then
           write(*,*)'Sorry divide by zero'
           write(*,*)m%grid(igrid)%unew(ind,1:nvar)
           write(*,*)m%grid(igrid)%uold(ind,1:nvar)
           write(*,*)m%grid(igrid)%refined(ind)
           call mdl_abort(mdl)
        endif
        x(1:ndim)=m%grid(igrid)%unew(ind,2:ndim+1)/mmm
        
        ! Compute total multipole
        if(ilevel==r%levelmin)then
           do idim=1,ndim+1
              g%multipole%q(idim)=g%multipole%q(idim)+m%grid(igrid)%unew(ind,idim)
           end do
        endif
#endif
        ! Rescale particle position at level ilevel
        do idim=1,ndim
           x(idim)=x(idim)/dx_loc
        end do
     
        ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
        do idim=1,ndim
           dd(idim)=x(idim)+0.5D0
           id(idim)=int(dd(idim))
           dd(idim)=dd(idim)-id(idim)
           dg(idim)=1.0D0-dd(idim)
           ig(idim)=id(idim)-1
        end do

        ! Periodic boundary conditions
        do idim=1,ndim
           if(ig(idim)<0)ig(idim)=m%ckey_max(ilevel+1)-1
           if(id(idim)==m%ckey_max(ilevel+1))id(idim)=0
        enddo

        ! Compute cloud volumes
#if NDIM==1
        vol(1)=dg(1)
        vol(2)=dd(1)
#endif
#if NDIM==2
        vol(1)=dg(1)*dg(2)
        vol(2)=dd(1)*dg(2)
        vol(3)=dg(1)*dd(2)
        vol(4)=dd(1)*dd(2)
#endif
#if NDIM==3
        vol(1)=dg(1)*dg(2)*dg(3)
        vol(2)=dd(1)*dg(2)*dg(3)
        vol(3)=dg(1)*dd(2)*dg(3)
        vol(4)=dd(1)*dd(2)*dg(3)
        vol(5)=dg(1)*dg(2)*dd(3)
        vol(6)=dd(1)*dg(2)*dd(3)
        vol(7)=dg(1)*dd(2)*dd(3)
        vol(8)=dd(1)*dd(2)*dd(3)
#endif

        ! Compute cells Cartesian key
#if NDIM==1
        ckey(1,1)=ig(1)
        ckey(1,2)=id(1)
#endif
#if NDIM==2
        ckey(1:2,1)=(/ig(1),ig(2)/)
        ckey(1:2,2)=(/id(1),ig(2)/)
        ckey(1:2,3)=(/ig(1),id(2)/)
        ckey(1:2,4)=(/id(1),id(2)/)
#endif
#if NDIM==3
        ckey(1:3,1)=(/ig(1),ig(2),ig(3)/)
        ckey(1:3,2)=(/id(1),ig(2),ig(3)/)
        ckey(1:3,3)=(/ig(1),id(2),ig(3)/)
        ckey(1:3,4)=(/id(1),id(2),ig(3)/)
        ckey(1:3,5)=(/ig(1),ig(2),id(3)/)
        ckey(1:3,6)=(/id(1),ig(2),id(3)/)
        ckey(1:3,7)=(/ig(1),id(2),id(3)/)
        ckey(1:3,8)=(/id(1),id(2),id(3)/)
#endif     

#ifdef GRAV
        ! Update mass density
        do inbor=1,twotondim
           hash_nbor(1:ndim)=ckey(1:ndim,inbor)
           ! Get parent cell using write-only cache
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.false.)
           if(associated(gridp))then
              gridp%rho(icell)=gridp%rho(icell)+mmm*vol(inbor)/vol_loc
#ifdef HYDRO
              if(r%ivar_refine>0)then
                 mask=m%grid(igrid)%uold(ind,r%ivar_refine)/m%grid(igrid)%uold(ind,1)
                 if(mask.gt.r%var_cut_refine)then
                    gridp%nref(icell)=gridp%nref(icell)+mmm*vol(inbor)/r%mass_sph
                 endif
              else
                 gridp%nref(icell)=gridp%nref(icell)+mmm*vol(inbor)/r%mass_sph
              endif
#endif
           end if
        end do
#endif     
     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache(s,m%grid_dict)

  end associate

end subroutine cic_multipole
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_cic_part(pst,input_array,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel,rtype
  integer,dimension(1:input_size)::input_array
  

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CIC_PART,pst%iUpper+1,input_size,0,input_array)
     call r_cic_part(pst%pLower,input_array,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     ilevel=input_array(1)
     rtype=input_array(2)
     if(rtype.ne.3)then
        call cic_part(pst%s,pst%s%p,ilevel)
     endif
     if((pst%s%r%star).and.(rtype.ne.2))then
        call cic_part(pst%s,pst%s%star,ilevel)
     endif
  endif

end subroutine r_cic_part
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine cic_part(s,p,ilevel)
  use amr_parameters, only: ndim,twotondim,dp
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use nbors_utils
  use cache_commons
  use cache
  use multigrid_fine_coarse, only:pack_fetch_phi,unpack_fetch_phi
  use hilbert
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  !
  ! Local variables
  real(dp),dimension(1:ndim)::x,dd,dg
  integer,dimension(1:ndim)::ig,id,ix
  real(dp),dimension(1:twotondim)::vol
  integer,dimension(1:ndim,1:twotondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::i,ipart,icell,ind,idim
  real(kind=8)::dx_loc,vol_loc
  type(oct),pointer::gridp
  type(msg_twin_realdp)::dummy_twin_realdp
  logical::star
  !integer::rtype ! rtype 1 all 2 dm 3 star 4 gas
  
  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Are particles stars?
  star = allocated(p%tp)
  
  ! Compute contribution to multipole
  if(ilevel==r%levelmin)then
     do i=1,p%npart
        g%multipole%q(1)=g%multipole%q(1)+p%mp(i)
     end do
     do idim=1,ndim
        do i=1,p%npart
           g%multipole%q(idim+1)=g%multipole%q(idim+1)+p%mp(i)*p%xp(i,idim)
        end do
     end do
  endif

  ! Sort particle according to current level Hilbert key
  do i=p%headp(ilevel),p%tailp(r%nlevelmax)
     p%sortp(i)=i
  end do
  ix=0
  call sort_hilbert(r,g,p,p%headp(ilevel),p%tailp(r%nlevelmax),ix,0,1,ilevel-1)

  ! Open write-only cache for array rho
  hash_nbor(0)=ilevel+1
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_twin_realdp)/32,&
                pack=pack_fetch_phi,unpack=unpack_fetch_phi,&
                init=init_flush_rho, flush=pack_flush_rho, combine=unpack_flush_rho)

  ! Loop over particles in Hilbert order
  do i=p%headp(ilevel),p%tailp(r%nlevelmax)
     ipart=p%sortp(i)

     ! Rescale particle position at level ilevel
     do idim=1,ndim
        x(idim)=p%xp(ipart,idim)/dx_loc
     end do
     
     ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
     do idim=1,ndim
        dd(idim)=x(idim)+0.5D0
        id(idim)=int(dd(idim))
        dd(idim)=dd(idim)-id(idim)
        dg(idim)=1.0D0-dd(idim)
        ig(idim)=id(idim)-1
     end do
     
     ! Periodic boundary conditions
     do idim=1,ndim
        if(ig(idim)<0)ig(idim)=m%ckey_max(ilevel+1)-1
        if(id(idim)==m%ckey_max(ilevel+1))id(idim)=0
     enddo

     ! Compute cloud volumes
#if NDIM==1
     vol(1)=dg(1)
     vol(2)=dd(1)
#endif
#if NDIM==2
     vol(1)=dg(1)*dg(2)
     vol(2)=dd(1)*dg(2)
     vol(3)=dg(1)*dd(2)
     vol(4)=dd(1)*dd(2)
#endif
#if NDIM==3
     vol(1)=dg(1)*dg(2)*dg(3)
     vol(2)=dd(1)*dg(2)*dg(3)
     vol(3)=dg(1)*dd(2)*dg(3)
     vol(4)=dd(1)*dd(2)*dg(3)
     vol(5)=dg(1)*dg(2)*dd(3)
     vol(6)=dd(1)*dg(2)*dd(3)
     vol(7)=dg(1)*dd(2)*dd(3)
     vol(8)=dd(1)*dd(2)*dd(3)
#endif

     ! Compute cells Cartesian key
#if NDIM==1
     ckey(1,1)=ig(1)
     ckey(1,2)=id(1)
#endif
#if NDIM==2
     ckey(1:2,1)=(/ig(1),ig(2)/)
     ckey(1:2,2)=(/id(1),ig(2)/)
     ckey(1:2,3)=(/ig(1),id(2)/)
     ckey(1:2,4)=(/id(1),id(2)/)
#endif
#if NDIM==3
     ckey(1:3,1)=(/ig(1),ig(2),ig(3)/)
     ckey(1:3,2)=(/id(1),ig(2),ig(3)/)
     ckey(1:3,3)=(/ig(1),id(2),ig(3)/)
     ckey(1:3,4)=(/id(1),id(2),ig(3)/)
     ckey(1:3,5)=(/ig(1),ig(2),id(3)/)
     ckey(1:3,6)=(/id(1),ig(2),id(3)/)
     ckey(1:3,7)=(/ig(1),id(2),id(3)/)
     ckey(1:3,8)=(/id(1),id(2),id(3)/)
#endif

#ifdef GRAV
     ! Update mass density
     do ind=1,twotondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        ! Get parent cell using write-only cache
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.false.)
        if(associated(gridp))then
           gridp%rho(icell)=gridp%rho(icell)+p%mp(ipart)*vol(ind)/vol_loc
           if(star)then
              gridp%nref(icell)=gridp%nref(icell)+p%mp(ipart)*vol(ind)/r%mass_sph
           else
              if(r%mass_cut_refine>0)then
                 if(p%mp(ipart)<r%mass_cut_refine)then
                    gridp%nref(icell)=gridp%nref(icell)+vol(ind)
                 endif
              else
                 gridp%nref(icell)=gridp%nref(icell)+vol(ind)
              endif
           endif
        endif
     end do
#endif

  end do
  ! End loop over particles

  call close_cache(s,m%grid_dict)

end associate
  
end subroutine cic_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_rho(grid,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  type(oct)::grid
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  
  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
#ifdef GRAV
  do ind=1,twotondim
     grid%rho(ind)=0.0
     grid%nref(ind)=0.0
  end do
#endif
  
end subroutine init_flush_rho
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_rho(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_twin_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_twin_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp_phi(ind)=grid%rho(ind)
     msg%realdp_dis(ind)=grid%nref(ind)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_rho
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_rho(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_twin_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_twin_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)
  
#ifdef GRAV
  do ind=1,twotondim
     grid%rho(ind)=grid%rho(ind)+msg%realdp_phi(ind)
     grid%nref(ind)=grid%nref(ind)+msg%realdp_dis(ind)
  end do
#endif

end subroutine unpack_flush_rho
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_split_part(pst,input_array,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel,rtype
  integer,dimension(1:input_size)::input_array

  integer::rID
  ilevel=input_array(1)
  rtype=input_array(2)
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SPLIT_PART,pst%iUpper+1,input_size,0,ilevel)
     call r_split_part(pst%pLower,input_array,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     if(rtype.ne.3)then
        call split_part(pst%s,pst%s%p,ilevel)
     endif
     if((pst%s%r%star).and.(rtype.ne.2))then
        call split_part(pst%s,pst%s%star,ilevel)
     endif
  endif

end subroutine r_split_part
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine pack_fetch_split(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_int4)::msg

  do ind=1,twotondim
     if(grid%refined(ind))then
        msg%int4(ind)=1
     else
        msg%int4(ind)=0
     endif
  enddo
  msg_array=transfer(msg,msg_array)
  
end subroutine pack_fetch_split
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine unpack_fetch_split(grid,msg_size,msg_array,hash_key)
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
     if(msg%int4(ind)==1)then
        grid%refined(ind)=.true.
     else
        grid%refined(ind)=.false.
     endif
  enddo

end subroutine unpack_fetch_split
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine split_part(s,p,ilevel)
  use amr_parameters, only: ndim,twotondim,dp,i8b
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use nbors_utils
  use cache_commons
  use hilbert
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  !
  ! Local variables
  real(dp),dimension(1:ndim)::x,xp_tmp,vp_tmp
  integer,dimension(1:ndim)::ii,ix,ix_ref
  integer(kind=8),dimension(0:ndim)::hash_key
  integer::i,ipart,jpart,idim,icell,ilev
  integer::npart_coarse,npart_fine
  real(kind=8)::dx_loc,vol_loc
  real(dp)::mp_tmp
  integer::levelp_tmp
  integer(i8b)::idp_tmp
  type(oct),pointer::gridp
  type(msg_int4)::dummy_int4

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Open read-only cache for array refined
  hash_key(0)=ilevel
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_int4)/32,&
                pack=pack_fetch_split,unpack=unpack_fetch_split)

  ! Loop over particles
  ix_ref=-1
  npart_coarse=0
  do i=p%headp(ilevel),p%tailp(r%nlevelmax)
     ipart=p%sortp(i)

     ! Acquire grid using read-only cache
     ix = int(p%xp(ipart,1:ndim)/(2*dx_loc))
     if(.NOT. ALL(ix.EQ.ix_ref))then
        hash_key(1:ndim)=ix(1:ndim)
        call get_grid(s,hash_key,m%grid_dict,gridp,flush_cache=.false.,fetch_cache=.true.)
        ix_ref=ix
     endif

     ! If particle sits outside current level,
     ! then it is clearly not in a refined cell.
     ! This can happen during second adaptive step
     if(.not.associated(gridp))then
        npart_coarse=npart_coarse+1
        p%levelp(ipart)=-p%levelp(ipart)
     else
        ! Rescale particle position at level ilevel
        do idim=1,ndim
           x(idim)=p%xp(ipart,idim)/dx_loc
        end do
        
        ! Shift particle position to to 2x2x2 grid corner
        do idim=1,ndim
           ii(idim)=int(x(idim)-2*ix_ref(idim))
        end do
        
        ! Compute parent cell index
#if NDIM==1
        icell=1+ii(1)
#endif
#if NDIM==2
        icell=1+ii(1)+2*ii(2)
#endif
#if NDIM==3
        icell=1+ii(1)+2*ii(2)+4*ii(3)
#endif
        ! Increase counter if cell is not refined
        if(.NOT.gridp%refined(icell))then
           npart_coarse=npart_coarse+1
           p%levelp(ipart)=-p%levelp(ipart)
        else
           p%sortp(i)=-p%sortp(i)
        endif
     endif

  end do
  ! End loop over particles

  call close_cache(s,m%grid_dict)

  p%tailp(ilevel)=p%headp(ilevel)+npart_coarse-1
  do ilev=ilevel+1,r%nlevelmax
     p%headp(ilev)=p%tailp(ilevel)+1
     p%tailp(ilev)=p%npart
  end do

  ! Loop over fine level particles
  ! This preserves the initial ordering after partioning
  npart_fine=0
  do ipart=p%headp(ilevel),p%tailp(r%nlevelmax)
     if(p%levelp(ipart)>0)then
        npart_fine=npart_fine+1
        p%workp(ipart)=p%headp(ilevel+1)+npart_fine-1
     endif
  end do

  ! Loop over coarse level particles
  ! This enforces Hilbert ordering after partioning
  npart_coarse=0
  do i=p%headp(ilevel),p%tailp(r%nlevelmax)
     ipart=p%sortp(i)
     if(ipart>0)then
        npart_coarse=npart_coarse+1
        p%workp(ipart)=p%headp(ilevel)+npart_coarse-1
        p%levelp(ipart)=-p%levelp(ipart)
     endif
  end do

  ! Swap particles using new index table
  do ipart=p%headp(ilevel),p%tailp(r%nlevelmax)
     do while(p%workp(ipart).NE.ipart)
        ! Swap new index
        jpart=p%workp(ipart)
        p%workp(ipart)=p%workp(jpart)
        p%workp(jpart)=jpart
        ! Swap positions
        xp_tmp(1:ndim)=p%xp(ipart,1:ndim)
        p%xp(ipart,1:ndim)=p%xp(jpart,1:ndim)
        p%xp(jpart,1:ndim)=xp_tmp(1:ndim)
        ! Swap velocities
        vp_tmp(1:ndim)=p%vp(ipart,1:ndim)
        p%vp(ipart,1:ndim)=p%vp(jpart,1:ndim)
        p%vp(jpart,1:ndim)=vp_tmp(1:ndim)
        ! Swap masses
        mp_tmp=p%mp(ipart)
        p%mp(ipart)=p%mp(jpart)
        p%mp(jpart)=mp_tmp
        ! Swap metallicity
        if(allocated(p%zp))then
           mp_tmp=p%zp(ipart)
           p%zp(ipart)=p%zp(jpart)
           p%zp(jpart)=mp_tmp
        endif
        ! Swap age
        if(allocated(p%tp))then
           mp_tmp=p%tp(ipart)
           p%tp(ipart)=p%tp(jpart)
           p%tp(jpart)=mp_tmp
        endif
        ! Swap ids
        idp_tmp=p%idp(ipart)
        p%idp(ipart)=p%idp(jpart)
        p%idp(jpart)=idp_tmp
        ! Swap levels
        levelp_tmp=p%levelp(ipart)
        p%levelp(ipart)=p%levelp(jpart)
        p%levelp(jpart)=levelp_tmp
     end do
  end do

  end associate

end subroutine split_part
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
recursive subroutine r_collect_multipole(pst,ilevel,input_size,multipole,output_size)
  use mdl_module
  use amr_parameters, only: ndim
  use ramses_commons, only: pst_t
  use amr_commons, only: multipole_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  type(multipole_t)::multipole,next_multipole

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_COLLECT_MULTIPOLE,pst%iUpper+1,input_size,output_size,ilevel)
     call r_collect_multipole(pst%pLower,ilevel,input_size,multipole,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_multipole)
     multipole%q = multipole%q+next_multipole%q
  else
     multipole%q = pst%s%g%multipole%q
  endif

end subroutine r_collect_multipole
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_broadcast_multipole(pst,multipole,input_size)
  use mdl_module
  use amr_parameters, only: ndim
  use ramses_commons, only: pst_t
  use amr_commons, only: multipole_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(multipole_t)::multipole

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BROADCAST_MULTIPOLE,pst%iUpper+1,input_size,0,multipole)
     call r_broadcast_multipole(pst%pLower,multipole,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     pst%s%g%multipole=multipole
     pst%s%g%rho_tot=pst%s%g%multipole%q(1)/pst%s%r%boxlen**ndim
!!!     pst%s%g%rho_tot=0d0 ! For non-periodic BC
  endif

end subroutine r_broadcast_multipole
#endif
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
recursive subroutine sort_hilbert(r,g,p,head_part, tail_part, ix_coarse, cstate_coarse, ilevel, final_level)
  use amr_parameters, only: dp, ndim, twotondim
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  use hilbert, only: next_state_diagram_reverse,one_digit_diagram
  implicit none
  
  type(run_t),intent(in)::r
  type(global_t),intent(in)::g
  type(part_t)::p
  integer, intent(in) :: ilevel, final_level
  integer, intent(in) :: head_part, tail_part
  integer, dimension(1:ndim), intent(in) :: ix_coarse
  integer, intent(in) :: cstate_coarse
  
  ! Description:
  ! This subroutine sort particles along the Hilbert key at the resolution
  ! set by final_level. It should be called first with ilevel=levelmin.
  
  ! Iputs: 
  ! - Head_part and tail_part are head and tail of particle distribution to work on.
  ! - Array sortp must be initialized with sortp(i)=i between head_part and tail_part.
  ! - Cartesian key of coarse cell in which these particles are contained.
  ! - State of the coarse cell for Hilbert ordering
  ! - Current and final level
  
  ! Example: 
  ! ix=(/0,0,0/)
  ! call sort_hilbert(1, npart, ix, 0, 1, nlevelmax) 
  ! will sort all particles according to their Hilbert key at levelmax.
  ! On output, array sortp is modified.
  
  ! Local variables
  integer :: ip, ind_part, idim, ipart, new_ipart
  integer :: ckey_max, cstate_fine, ind_cart_part, head_fine, tail_fine
  real(dp) :: ckey_factor
  integer, dimension(1:ndim) :: ix_fine, ix_ref, ix_part
  integer, dimension(0:twotondim-1,1:ndim) :: ix, ix_child
  integer, dimension(0:twotondim-1) :: nstate, sdigit, ind, ind_cart, ind_hilbert
  integer, dimension(0:twotondim-1) :: numb_part, offset
  
  ! Compute particle position to cartesian key factor
  ckey_max = 2**ilevel
  ckey_factor = 2.0**ilevel / dble(r%boxlen)
  
  ! Initial Cartesian offset for fine cells
  do idim = 1, ndim
     ix_ref(idim) = ISHFT(ix_coarse(idim),1)
  end do
  
  ! Compute the Hilbert index for fine cells
  do ip = 0, twotondim-1
     sdigit(ip) = ip
  end do
  
  ! Compute lookup index in state diagrams
  do ip = 0, twotondim-1
     ind(ip) = cstate_coarse * twotondim + sdigit(ip)
  end do
  
  ! Save next state
  do ip = 0, twotondim-1
     nstate(ip) = next_state_diagram_reverse(ind(ip))
  end do
  
  ! Add one integer key digit each
  do idim = 1, ndim
     do ip = 0, twotondim-1
        ix(ip, idim) = one_digit_diagram(ind(ip), idim)
     end do
  end do
  
  ! Compute Cartesian index for children cells
  ind_cart = 0
  do idim = 1, ndim
     do ip = 0, twotondim-1
        ix_child(ip, idim) = ix_ref(idim) + ix(ip, idim)
        ind_cart(ip) = ind_cart(ip) + ix(ip, idim) * 2**(idim-1)
     end do
  end do
  
  ! Compute mapping from Cartesian to Hilbert order
  ind_hilbert = 0
  do ip = 0, twotondim-1
     ind_hilbert(ind_cart(ip))=ip
  end do
  
  ! Count particles per children cell
  numb_part = 0
  do ipart = head_part, tail_part
     ind_part = p%sortp(ipart)
     ind_cart_part = 0
     do idim = 1,ndim
        ix_part(idim) = int(p%xp(ind_part,idim)*ckey_factor) - ix_ref(idim)
        ind_cart_part = ind_cart_part + ix_part(idim) * 2**(idim-1)
     end do
     ip = ind_hilbert(ind_cart_part)
     numb_part(ip) = numb_part(ip) + 1
  end do
  
  offset = head_part-1
  do ip = 1, twotondim-1
     offset(ip) = offset(ip-1) + numb_part(ip-1)
  end do
  
  ! Compute new sortp array
  numb_part = 0
  do ipart = head_part, tail_part
     ind_part = p%sortp(ipart)
     ind_cart_part = 0
     do idim = 1,ndim
        ix_part(idim) = int(p%xp(ind_part,idim)*ckey_factor) - ix_ref(idim)
        ind_cart_part = ind_cart_part + ix_part(idim) * 2**(idim-1)
     end do
     ip = ind_hilbert(ind_cart_part)
     numb_part(ip) = numb_part(ip) + 1
     new_ipart = offset(ip) + numb_part(ip)
     p%workp(new_ipart) = ind_part
  end do
  do ipart = head_part,tail_part
     p%sortp(ipart) = p%workp(ipart)
  end do
  
  ! Recursive call
  if(ilevel < final_level)then
     do ip = 0, twotondim-1
        if(numb_part(ip) > 0)then
           head_fine = offset(ip) + 1
           tail_fine = offset(ip) + numb_part(ip)
           ix_fine(1:ndim) = ix_child(ip,1:ndim)
           cstate_fine = nstate(ip)
           call sort_hilbert(r,g,p,head_fine,tail_fine,ix_fine,cstate_fine,ilevel+1,final_level)
        endif
     end do
  endif
  
end subroutine sort_hilbert
!###############################################
!###############################################
!###############################################
!###############################################
end module rho_fine_module
