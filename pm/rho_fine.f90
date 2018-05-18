!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine m_rho_fine(pst,ilevel)
  use amr_parameters, only: dp,ndim
  use ramses_commons, only: pst_t
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
  real(kind=8),dimension(1:ndim+1)::multipole_tot
  integer,dimension(:),allocatable::input_array,output_array
  integer::i,input_size

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)

  if(.not. r%poisson)return
  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,'(" Entering rho_fine for level ",I2)')ilevel

  !---------------------------
  ! Reset multipole to zero
  !---------------------------
  if(ilevel==r%levelmin)then
     multipole_tot=0d0

     input_size=(storage_size(multipole_tot)/32)*(ndim+1)
     allocate(input_array(1:input_size))
     input_array=transfer(multipole_tot,input_array)
     call r_broadcast_multipole(pst,input_size,0,input_array)
     deallocate(input_array)

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
           call r_multipole_leaf_cells(pst,1,0,i)
        endif

        ! Average down multipoles in all split cells
        if(i<r%nlevelmax)then
           if(m%noct_tot(i+1)>0)then
              if(r%verbose)write(*,'(" Compute split multipoles for level ",I2)')i        
              call r_multipole_split_cells(pst,1,0,i)
           endif
        endif

     endif

     ! Reset array rho to zero
     if(m%noct_tot(i)>0)then
        call r_reset_rho(pst,1,0,i)
     endif

        ! Gas mass deposition using pseudo-particles
     if(r%hydro.AND.m%noct_tot(i)>0)then
        if(r%verbose)write(*,'(" Compute rho from multipoles for level ",I2)')i
        call r_cic_multipole(pst,1,0,i)
     endif

  end do
  ! End loop over finer levels

  !-------------------------------------------------------
  ! Compute particle contribution to density field
  !-------------------------------------------------------
  if(r%pic)then
     do i=ilevel,r%nlevelmax
        if(m%noct_tot(i)>0)then
           call r_cic_part(pst,1,0,i)
        endif
        if(m%noct_tot(i)>0.AND.i<r%nlevelmax)then
           call r_split_part(pst,1,0,i)
        endif
     end do
  endif

  !---------------------------------------------------------------------
  ! Collect multipole contribution from all CPU and broadcast rho_tot
  !---------------------------------------------------------------------
  if(ilevel==r%levelmin)then

     ! Collect local multipole from all CPU
     input_size=(storage_size(multipole_tot)/32)*(ndim+1)
     allocate(output_array(1:input_size))
     call r_collect_multipole(pst,1,input_size,ilevel,output_array)
     multipole_tot=transfer(output_array,multipole_tot)
     deallocate(output_array)

     ! Broadcast total multipole to all CPU
     input_size=(storage_size(multipole_tot)/32)*(ndim+1)
     allocate(input_array(1:input_size))
     input_array=transfer(multipole_tot,input_array)
     call r_broadcast_multipole(pst,(ndim+1)*(storage_size(multipole_tot)/32),0,input_array)     
     deallocate(input_array)

     if(r%verbose)write(*,*)'rho_average=',g%rho_tot
  endif  

  end associate

end subroutine m_rho_fine
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_multipole_leaf_cells(pst,input_size,output_size,ilevel)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer::ilevel

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_MULTIPOLE_LEAF_CELLS,pst%iUpper+1,input_size,output_size,ilevel)
     call r_multipole_leaf_cells(pst%pLower,input_size,output_size,ilevel)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
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
  integer::igrid,ind,idim,ivar,nstride,ioct,icell
  real(dp),dimension(1:ndim),save::xx
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
              xx(idim)=(2*m%grid(igrid)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx_loc
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
recursive subroutine r_multipole_split_cells(pst,input_size,output_size,ilevel)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer::ilevel

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_MULTIPOLE_SPLIT_CELLS,pst%iUpper+1,input_size,output_size,ilevel)
     call r_multipole_split_cells(pst%pLower,input_size,output_size,ilevel)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
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
  use ramses_commons, only: ramses_t
  use cache_commons
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !-------------------------------------------------------------------
  ! This routine compute the monopole and dipole of the gas mass and
  ! the analytical profile (if any) within each cell.
  ! For pure particle runs, this is not necessary and the
  ! routine is not even called.
  !-------------------------------------------------------------------
  integer::igrid,ind,idim,ivar,nstride,ioct,icell
  integer::parent_cell
  integer,external::get_parent_cell
  real(kind=8)::average
  integer(kind=8),dimension(0:ndim)::hash_key
  logical::leaf_cell

  associate(r=>s%r,g=>s%g,m=>s%m)
  
  !-------------------------------------------------------
  ! Perform octree restriction from level ilevel+1
  !-------------------------------------------------------
  call open_cache(s,operation_multipole,domain_decompos_amr)

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ioct=m%head(ilevel+1),m%tail(ilevel+1)
     hash_key(1:ndim)=m%grid(ioct)%ckey(1:ndim)
     ! Get parent cell using a write-only cache
     parent_cell=get_parent_cell(s,hash_key,m%grid_dict,.true.,.false.)
     igrid=(parent_cell-1)/twotondim+1
     icell=parent_cell-(igrid-1)*twotondim
#ifdef HYDRO
     ! Average conservative variables
     do ivar=1,ndim+1
        average=0.0d0
        do ind=1,twotondim
           average=average+m%grid(ioct)%unew(ind,ivar)
        end do
        ! Scatter result to cell
        m%grid(igrid)%unew(icell,ivar)=average
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
subroutine init_flush_multipole(grid,msg_size)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  type(oct)::grid
  integer::msg_size

  integer::ind,ivar
  
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
  integer,dimension(1:msg_size)::msg_array

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
subroutine unpack_flush_multipole(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size)::msg_array

  integer::ind,ivar
  type(msg_realdp)::msg

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
recursive subroutine r_reset_rho(pst,input_size,output_size,ilevel)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer::ilevel

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_RESET_RHO,pst%iUpper+1,input_size,output_size,ilevel)
     call r_reset_rho(pst%pLower,input_size,output_size,ilevel)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
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
     end do
  end do
#endif  

end subroutine reset_rho
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_cic_multipole(pst,input_size,output_size,ilevel)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer::ilevel

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_CIC_MULTIPOLE,pst%iUpper+1,input_size,output_size,ilevel)
     call r_cic_multipole(pst%pLower,input_size,output_size,ilevel)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     call cic_multipole(pst%s,ilevel)
  endif

end subroutine r_cic_multipole
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cic_multipole(s,ilevel)
  use amr_parameters, only: ndim,twotondim,dp
  use ramses_commons, only: ramses_t
  use cache_commons
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !
  ! Local variables
  real(dp),dimension(1:ndim),save::x,dd,dg
  integer,dimension(1:ndim),save::ig,id
  real(dp),dimension(1:twotondim),save::vol
  integer,dimension(1:ndim,1:twotondim)::ckey
  integer(kind=8),dimension(0:ndim),save::hash_nbor
  integer::inbor,igrid,ind,idim
  integer::ioct,icell,parent_cell
  integer,external::get_parent_cell
  real(kind=8)::dx_loc,vol_loc,mmm

  associate(r=>s%r,g=>s%g,m=>s%m)
    
  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Use hash table directly for cells (not for grids)
  hash_nbor(0)=ilevel+1

  call open_cache(s,operation_rho,domain_decompos_amr)

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
           call mdl_abort
        endif
        x(1:ndim)=m%grid(igrid)%unew(ind,2:ndim+1)/mmm
        
        ! Compute total multipole
        if(ilevel==r%levelmin)then
           do idim=1,ndim+1
              g%multipole(idim)=g%multipole(idim)+m%grid(igrid)%unew(ind,idim)
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
           parent_cell=get_parent_cell(s,hash_nbor,m%grid_dict,.true.,.false.)
           if(parent_cell>0)then
              ioct=(parent_cell-1)/twotondim+1
              icell=parent_cell-(ioct-1)*twotondim
              m%grid(ioct)%rho(icell)=m%grid(ioct)%rho(icell)+mmm*vol(inbor)/vol_loc
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
recursive subroutine r_cic_part(pst,input_size,output_size,ilevel)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer::ilevel

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_CIC_PART,pst%iUpper+1,input_size,output_size,ilevel)
     call r_cic_part(pst%pLower,input_size,output_size,ilevel)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     call cic_part(pst%s,ilevel)
  endif

end subroutine r_cic_part
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine cic_part(s,ilevel)
  use amr_parameters, only: ndim,twotondim,dp
  use ramses_commons, only: ramses_t
  use cache_commons
  use hilbert
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !
  ! Local variables
  real(dp),dimension(1:ndim),save::x,dd,dg
  integer,dimension(1:ndim),save::ig,id,ix
  real(dp),dimension(1:twotondim),save::vol
  integer,dimension(1:ndim,1:twotondim),save::ckey
  integer(kind=8),dimension(0:ndim),save::hash_nbor
  integer::i,ipart,igrid,ind,idim
  integer::icell,parent_cell
  integer,external::get_parent_cell
  real(kind=8)::dx_loc,vol_loc,vol2
  
  associate(r=>s%r,g=>s%g,m=>s%m,p=>s%p)

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Compute contribution to multipole
  if(ilevel==r%levelmin)then
     do i=1,p%npart
        g%multipole(1)=g%multipole(1)+p%mp(i)
     end do
     do idim=1,ndim
        do i=1,p%npart
           g%multipole(idim+1)=g%multipole(idim+1)+p%mp(i)*p%xp(i,idim)
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
  call open_cache(s,operation_rho,domain_decompos_amr)

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
        parent_cell=get_parent_cell(s,hash_nbor,m%grid_dict,.true.,.false.)
        if(parent_cell>0)then
           igrid=(parent_cell-1)/twotondim+1
           icell=parent_cell-(igrid-1)*twotondim
           vol2=p%mp(ipart)*vol(ind)/vol_loc
           m%grid(igrid)%rho(icell)=m%grid(igrid)%rho(icell)+vol2
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
subroutine init_flush_rho(grid,msg_size)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  type(oct)::grid
  integer::msg_size

  integer::ind
  
#ifdef GRAV
  do ind=1,twotondim
     grid%rho(ind)=0.0
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
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size)::msg_array

  integer::ind
  type(msg_small_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp(ind)=grid%rho(ind)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_rho
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_rho(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size)::msg_array

  integer::ind
  type(msg_small_realdp)::msg

  msg=transfer(msg_array,msg)
  
#ifdef GRAV
  do ind=1,twotondim
     grid%rho(ind)=grid%rho(ind)+msg%realdp(ind)
  end do
#endif

end subroutine unpack_flush_rho
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_split_part(pst,input_size,output_size,ilevel)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer::ilevel

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_SPLIT_PART,pst%iUpper+1,input_size,output_size,ilevel)
     call r_split_part(pst%pLower,input_size,output_size,ilevel)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     call split_part(pst%s,ilevel)
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
  integer,dimension(1:msg_size)::msg_array

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
subroutine unpack_fetch_split(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size)::msg_array

  integer::ind
  type(msg_int4)::msg

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
subroutine split_part(s,ilevel)
  use amr_parameters, only: ndim,twotondim,dp,i8b
  use ramses_commons, only: ramses_t
  use cache_commons
  use hilbert
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !
  ! Local variables
  real(dp),dimension(1:ndim),save::x,xp_tmp,vp_tmp
  integer,dimension(1:ndim),save::ii,ix,ix_ref
  integer(kind=8),dimension(0:ndim),save::hash_key
  integer::i,ipart,jpart,igrid,idim,icell
  integer,external::get_grid
  integer::npart_coarse,npart_fine
  real(kind=8)::dx_loc,vol_loc
  real(dp)::mp_tmp
  integer::levelp_tmp
  integer(i8b)::idp_tmp

  associate(r=>s%r,g=>s%g,m=>s%m,p=>s%p,mdl=>s%mdl)
    
  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Open read-only cache for array refined
  hash_key(0)=ilevel
  call open_cache(s,operation_split,domain_decompos_amr)

  ! Loop over particles
  ix_ref=-1
  igrid=0
  npart_coarse=0
  do i=p%headp(ilevel),p%tailp(r%nlevelmax)
     ipart=p%sortp(i)

     ! Acquire grid using read-only cache
     ix = int(p%xp(ipart,1:ndim)/(2*dx_loc))
     if(.NOT. ALL(ix.EQ.ix_ref))then
        hash_key(1:ndim)=ix(1:ndim)
        igrid=get_grid(s,hash_key,m%grid_dict,.false.,.true.)
        ix_ref=ix
     endif

     ! If particle sits outside current level,
     ! then it is clearly not in a refined cell.
     ! This can happen during second adaptive step
     if(igrid==0)then
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
        
        ! Compute parent cell
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
        if(.NOT.m%grid(igrid)%refined(icell))then
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
  p%headp(ilevel+1)=p%tailp(ilevel)+1

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
  ! set by final_level. It should be called first with ilevel=1
  ! arrays in memory and that positions, 3-integer hilbert keys
  ! and next_state are allocated as particle-based quantities.
  
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
recursive subroutine r_collect_multipole(pst,input_size,output_size,ilevel,output_array)
  use amr_parameters, only: ndim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer::ilevel
  integer,dimension(1:output_size)::output_array

  integer,dimension(1:output_size)::next_output_array
  real(kind=8),dimension(1:ndim+1)::multipole,next_multipole

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_COLLECT_MULTIPOLE,pst%iUpper+1,input_size,output_size,ilevel)
     call r_collect_multipole(pst%pLower,input_size,output_size,ilevel,output_array)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size,next_output_array)
     multipole=transfer(output_array,multipole)
     next_multipole=transfer(next_output_array,next_multipole)
     multipole=multipole+next_multipole
     output_array=transfer(multipole,output_array)
  else
     output_array=transfer(pst%s%g%multipole,output_array)
  endif

end subroutine r_collect_multipole
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_broadcast_multipole(pst,input_size,output_size,input_array)
  use amr_parameters, only: ndim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_BROADCAST_MULTIPOLE,pst%iUpper+1,input_size,output_size,input_array)
     call r_broadcast_multipole(pst%pLower,input_size,output_size,input_array)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     pst%s%g%multipole=transfer(input_array,pst%s%g%multipole)
     pst%s%g%rho_tot=pst%s%g%multipole(1)/pst%s%r%boxlen**ndim
!!!     pst%s%g%rho_tot=0d0 ! For non-periodic BC
  endif

end subroutine r_broadcast_multipole
!###############################################
!###############################################
!###############################################
!###############################################

