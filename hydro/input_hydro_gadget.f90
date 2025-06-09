module input_hydro_gadget_module
contains
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_input_hydro_gadget(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size

  integer::ilevel
  integer::rID
  
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_HYDRO_GADGET,pst%iUpper+1,input_size,0,ilevel)
     call r_input_hydro_gadget(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call input_hydro_gadget(pst%s,ilevel)
  endif
  
end subroutine r_input_hydro_gadget
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_hydro_gadget(s,ilevel)
  use amr_parameters, only: ndim,twotondim,dp
  use amr_commons, only: oct
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use godunov_fine_module, only: init_flush_godunov,pack_flush_godunov,unpack_flush_godunov
  use input_hydro_condinit_module, only: input_hydro_vecpot,cons_from_prim
  use marshal, only: pack_fetch_refine,unpack_fetch_refine
  use boundaries, only: init_bound_refine
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !-----------------------------------------------------
  ! Compute initial conditions from gadget gas particles
  !-----------------------------------------------------
  real(dp),dimension(1:ndim)::x,dd,dg
  integer,dimension(1:ndim)::ig,id,ix
  real(dp),dimension(1:twotondim)::vol
  integer,dimension(1:ndim,1:twotondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::igrid,ipart,icell,ind,idim,ivar
  real(kind=8)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,scale_m
  real(kind=8)::dx_loc,vol_loc,ekin
  type(oct),pointer::gridp
  type(msg_large_realdp)::dummy_large_realdp

  associate(r=>s%r,g=>s%g,m=>s%m,p=>s%gas)

  if(m%noct(ilevel)==0)return

#ifdef HYDRO
  !------------------------------------
  ! Reset conservative variables unew
  !------------------------------------
  do igrid=m%head(ilevel),m%tail(ilevel)
     m%grid(igrid)%unew=0.0
  end do

  !----------------------------------------------
  ! Deposit gas particle variables to the grid
  ! using cloud-in-cell. We deposit only mass,
  ! momentum, internal energy and metallicity.
  !----------------------------------------------

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel
  vol_loc=dx_loc**ndim

  ! Open write-only cache for array rho
  hash_nbor(0)=ilevel+1
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
       & hilbert=m%domain,pack_size=storage_size(dummy_large_realdp)/32,&
       & pack=pack_fetch_refine,unpack=unpack_fetch_refine,&
       & init=init_flush_godunov,flush=pack_flush_godunov,&
       & combine=unpack_flush_godunov,bound=init_bound_refine)

  ! Loop over gas particles
  do ipart=1,p%npart

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
     ! Update mass, momentum, specific internal energy and metallicity
     do ind=1,twotondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        ! Get parent cell using write-only cache
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.false.)
        if(associated(gridp))then
           gridp%unew(icell,1)=gridp%unew(icell,1)+p%mp(ipart)*vol(ind)/vol_loc
           do idim=1,3
              gridp%unew(icell,idim+1)=gridp%unew(icell,idim+1)+p%mp(ipart)*p%vp(ipart,idim)*vol(ind)/vol_loc
           end do
           gridp%unew(icell,5)=gridp%unew(icell,5)+p%mp(ipart)*p%up(ipart)*vol(ind)/vol_loc
           if(r%metal) then
              gridp%unew(icell,r%imetal)=gridp%unew(icell,r%imetal)+p%mp(ipart)*p%zp(ipart)*vol(ind)/vol_loc
           endif
        endif
     end do
     ! End loop over cloud
  end do
  ! End loop over particles

  call close_cache(s,m%grid_dict)

  !--------------------
  ! Set uold to unew
  !--------------------
  do igrid=m%head(ilevel),m%tail(ilevel)
     m%grid(igrid)%uold=m%grid(igrid)%unew
  end do

  !----------------------------------------------
  ! Set empty cells to minumum values from
  ! the namelist IG medium variables.
  ! Convert primitive to conservative variables.
  !----------------------------------------------
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Deal with empty cells
  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim
        do ivar=nvar,1,-1
           if(m%grid(igrid)%uold(ind,1)<r%IG_rho/scale_nH)then
              m%grid(igrid)%uold(ind,ivar)=0.0
              if(ivar.eq.1)m%grid(igrid)%uold(ind,ivar)=max(dble(r%IG_rho)/scale_nH,dble(r%smallr))
              if(ivar.eq.5)m%grid(igrid)%uold(ind,ivar)=r%IG_T2/scale_T2/(r%gamma-1)*max(dble(r%IG_rho)/scale_nH,dble(r%smallr))
              if(r%metal)then
                 if(ivar.eq.r%imetal)m%grid(igrid)%uold(ind,ivar)=r%IG_metal*max(dble(r%IG_rho)/scale_nH,dble(r%smallr))
              endif
           endif
        end do
     end do
  end do

  !----------------------------------------------
  ! Compute proper primitive variables.
  !----------------------------------------------
  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim
        ! From momentum to velocity
        do idim=1,3
           m%grid(igrid)%uold(ind,idim+1)=m%grid(igrid)%uold(ind,idim+1)/m%grid(igrid)%uold(ind,1)
        end do
        ! From internal energy density to pressure
        m%grid(igrid)%uold(ind,5)=(r%gamma-1)*m%grid(igrid)%uold(ind,5)
        ! From metal mass density to metal mass fraction
        if(r%metal) then
           m%grid(igrid)%uold(ind,r%imetal)=m%grid(igrid)%uold(ind,r%imetal)/m%grid(igrid)%uold(ind,1)
        endif
        ! Compute entropy from pressure
        if(r%entropy) then
           m%grid(igrid)%uold(ind,r%ientropy)=m%grid(igrid)%uold(ind,5)/m%grid(igrid)%uold(ind,1)**r%gamma
        endif
     end do
  end do

  ! Compute initial magnetic field
  call input_hydro_vecpot(r,g,m,ilevel)

  ! Convert primitive to conservative
  call cons_from_prim(r,g,m,ilevel)
#endif

  end associate

end subroutine input_hydro_gadget
!################################################################
!################################################################
!################################################################
!################################################################
end module input_hydro_gadget_module
