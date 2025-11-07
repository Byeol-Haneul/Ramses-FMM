module source_hydro_fine_module
contains
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_source_hydro_fine(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use hydro_parameters, only: nener
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID
  logical::ok

  ! Check if hydro source terms are required
  ok = pst%s%r%entropy .and.  pst%s%r%dual_energy .GE. 0
  ok = ok .or. nener>0
  ok = ok .or. pst%s%r%sgs_turb
  if(.not. ok)return
  
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SOURCE_HYDRO_FINE,pst%iUpper+1,input_size,0,ilevel)
     call r_source_hydro_fine(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call source_hydro_fine(pst%s,ilevel)
  endif

end subroutine r_source_hydro_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine source_hydro_fine(s,ilevel)
  use amr_parameters, only: ndim, twondim, twotondim, threetondim, nvector
  use amr_commons, only: oct, nbor
  use ramses_commons, only: ramses_t
  use hydro_parameters, only: nvar, nener
  use hydro_flag_module, only: pack_fetch_hydro, unpack_fetch_hydro
  use cache_commons
  use cache
  use nbors_utils
  use boundaries, only: init_bound_refine
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !--------------------------------------------------------------
  ! Add all hydro source terms to unew. This includes pdV terms
  ! for non-thermal energies and for turbulent energy if present.
  ! The dual energy fix for hich Mach number flow is also
  ! peerformed here if required.
  !--------------------------------------------------------------
  integer,dimension(1:3,1:8),save::iii=reshape(&
       & (/0,0,0,1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1/),(/3,8/))
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer::igrid,ind,idim,jdim,ivar,irad,i_nbor
  integer::igridd,igridg,icelld,icellg,igridp,icellp
  integer,dimension(1:twondim)::igridn,icelln
  real(kind=8),dimension(1:twondim)::dxn
  integer(kind=8),dimension(0:ndim)::hash_key,hash_nbor
  real(kind=8),dimension(1:ndim,1:ndim)::vg,vd,gradu,E
  real(kind=8),dimension(1:ndim)::dxg,dxd
  type(nbor),dimension(1:twondim)::gridn
  type(oct),pointer::gridp
  type(msg_realdp)::dummy_realdp
  real(kind=8)::dx,phi_diss,div,divu
  real(kind=8)::d,u,v,w,bx,by,bz,d_old,sigma
  real(kind=8)::e_kin,e_mag,e_cons,e_prim,e_turb,e_trunc,T2_cons,T2_fix
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

#ifdef HYDRO

  associate(r=>s%r,g=>s%g,m=>s%m)    

  if(r%verbose.and.g%myid==1)write(*,'("   Entering source_hydro_fine for level ",I2)')ilevel

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum temperature for dual energy switch
  if(r%T2_fix>0)then
     T2_fix=r%T2_fix/scale_T2
  else
     T2_fix=1d100
  endif

  hash_key(0)=ilevel+1
  dx=r%boxlen/2**ilevel

  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
       hilbert=m%domain, pack_size=storage_size(dummy_realdp)/32,&
       pack=pack_fetch_hydro,unpack=unpack_fetch_hydro,&
       bound=init_bound_refine)

  ! Loop over active grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Loop over neighboring cells
     do ind=1,twotondim

        ! Compute neughbring cell hash key
        hash_key(1:ndim)=2*m%grid(igrid)%ckey(1:ndim)+iii(1:ndim,ind)

        ! If a neighbor cell does not exist, replace it by its father cell
        do i_nbor=1,twondim
           hash_nbor(0)=hash_key(0)
           ! Periodic boundary conditions
           do idim=1,ndim
              hash_nbor(idim)=hash_key(idim)+shift(idim,i_nbor)
              if(r%periodic(idim))then
                 if(hash_nbor(idim)< m%box_ckey_min(idim,ilevel+1))hash_nbor(idim)=m%box_ckey_max(idim,ilevel+1)-1
                 if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel+1))hash_nbor(idim)=m%box_ckey_min(idim,ilevel+1)
              endif
           enddo
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icellp,flush_cache=.false.,fetch_cache=.true.)
           if(associated(gridp))then
              gridn(i_nbor)%p=>gridp
              icelln(i_nbor)=icellp
              dxn(i_nbor)=dx
           else
              hash_nbor(0)=hash_nbor(0)-1
              hash_nbor(1:ndim)=hash_nbor(1:ndim)/2
              call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icellp,flush_cache=.false.,fetch_cache=.true.)
              gridn(i_nbor)%p=>gridp
              icelln(i_nbor)=icellp
              dxn(i_nbor)=1.5*dx
           endif
           ! Lock grid in cache
           call lock_cache(s,gridn(i_nbor)%p)
        end do

        ! Loop over dimensions
        do idim=1,ndim
           ! Gather 3-velocity
           icellg=icelln(2*idim-1)
           icelld=icelln(2*idim  )
           vg(idim,1:ndim)=gridn(2*idim-1)%p%uold(icellg,2:ndim+1)/gridn(2*idim-1)%p%uold(icellg,1)
           vd(idim,1:ndim)=gridn(2*idim  )%p%uold(icelld,2:ndim+1)/gridn(2*idim  )%p%uold(icelld,1)
           dxg(idim)=dxn(2*idim-1)
           dxd(idim)=dxn(2*idim  )
        end do
        
        do i_nbor=1,twondim
           call unlock_cache(s,gridn(i_nbor)%p)
        end do

        ! Compute divu = Trace G
        divu=0.0d0
        do idim=1,ndim
           divu = divu + (vd(idim,idim)-vg(idim,idim)) / (dxg(idim)+dxd(idim))
        end do

        ! Add -pdV term to non-thermal energies
#if NENER>0
        do irad=1,nener
           m%grid(igrid)%unew(ind,5+irad)=m%grid(igrid)%unew(ind,5+irad) &
                & -(r%gamma_rad(irad)-1.0d0)*m%grid(igrid)%uold(ind,5+irad)*divu*g%dtnew(ilevel)
        end do
#endif     
        ! Correct total energy if internal energy is too small
        if(r%entropy.and.r%dual_energy.GE.0)then
           d=max(dble(m%grid(igrid)%unew(ind,1)),r%smallr)
           u=m%grid(igrid)%unew(ind,2)/d
           v=m%grid(igrid)%unew(ind,3)/d
           w=m%grid(igrid)%unew(ind,4)/d
           e_kin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e_kin=e_kin+m%grid(igrid)%unew(ind,5+irad)
           end do
#endif
           e_mag=0.0d0
#ifdef MHD
           bx=0.5d0*(m%grid(igrid)%bnew(ind,1)+m%grid(igrid)%bnew(ind,4))
           by=0.5d0*(m%grid(igrid)%bnew(ind,2)+m%grid(igrid)%bnew(ind,5))
           bz=0.5d0*(m%grid(igrid)%bnew(ind,3)+m%grid(igrid)%bnew(ind,6))
           e_mag=0.5d0*(bx**2+by**2+bz**2)
#endif
           e_cons=m%grid(igrid)%unew(ind,5)-e_kin-e_mag
           T2_cons=(r%gamma-1)*e_cons/d
           e_prim=m%grid(igrid)%unew(ind,r%ientropy)*d**(r%gamma-1.0)/(r%gamma-1.0)
           div=abs(divu)*dx
           e_trunc=r%dual_energy*d*max(abs(divu)*dx,3.0d0*g%hexp*dx)**2
           if(e_cons<e_trunc.AND.T2_cons<T2_fix)then
              m%grid(igrid)%unew(ind,5)=e_prim+e_kin+e_mag
           else
              m%grid(igrid)%unew(ind,r%ientropy)=e_cons/d**(r%gamma-1.0)*(r%gamma-1.0)
           end if
        end if

        ! Add source terms for subgrid turbulence model
        if(r%sgs_turb)then
           ! Compute gradu = G
           gradu(1:ndim,1:ndim)=0.0d0
           do idim=1,ndim
              do jdim=1,ndim
                 gradu(idim,jdim) = (vd(idim,jdim)-vg(idim,jdim)) / (dxg(idim)+dxd(idim))
              enddo
           end do
           
           ! Compute 0.5*(gradu+gradu^T) = E
           E(1:ndim,1:ndim)=0.0d0
           do idim=1,ndim
              do jdim=1,ndim
                 E(idim,jdim) = 0.5*(gradu(idim,jdim)+gradu(jdim,idim))
              enddo
           end do
           
           ! Compute turbulent viscosity dissipation function
           phi_diss=0
           do idim=1,ndim
              do jdim=1,ndim
                 phi_diss = phi_diss+2.0*E(idim,jdim)**2
              enddo
           enddo
           phi_diss=phi_diss-2.0/3.0*divu**2
           
           ! Add -pdV term to subgrid turbulent energy
           m%grid(igrid)%unew(ind,r%iturb)=m%grid(igrid)%unew(ind,r%iturb) &
                & -(2.0/3.0)*m%grid(igrid)%uold(ind,r%iturb)*divu*g%dtnew(ilevel)

           if(g%nstep_coarse==0) then
              ! Initialize subgrid turbulent energy with stationary solution
              m%grid(igrid)%unew(ind,r%iturb)=m%grid(igrid)%uold(ind,1)*dx**2*phi_diss
           else
              ! Implicit solution wrt to decay term only
              d_old=max(dble(m%grid(igrid)%uold(ind,1)),r%smallr)
              e_turb=m%grid(igrid)%uold(ind,r%iturb)
              sigma=sqrt(max(2.0*e_turb/d_old,dble(r%smallc)**2))
              m%grid(igrid)%unew(ind,r%iturb)=(m%grid(igrid)%unew(ind,r%iturb) &
                   &  +d_old*dx*sigma*phi_diss*g%dtnew(ilevel)) &
                   & /(1.0+sigma/dx*g%dtnew(ilevel))
           end if

        endif
        
     end do
     ! End loop over cells
  end do
  ! End loop over grids

  call close_cache(s,m%grid_dict)

  end associate

#endif

end subroutine source_hydro_fine
!################################################################
!################################################################
!################################################################
!################################################################
end module source_hydro_fine_module
