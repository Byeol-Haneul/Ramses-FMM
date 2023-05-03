module feedback_module
contains
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
recursive subroutine r_feedback(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_FEEDBACK,pst%iUpper+1,input_size,0,ilevel)
     call r_feedback(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call thermal_feedback(pst%s,pst%s%s,ilevel)
  endif

end subroutine r_feedback
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine thermal_feedback(s,p,ilevel)
  use amr_parameters, only: ndim,twotondim,dp
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use nbors_utils
  use cache_commons
  use cache
  use marshal, only: pack_fetch_refine,unpack_fetch_refine
  use boundaries, only: init_bound_refine
  use godunov_fine_module, only: init_flush_godunov,pack_flush_godunov,unpack_flush_godunov
  use hilbert
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  ! Local variables
  integer,dimension(1:ndim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_cell
  integer::i,ipart,icell,ind,idim
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx_loc,vol_loc,vol_cell
  real(dp)::mejecta,mloss,mzloss,zloss,ekinetic,ethermal
  real(dp)::birth_time,t_sn,e_sn,dteff,dold
  type(oct),pointer::gridp
  type(msg_large_realdp)::dummy_large_realdp
  logical::ok_level,ok_leaf

#ifdef HYDRO
#if NDIM==3
  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Supernovae progenitors life time from Myr to proper time in code units
  t_sn=r%t_SN*1d6*(365.*24.*3600.)/(scale_t/g%aexp**2)

  ! Supernovae specific energy from cgs to code units
  e_sn=r%e_SN/(10d0*2d33)/scale_v**2

  ! Open cache for array uold (fetch) and unew (flush)
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_large_realdp)/32,&
                pack=pack_fetch_refine,unpack=unpack_fetch_refine,&
                init=init_flush_godunov, flush=pack_flush_godunov,&
                combine=unpack_flush_godunov, bound=init_bound_refine)

  ! Loop over particles in Hilbert order
  do ipart=p%headp(ilevel),p%tailp(ilevel)

     ! Compute time step for that particle
     dteff=g%dtnew(p%levelp(ipart))*g%aexp**2

     ! Select only recently formed stars
     birth_time=p%tp(ipart) ! Proper time
     if(birth_time.lt.(g%texp-t_sn-dteff))cycle ! Already exploded
     if(birth_time.ge.(g%texp-t_sn))cycle ! Not old enough

     ok_level=.true.

     ! Find parent cell at level ilevel
     do idim=1,ndim
        ckey(idim)=int(p%xp(ipart,idim)/dx_loc)
     end do

     ! Cell volume at level ilevel
     vol_cell=vol_loc

     ! Get parent cell at level ilevel using cache
     hash_cell(0)=ilevel+1
     hash_cell(1:ndim)=ckey(1:ndim)
     call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.true.)

     ! If cell does not exist at current level, then find cell at coarser level
     if(.not.associated(gridp))then

        ! NGP at level ilevel-1
        do idim=1,ndim
           ckey(idim)=int(p%xp(ipart,idim)/dx_loc/2)
        end do

        ! Cell volume at level ilevel-1
        vol_cell=vol_loc*2**ndim

        ! Get parent cell at level ilevel-1 using cache
        hash_cell(0)=ilevel
        hash_cell(1:ndim)=ckey(1:ndim)
        call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.true.)
        if(.not.associated(gridp))ok_level=.false.

     end if

     if(.not. ok_level)then
        write(*,*)"Something went wrong in thermal_feedback"
        write(*,*)"Current level grid and coarser grid both dont exist..."
        stop
     endif

     ok_leaf = .not. gridp%refined(icell)
     if(.not. ok_leaf)then
        write(*,*)"Something went wrong in thermal_feedback"
        write(*,*)"Cell should be a leaf cell..."
        stop
     endif

     ! Compute supernova properties
     mejecta=r%eta_SN*p%mp(ipart)
     mloss=mejecta/vol_cell
     ethermal=mloss*e_sn
     ekinetic=mloss*0.5d0*(p%vp(ipart,1)**2+p%vp(ipart,2)**2+p%vp(ipart,3)**2)
     zloss=r%yield_SN+(1d0-r%yield_SN)*p%zp(ipart)
     mzloss=mloss*zloss

     ! Update unew
     gridp%unew(icell,1)=gridp%unew(icell,1)+mloss
     gridp%unew(icell,2)=gridp%unew(icell,2)+mloss*p%vp(ipart,1)
     gridp%unew(icell,3)=gridp%unew(icell,3)+mloss*p%vp(ipart,2)
     gridp%unew(icell,4)=gridp%unew(icell,4)+mloss*p%vp(ipart,3)
     gridp%unew(icell,5)=gridp%unew(icell,5)+ekinetic+ethermal
     if(r%metal)gridp%unew(icell,r%imetal)=gridp%unew(icell,r%imetal)+mzloss

     ! If dual energy scheme is activated, update entropy
     if(r%entropy.and.r%dual_energy.GE.0)then
        dold = gridp%uold(icell,1)
        gridp%unew(icell,r%ientropy)=gridp%unew(icell,r%ientropy)+ethermal/dold**(r%gamma-1)*(r%gamma-1)
     endif

     ! Update particle mass
     p%mp(ipart)=p%mp(ipart)-mejecta

  end do
  ! End loop over particles

  call close_cache(s,m%grid_dict)

end associate
#endif  
#endif
end subroutine thermal_feedback
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module feedback_module
