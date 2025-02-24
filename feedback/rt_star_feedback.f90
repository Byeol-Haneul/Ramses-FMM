!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Module for RT stellar feedback
!_________________________________________________________________________
!
MODULE rt_star_feedback
  !_________________________________________________________________________
  use amr_commons, only: run_t, global_t
  use constants, only: L_sun, m_sun, clight, eV2erg, hplanck, sec2Gyr
  use hydro_parameters, only: nion
  use rt_parameters, only: nrtgrp
  implicit none

CONTAINS
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
recursive subroutine r_star_rt_feedback(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_STAR_RT_FEEDBACK,pst%iUpper+1,input_size,0,ilevel)
     call r_star_rt_feedback(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call star_rt_feedback(pst%s, pst%s%star, ilevel)
  endif

end subroutine r_star_rt_feedback
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine star_rt_feedback(s, p, ilevel)
  use amr_parameters, only: ndim, twotondim, dp
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use rt_parameters, only: nrtgrp
  use SED_module, only: getNPhotonsEmitted
  use nbors_utils
  use cache_commons
  use cache
  use hilbert
  implicit none
  type(ramses_t) :: s
  type(part_t) :: p
  integer :: ilevel
  !==================================================================
  ! This is the RAMSES routine for stellar radiation feedback.
  ! The emissivity grid variable is updated using each star particle
  ! luminosity. Energy will be deposited later during RT subcycles.
  !==================================================================
  ! Local variables
  integer,dimension(1:ndim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_cell
  integer::ipart,icell,idim
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::scale_Np,scale_Fp,scale_inp,scale_inp_cell,scale_msun
  real(dp)::dx_loc,vol_loc,vol_cell
  real(dp)::z,mass,age,code2Gyr,dt_Gyr,dt_loc_Gyr,t_SN_Gyr
  type(oct),pointer::gridp
  type(msg_rt_emissivity_realdp)::dummy_rt_emissivity_realdp
  logical::ok_level
  real(dp),dimension(nrtgrp)::part_NpInp, lum


#ifdef RT
#if NDIM==3
  associate(r=>s%r, g=>s%g, m=>s%m)

  ! Mesh spacing in that level
  dx_loc = r%boxlen / 2**ilevel 
  vol_loc = dx_loc**ndim

  ! Conversion factor from user units to cgs units
  call units(r, g, scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
  call rt_units(r, g, scale_Np, scale_Fp)
  scale_inp = r%rt_esc_frac * scale_d / scale_Np / vol_loc / M_sun
  scale_msun = scale_d * scale_l**ndim / M_sun

  ! Proper time (codeunits) to Gyr
  code2Gyr = scale_t * sec2Gyr / g%aexp**2

  ! Time step from code units to Gyr
  dt_Gyr = g%dtnew(ilevel) * scale_t * sec2Gyr

  ! Supernovae progenitors life time from Myr to Gyr
  t_SN_Gyr = r%t_SNII * 1d-3

  ! Open cache for array emissivity (flush)
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32, &
                hilbert=m%domain,pack_size=storage_size(dummy_rt_emissivity_realdp)/32, &
                pack=pack_fetch_emissivity, unpack=unpack_fetch_emissivity, &
                init=init_flush_emissivity, flush=pack_flush_emissivity, &
                combine=unpack_flush_emissivity)

  ! Loop over particles in Hilbert order
  do ipart=p%headp(ilevel),p%tailp(ilevel)

     ok_level=.true.

     ! Find parent cell at level ilevel
     do idim=1,ndim
        ckey(idim)=int(p%xp(ipart,idim)/dx_loc)
     end do

     ! Cell volume at level ilevel
     vol_cell=vol_loc
     scale_inp_cell=scale_inp

     ! Get parent cell at level ilevel using cache
     hash_cell(0)=ilevel+1
     hash_cell(1:ndim)=ckey(1:ndim)
     call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.false.)

     ! If cell does not exist at current level, then find cell at coarser level
     if(.not.associated(gridp))then

        ! NGP at level ilevel-1
        do idim=1,ndim
           ckey(idim)=int(p%xp(ipart,idim)/dx_loc/2)
        end do

        ! Cell volume at level ilevel-1
        vol_cell=vol_loc*2**ndim
        scale_inp_cell=scale_inp/2**ndim
        
        ! Get parent cell at level ilevel-1 using cache
        hash_cell(0)=ilevel
        hash_cell(1:ndim)=ckey(1:ndim)
        call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icell &
                                ,flush_cache=.true.,fetch_cache=.false.)
        if(.not.associated(gridp))ok_level=.false.

     end if

     if(.not. ok_level)then
        write(*,*)"Something went wrong in star_rt_feedback"
        write(*,*)"Current level grid and coarser grid both dont exist..."
        stop
     endif

     ! Compute star particle properties
     if(r%metal) then
        z = max(p%zp(ipart), 1d-5)                     ! [m_metals/m_tot]
     else
        z = max(r%z_ave*0.02, 1d-5)
     endif
     age = (g%texp - p%tp(ipart)) * code2Gyr

     ! Possibilities: Born i) before dt, ii) within dt, iii) after dt:
     dt_loc_Gyr = max(min(dt_Gyr, age), 0.)
     call getNPhotonsEmitted(r, s%SED, age, dt_loc_Gyr, z, part_NpInp(1:nrtgrp))

     mass = p%mp(ipart)
     if(age.gt.t_SN_Gyr) then
        mass = mass / (1d0 - r%eta_SNII)
     endif

     part_NpInp(1:nrtgrp) = part_NpInp(1:nrtgrp) * mass * scale_inp_cell ! #photons cm-3

     if(r%rt_emission_stats) then
        g%step_nPhot = g%step_nPhot + part_NpInp(1) / scale_inp_cell * scale_msun
        g%step_nStar = g%step_nStar + dt_loc_Gyr/sec2Gyr/scale_t
        g%step_mStar = g%step_mStar + p%mp(ipart) * scale_msun &
                                    * dt_loc_Gyr /sec2Gyr / scale_t
     endif

     lum(1:nrtgrp) = 0.
     if(dt_loc_Gyr > 0.)then
        lum(1:nrtgrp) = part_NpInp(1:nrtgrp) / dt_Gyr * sec2Gyr ! #photons cm-3 s-1
     endif
     lum(1:nrtgrp) = lum(1:nrtgrp) * scale_t ! back to code units

     ! Update parent cell emissivity
     gridp%emissivity(icell,1:nrtgrp) = gridp%emissivity(icell,1:nrtgrp) + lum(1:nrtgrp)

  end do
  ! End loop over particles

  call close_cache(s,m%grid_dict)

end associate
#endif
#endif
end subroutine star_rt_feedback
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine pack_fetch_emissivity(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use rt_parameters, only: nrtgrp
  use amr_commons, only: oct
  use cache_commons, only: msg_rt_emissivity_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  type(msg_rt_emissivity_realdp)::msg

#ifdef RT  
  msg%realdp=grid%emissivity
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_emissivity
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine unpack_fetch_emissivity(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use rt_parameters, only: nrtgrp
  use amr_commons, only: oct
  use cache_commons, only: msg_rt_emissivity_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key
  type(msg_rt_emissivity_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef RT
  grid%emissivity=msg%realdp
#endif

end subroutine unpack_fetch_emissivity
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine init_flush_emissivity(grid,hash_key)
  use amr_parameters, only: ndim,twotondim
  use rt_parameters, only: nrtvar
  use amr_commons, only: oct
  type(oct)::grid
  integer(kind=8),dimension(0:ndim)::hash_key

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
#ifdef RT
  grid%emissivity=0.0d0
#endif

end subroutine init_flush_emissivity
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine pack_flush_emissivity(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use rt_parameters, only: nrtgrp
  use amr_commons, only: oct
  use cache_commons, only: msg_rt_emissivity_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  type(msg_rt_emissivity_realdp)::msg

#ifdef RT
  msg%realdp=grid%emissivity
#endif
  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_emissivity
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine unpack_flush_emissivity(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use rt_parameters, only: nrtgrp
  use amr_commons, only: oct
  use cache_commons, only: msg_rt_emissivity_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  type(msg_rt_emissivity_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)
#ifdef RT
  grid%emissivity=grid%emissivity+msg%realdp
#endif

end subroutine unpack_flush_emissivity
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
END MODULE rt_star_feedback
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
