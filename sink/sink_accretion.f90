module sink_accretion_module

   type :: out_accretion_t
      real(kind=8)::mass
   end type out_accretion_t

contains
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################`
recursive subroutine r_sink_accretion(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(out_accretion_t)::output,next_output

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SINK_ACCRETION,pst%iUpper+1,input_size,output_size,ilevel)
     call r_sink_accretion(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
  else
     call sink_accretion(pst%s,pst%s%sink,ilevel,output%mass)
  endif

end subroutine r_sink_accretion
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine sink_accretion(s,p,ilevel,macc_loc)
  use constants
  use amr_parameters, only: ndim,twotondim,dp
  use amr_commons, only: nbor,oct
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t,cross
  use params_module
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
  real(kind=8)::macc_loc
  !==================================================================
  ! This is the RAMSES routine for sink (black hole) particle accretion.
  ! For now, it is focused on a simple mass-weighted Bondi-Hoyle-Lyttleton accretion scheme
  ! The routine modifies hydro variables unew, as well as sink particle properties.
  ! Written by Nicholas Choustikov (Feb 2025)
  ! TODO:
  ! 1. Test
  ! 2. Add other accretion methods, using r%accretion_type
  ! 3. Add more weightings (e.g. Yohan's gaussian)
  ! 4. Add a timestep control based on the accretion rate
  ! 5. Add more safeties as in RAMSES OG
  ! 6. Think of more things
  !==================================================================
  ! Local variables
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,factG ! Units
  logical::ok
  real(dp)::rr,x,y,z,rrad
  integer::nBHnei,iBHnei
  real(dp),dimension(:,:),allocatable::xBHnei
  real(dp)::dx_loc,vol_loc,vol_cell
  real(dp),dimension(1:ndim)::xcen,xnei
  integer,dimension(1:ndim)::ckey,ckey_nbor
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::i,j,k,ipart,icelln,ind,idim,nrad
  real(dp)::d,e,ethermal,r2_sink,v_bondi,cs_sink,cs,rho_sink,velocity
  real(dp)::weight
  real(dp),dimension(1:ndim)::vv,v_rel,x_acc,p_acc,l_acc,vel_sink
  type(oct),pointer::gridn
  real(dp)::dMBH_overdt,dMED_overdt,m_acc,d_acc
  type(msg_large_realdp)::dummy_large_realdp

#ifdef HYDRO
#if NDIM==3

  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Gravitational constant
  factG=1
  if(r%cosmo)factG=3d0/4d0/twopi*g%omega_m*g%aexp

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Compute number of cells within sink sphere
  rrad = r%sink_accretion_radius/dx_loc
  nrad = int(rrad)
  nBHnei = 0
  do k = -nrad, nrad
     z = dble(k) / dble(nrad) * rrad
     do j = -nrad, nrad
        y = dble(j) / dble(nrad) * rrad
        do i = -nrad, nrad
           x = dble(i) / dble(nrad) * rrad
           rr = sqrt(dble(x*x+y*y+z*z))
           if(rr .le. rrad) nBHnei = nBHnei + 1
        enddo
     enddo
  enddo
  allocate(xBHnei(1:ndim,1:nBHnei))
  iBHnei = 0
  do k = -nrad, nrad
     z = dble(k) / dble(nrad) * rrad
     do j = -nrad, nrad
        y = dble(j) / dble(nrad) * rrad
        do i = -nrad, nrad
           x = dble(i) / dble(nrad) * rrad
           rr = sqrt(dble(x*x+y*y+z*z))
           if(rr .le. rrad)then
              iBHnei = iBHnei + 1
              xBHnei(1, iBHnei) = x
              xBHnei(2, iBHnei) = y
              xBHnei(3, iBHnei) = z
           endif
        enddo
     enddo
  enddo
  weight=1d0/dble(nBHnei)
  
  ! Open cache for array uold (fetch) and unew (flush)
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_large_realdp)/32,&
                pack=pack_fetch_refine,unpack=unpack_fetch_refine,&
                init=init_flush_godunov, flush=pack_flush_godunov,&
                combine=unpack_flush_godunov, bound=init_bound_refine)

  ! Loop over particles in Hilbert order
  hash_nbor(0) = ilevel+1
  do ipart = p%headp(ilevel), p%tailp(ilevel)
     
     ! Black hole position
     xcen(1:ndim) = p%xp(ipart,1:ndim) / dx_loc

     ! Initialise sink information at zero
     rho_sink=0d0; vel_sink=0d0; cs_sink=0d0

     ! Loop over all cells in the accretion region, collecting physical properties and hash IDs
     do j = 1,nBHnei

        ! Compute neighbouring cell coordinates
        xnei(1:ndim) = xcen(1:ndim) + xBHnei(1:ndim,j)
        ! Periodic boundary conditions
        do idim=1,ndim
           if(xnei(idim)<                0.0d0)xnei(idim)=xnei(idim)+m%ckey_max(ilevel+1)
           if(xnei(idim)>=m%ckey_max(ilevel+1))xnei(idim)=xnei(idim)-m%ckey_max(ilevel+1)
        end do
        ! Get neighboring cell at current level
        hash_nbor(1:ndim)  = int(xnei(1:ndim))
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.)
        ! If missing cycle
        if(.not.associated(gridn))cycle

        ! Get physical information
        d                = max(gridn%uold(icelln,1),r%smallr)
        vv(1)            =     gridn%uold(icelln,2)/d
        vv(2)            =     gridn%uold(icelln,3)/d
        vv(3)            =     gridn%uold(icelln,4)/d
        e                =     gridn%uold(icelln,5)
        ethermal         = (e - 0.5d0*d*sum(vv(:)**2)) / d
        cs               = max((r%gamma-1.0d0)*ethermal,r%smallc**2)!*boost**(-2d0/3d0)
        ! TODO: Add boost above
        ! Add to average (weighted) information
        rho_sink         = rho_sink + d * weight
        vel_sink(1:ndim) = vel_sink(1:ndim) + vv(1:ndim) * weight
        cs_sink          = cs_sink + cs * weight
     end do

     !write(*,*)'test A',nBHnei,d,vv(:),cs
     !!! Compute BHL accretion rate
     ! Bondi velocity
     v_rel(1:ndim) = vel_sink(1:ndim) - p%vp(ipart,1:ndim)
     if(r%bondi_use_vrel)then
        velocity   = sqrt(sum(v_rel(:)**2))
        v_bondi    = sqrt(velocity**2 + cs**2)
     else
        v_bondi    = sqrt(cs**2)
     end if
     ! Sink radius
     r2_sink = (factG * p%mp(ipart) / v_bondi**2)**2
     ! Density at infinity (using extrapolation)
     !rho_inf(isink)=density/(bondi_alpha(ir_cloud*0.5d0*dx_min/(r2(isink)+tiny(0.0_dp))**0.5d0))
     
     ! Bondi-Hoyle-Lyttleton accretion rate
     dMBH_overdt = 4.0d0 * pi * rho_sink * r2_sink * v_bondi
     ! Eddinton accretion rate
     dMED_overdt = 4.0d0*pi*factG_in_cgs*p%mp(ipart)*mH/(0.1d0*sigma_T*c_cgs)*scale_t
     if(r%eddington_cap)dMBH_overdt = min(dMBH_overdt,dMED_overdt)

     write(*,*)'Sink properties:',rho_sink,cs_sink,v_bondi,r2_sink
     write(*,*)'Accretion:',ipart,dMBH_overdt,dMED_overdt,p%mp(ipart)

     ! Loop over all cells in the accretion region, proceeding with accretion
     m_acc=0.0d0; x_acc=0.0d0; p_acc=0.0d0; l_acc=0.0d0
     do j = 1, nBHnei

        ! Compute neighbouring cell coordinates
        xnei(1:ndim) = xcen(1:ndim) + xBHnei(1:ndim,j)
        ! Periodic boundary conditions
        do idim=1,ndim
           if(xnei(idim)<                0.0d0)xnei(idim)=xnei(idim)+m%ckey_max(ilevel+1)
           if(xnei(idim)>=m%ckey_max(ilevel+1))xnei(idim)=xnei(idim)-m%ckey_max(ilevel+1)
        end do
        ! Get neighboring cell at current level
        hash_nbor(1:ndim)  = int(xnei(1:ndim))
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.)
        ! If missing cycle
        if(.not.associated(gridn))cycle

        ! Get physical information
        d          = max(gridn%uold(icelln,1),r%smallr)
        vv(1:ndim) =     gridn%uold(icelln,2:ndim+1)/d
        e          =     gridn%uold(icelln,5)/d
        ! We need to remove all non-thermal energies as they are not accreted

        ! Get accreted mass for this cell
        d_acc = dMBH_overdt * g%dtnew(p%levelp(ipart)) * weight / vol_loc

        ! Compute weighting (TODO: Add something more intelligent)
        if(r%accretion_method=='mass')then
           d_acc = d_acc * d / rho_sink
        end if

        ! TODO: Add other limiters here
        d_acc = min(d_acc, 0.25d0 * d)
        d_acc = max(d_acc, 0.0_dp)

        ! Accrete from the cell
        gridn%unew(icelln,1)         = gridn%unew(icelln,1)      - d_acc
        do idim=1,ndim
           gridn%unew(icelln,idim+1) = gridn%unew(icelln,idim+1) - d_acc * vv(idim)
        end do
        gridn%unew(icelln,5)         = gridn%unew(icelln,5)      - d_acc * e
        ! TODO: Add passive scalar accretion here

        !!! Accrete onto the black hole
        ! Accreted mass
        m_acc = m_acc + d_acc * vol_loc
        ! Accreted relative center of mass
        x_acc(1:ndim) = x_acc(1:ndim) + d_acc * xBHnei(1:ndim,j) * vol_loc
        ! Accreted relative momentum
        p_acc(1:ndim) = p_acc(1:ndim) + d_acc + v_rel(1:ndim) * vol_loc
        ! Accreted relative angular momentum
        l_acc(1:ndim) = l_acc(1:ndim) + d_acc * cross(xBHnei(1:ndim,j), v_rel(1:ndim)) * vol_loc

     end do ! End loop over j

     ! Add accreted properties to sink variables
     p%xp(ipart,1:ndim) = ( p%mp(ipart) * p%xp(ipart,1:ndim) + x_acc(1:ndim) ) / ( p%mp(ipart) + m_acc )
     p%vp(ipart,1:ndim) = ( p%mp(ipart) * p%vp(ipart,1:ndim) + p_acc(1:ndim) ) / ( p%mp(ipart) + m_acc )
     p%jp(ipart,1:ndim) = ( p%mp(ipart) * p%jp(ipart,1:ndim) + l_acc(1:ndim) ) / ( p%mp(ipart) + m_acc )
     p%mp(ipart)        = p%mp(ipart) + m_acc

  end do ! End loop over ipart

  call close_cache(s,m%grid_dict)

  deallocate(xBHnei)

  end associate
#endif
#endif
end subroutine sink_accretion

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module sink_accretion_module
