module BZ_sink_module
   use rho_fine_module, only: cic_weight, cic_index, tsc_weight, tsc_index, pcs_weight, pcs_index
contains
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine prepare_blandford_znajek_jet(s,p,ipart,f_edd,eta_rad,eta_BZ,scale_t,scale_v,edot_jet,pdot_jet,mdot_jet)
   use constants
   use amr_parameters, only: ndim,dp
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ipart
   real(dp)::f_edd,eta_rad,eta_BZ,scale_t,scale_v
   real(dp)::edot_jet,pdot_jet,mdot_jet
   !==================================================================
   ! This is the RAMSES routine to compute all necessary feedback quantities
   ! for the Blandford-Znajek jet. See Talbot+2021 for details.
   ! Written by Nicholas Choustikov (May 2025)
   !==================================================================
   real(dp)::t_salp,m_dot_internal

   ! Compute the internal accretion rate
   t_salp = (sigma_T * c_cgs / (4.0d0 * pi * factG_in_cgs * mH)) * (eta_rad/0.1d0) / scale_t
   m_dot_internal = f_edd * p%mBH(ipart) / t_salp

   ! Compute the jet mass flux
   mdot_jet = m_dot_internal * s%r%jet_mass_loading / (1.0d0 + s%r%jet_mass_loading)

   ! Compute the jet energy flux
   edot_jet = eta_BZ/(1.0d0 + s%r%jet_mass_loading) * m_dot_internal * (c_cgs/scale_v)**2

   ! Compute the jet momentum flux
   pdot_jet = sqrt(2.0d0*s%r%jet_mass_loading*eta_BZ)*(c_cgs/scale_v)/(1.0d0 + s%r%jet_mass_loading) * m_dot_internal

end subroutine prepare_blandford_znajek_jet

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine launch_blandford_znajek_jet(s,p,ipart,ilevel,edot_jet,pdot_jet,mdot_jet,tan_theta,nBH_fb_nei,dx_loc,vol_loc)
   use constants
   use amr_parameters, only: ndim,twotondim,dp
   use hydro_parameters, only: nvar, nener
   use amr_commons, only: nbor,oct
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t,cross
   use params_module
   use nbors_utils
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ilevel,ipart,nBH_fb_nei
   real(dp)::edot_jet,pdot_jet,mdot_jet,tan_theta
   real(dp)::dx_loc,vol_loc
   !==================================================================
   ! This is the RAMSES routine which launches the Blandford-Znajek jet.
   ! Written by Nicholas Choustikov (Apr 2025)
   !==================================================================
   real(dp)::rr,x,y,z,rrad
   real(dp),dimension(1:ndim,1:nBH_fb_nei)::xBH_fb_nei
   integer,dimension(1:ndim,1:nBH_fb_nei)::ckey_fb_nei
   real(dp),dimension(1:nBH_fb_nei)::weight_fb_nei
   real(dp),dimension(1:ndim)::xcen,xnei,x_rel
   integer(kind=8),dimension(0:ndim)::hash_nbor
   real(dp),dimension(1:ndim)::jet_direction
   integer::i,j,k,ii,jj,kk,icelln,ind,idim,ivar,iBHnei
   real(dp)::d,e,ethermal,r_rel,rho_gas_fb
   real(dp),dimension(1:ndim)::vv
   real(dp)::cone_dist,orth_dist,weight,local_weight,total_weight
   type(oct),pointer::gridn
   logical::ok
   real(dp)::fbk_mass_agn_loc,fbk_mom_agn_loc,fbk_ener_agn_loc
   real(dp),dimension(1:ndim,1:twotondim)::xCIC
   integer,dimension(1:ndim,1:twotondim)::ckeyCIC
   real(dp),dimension(1:twotondim)::volCIC

   associate(r=>s%r,g=>s%g,m=>s%m)

   if(r%verbose)write(*,*)'Entering launch_blandford_znajek_jet...'
   if(.not.r%agn)return

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! AGN Feedback: Set everything up
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   
   hash_nbor(0) = ilevel+1
   xBH_fb_nei=0d0; ckey_fb_nei=0d0; weight_fb_nei=0d0

   ! Black hole position
   xcen(1:ndim) = p%xp(ipart,1:ndim) / dx_loc

   ! Find the jet direction
   jet_direction(1:ndim) = p%jp(ipart,1:ndim) / (norm2(p%jp(ipart,:)) + tiny(0.0_dp)) 

   !!! Compute all of the necessary weights
   ! Loop over all possible cells within the feedback region
   iBHnei = 0
   weight_fb_nei = 0d0; xBH_fb_nei = 0d0; total_weight=0d0
   do kk=-r%agn_feedback_radius,r%agn_feedback_radius
      x_rel(3)=dble(kk)
      do jj=-r%agn_feedback_radius,r%agn_feedback_radius
         x_rel(2)=dble(jj)
         do ii=-r%agn_feedback_radius,r%agn_feedback_radius
            x_rel(1)=dble(ii)
            r_rel=norm2(x_rel(:))
            if(r_rel.lt.dble(r%agn_feedback_radius))then
               iBHnei = iBHnei + 1

               !!! Collect all of the necessary positions and cartesian keys
               do idim=1,ndim
                  ! New CIC version
                  xBH_fb_nei(idim,iBHnei) = x_rel(idim) + xcen(idim)
                  ckey_fb_nei(idim,iBHnei) = int(xBH_fb_nei(idim,iBHnei))
               end do

               !!! Compute the weight of the cell in question
               ok=.false.

               cone_dist = dot_product(x_rel(1:ndim),jet_direction(1:ndim))
               orth_dist = norm2((x_rel(1:ndim) - cone_dist*jet_direction(1:ndim)))
               if(orth_dist.le.abs(cone_dist)*tan_theta)ok=.true.
               if(r_rel.lt.1)ok=.false. ! Exclude the central cell in jet mode

               if(ok)then
                  call BZ_psy_function(r_rel,local_weight)
                  weight_fb_nei(iBHnei) = local_weight
                  total_weight = total_weight + local_weight

                  hash_nbor(1:ndim)  = ckey_fb_nei(1:ndim,iBHnei)
                  call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.)
                  ! If missing cycle
                  if(.not.associated(gridn))cycle
                  d = max(gridn%uold(icelln,1), r%smallr)
                  rho_gas_fb = rho_gas_fb + d*local_weight
               end if

            end if
         end do ! End loop over ii
      end do ! End loop over jj
   end do ! End loop over kk

   ! Normalise the weights
   if(total_weight==0.0d0)write(*,*)'PROBLEM: total weight 0'
   weight_fb_nei = weight_fb_nei / total_weight
   rho_gas_fb = rho_gas_fb / total_weight

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Administer the AGN Feedback
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Loop over the affected cells
   do iBHnei=1,nBH_fb_nei
      ! Skip cells with zero weight
      if(weight_fb_nei(iBHnei)==0.0d0)cycle
      call BZ_sink_B_spline_weights_CIC(s,xBH_fb_nei(1:ndim,iBHnei),xCIC,ckeyCIC,volCIC,ilevel)
      
      do j=1,twotondim              
         ! Compute neighbouring cell coordinates
         xnei(1:ndim) = xCIC(1:ndim,j) 
         !xnei(1:ndim) = xBH_fb_nei(1:ndim,iBHnei) ! For now we use this, as otherwise some of the cells do overlap with the central cell
         x_rel(1:ndim) = xnei(1:ndim) - xcen(1:ndim)

         ! Periodic boundary conditions
         do idim=1,ndim
            ! Note, periodic BCs for xCIC are already enforced in sink_B_spline_weights_CIC
            if(x_rel(idim)<-r%boxlen/2d0)x_rel(idim)=x_rel(idim)+r%boxlen
            if(x_rel(idim)> r%boxlen/2d0)x_rel(idim)=x_rel(idim)-r%boxlen
         end do
         r_rel = norm2(x_rel(:))

         ! Get neighboring cell at current level
         hash_nbor(1:ndim)  = ckeyCIC(1:ndim,j)
         call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.)
         ! If missing cycle
         if(.not.associated(gridn))cycle

         ! Get the local gas properties
         d          = max(gridn%uold(icelln,1),r%smallr)
         vv(1:ndim) =     gridn%uold(icelln,2:ndim+1)/d
         !e          =     gridn%uold(icelln,5)/d

         ! Get the weight
         weight = volCIC(j) * weight_fb_nei(iBHnei)

         !!! Do the feedback in this cell
         ! Get the local feedback properties
         fbk_mass_agn_loc = mdot_jet * g%dtnew(ilevel) * weight / vol_loc
         fbk_mom_agn_loc  = pdot_jet * g%dtnew(ilevel) * weight / vol_loc
         if(.not.r%BZ_momentum_conserving_jet)fbk_ener_agn_loc = edot_jet*g%dtnew(ilevel)*weight*(d+fbk_mass_agn_loc)/vol_loc
   
         ! Mass-weigh if required
         if(r%agn_use_mass_weighting)then
            fbk_mass_agn_loc = fbk_mass_agn_loc * (d/rho_gas_fb)
            fbk_mom_agn_loc  = fbk_mom_agn_loc  * (d/rho_gas_fb)
            if(.not.r%BZ_momentum_conserving_jet)fbk_ener_agn_loc=fbk_ener_agn_loc*(d/rho_gas_fb)
         end if

         ! Now we inject the actual feedback
         gridn%unew(icelln,1)     = gridn%unew(icelln,1)   + fbk_mass_agn_loc
         gridn%unew(icelln,2:4)   = gridn%unew(icelln,2:4) + fbk_mom_agn_loc*dot_product(jet_direction(:),x_rel(:))*jet_direction(1:ndim)/(r_rel+tiny(0.0_dp))
         if(r%BZ_momentum_conserving_jet)then
            gridn%unew(icelln,5)  = gridn%unew(icelln,5)   + fbk_mom_agn_loc*dot_product(jet_direction(:),x_rel(:)/(r_rel+tiny(0.0_dp)))*dot_product(jet_direction(1:ndim), vv(1:ndim))
         else
            gridn%unew(icelln,5)  = gridn%unew(icelln,5)   + fbk_ener_agn_loc
         end if
      end do ! End loop over j
   end do ! End loop over iBHnei

   end associate

end subroutine launch_blandford_znajek_jet

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine evolve_BH_disc_system(s,p,ipart,ilevel,factG,scale_d,scale_l,scale_t,m_acc,l_acc,f_edd,eta_rad,eta_BZ)
   use constants
   use amr_parameters, only: ndim,dp
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t

   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ipart,ilevel
   real(dp)::factG,scale_d,scale_l,scale_t,scale_v
   real(dp)::f_edd,eta_rad,eta_BZ
   real(dp),dimension(1:ndim)::l_acc
   !==================================================================
   ! This is the RAMSES routine to drive the evolution of the internal
   ! Black hole/disk sub-grid model. This is designed to be modular,
   ! allowing easy coupling to any desired disk model, as long as the
   ! structure is followed.
   ! Written by Nicholas Choustikov (May 2025)
   !==================================================================
   real(dp)::a,J_BH_mag,J_Disc_mag,grade
   real(dp)::scale_m_msun,m_acc,m_acc_interior

   associate(r=>s%r,g=>s%g,m=>s%m)

   if(r%verbose)write(*,*)'Entering evolve_BH_disk_system...'

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Compute necessary quantities from last timestep
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Get the black hole spin parameter
   ! TODO: Figure out exactly when this should be pro/retro-grade
   J_BH_mag = norm2(p%jp(ipart,1:ndim))
   
   a = (c_cgs*scale_l/scale_t) * J_BH_mag / (factG * p%mBH(ipart)**2)
   !a = min(a,0.998d0) ! Limit spin parameter Thorne+1974

   ! Compute if the accretion is pro- or retro-grade
   grade = sign(1.0d0, dot_product(p%jp(ipart,1:ndim),p%jD(ipart,1:ndim)))

   ! Useful units
   scale_m_msun = scale_d * scale_l**3 / m_sun
   scale_v = scale_l/scale_t

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Use prescribed disk structure to solve for accretion rate
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   call solve_for_internal_accretion_rate(s,p,ipart,a,grade,scale_m_msun,f_edd)

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Compute all efficiencies
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   call radiative_efficiency(a,grade,f_edd,r%fedd_ADAF,r%fedd_Edd,eta_rad)
   call BZ_efficiency(a,f_edd,r%fedd_ADAF,eta_BZ)

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Compute mass evolution of BH/Disc system
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   call drive_mass_evolution(s,p,ipart,ilevel,f_edd,a,m_acc,eta_rad,eta_BZ,scale_t,scale_m_msun,m_acc_interior)

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Compute angular momentum evolution for BH/Disc system
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   call drive_angular_momentum_evolution(s,p,ipart,ilevel,factG,scale_v,f_edd,a,grade,m_acc,m_acc_interior,l_acc,eta_rad,eta_BZ,scale_t,scale_m_msun)

   end associate
   
end subroutine evolve_BH_disc_system

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine drive_mass_evolution(s,p,ipart,ilevel,f_edd,a,m_acc,eta_rad,eta_BZ,scale_t,scale_m_msun,m_acc_interior)
   use constants
   use amr_parameters, only: dp
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ipart,ilevel
   real(dp)::f_edd,a,m_acc,eta_rad,eta_BZ
   real(dp)::scale_t,scale_m_msun,scale_v,m_acc_interior
   !==================================================================
   ! Update the black hole and disc masses based on the internal evolution
   ! of the disc model. Here, we use a second-order scheme (Fiacconi+18;Talbot+21).
   ! Written by Nicholas Choustikov (May 2025)
   !==================================================================
   real(dp)::m_disc_temp,m_bh_temp
   real(dp)::t_salp,self_gravity_mass

   ! Compute the Salpeter time
   t_salp       = (sigma_T * c_cgs / (4.0d0 * pi * factG_in_cgs * mH)) * (eta_rad/0.1) / scale_t
   
   ! Internal mass accretion
   m_acc_interior = (f_edd / t_salp) * p%mBH(ipart) * s%g%dtnew(ilevel)

   ! Do the predicter step of the integration
   m_bh_temp    = p%mBH(ipart) + ((1.0d0 - eta_rad - eta_BZ)/(1 + s%r%jet_mass_loading)) * m_acc_interior
   m_disc_temp  = p%mD(ipart)  - m_acc_interior + m_acc
   
   ! Do the corrector step of the integration
   p%mBH(ipart) = p%mBH(ipart) + ((1.0d0 - eta_rad - eta_BZ)/(1 + s%r%jet_mass_loading)) * (f_edd / t_salp) * 0.5d0 * (p%mBH(ipart) + m_bh_temp) * s%g%dtnew(ilevel)
   p%mD(ipart)  = p%mD(ipart)  - (f_edd/t_salp)*0.5d0*(p%mBH(ipart) + m_bh_temp)*s%g%dtnew(ilevel)

   ! Update the dynamical mass of the particle (i.e. BH + Disc mass)
   p%mp(ipart)  = p%mBH(ipart) + p%mD(ipart)

end subroutine drive_mass_evolution

subroutine drive_angular_momentum_evolution(s,p,ipart,ilevel,factG,scale_v,f_edd,a,grade,m_acc,m_acc_interior,l_acc,eta_rad,eta_BZ,scale_t,scale_m_msun)
   use constants
   use amr_parameters, only: dp,ndim
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t,cross
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ipart,ilevel
   real(dp)::factG,f_edd,a,grade,m_acc,eta_rad,eta_BZ
   real(dp)::scale_t,scale_v,scale_m_msun,m_acc_interior
   real(dp),dimension(1:ndim)::l_acc
   !==================================================================
   ! Update the black hole and disc angular momenta based on the internal
   ! model of the disc. Here, we account for both the diffusive and wave-like
   ! regimes (Ingram & Motta 2019). See also Fiacconi+18;Talbot+21;Koudmani+24;Kao+25.
   ! Written by Nicholas Choustikov (May 2025)
   !==================================================================
   real(dp)::L_ISCO,L_BZ,L,phi_BH,omega_BH
   real(dp)::J_BH_mag,J_BH_mag_new
   real(dp),dimension(1:ndim)::J_BH_temp,J_D_temp,J_BH_hat,J_cons_temp,J_BH_temp_hat,J_D_temp_hat
   real(dp)::eta_rad_thin,m_bh_warp,f_edd_crit
   real(dp)::t_GM,t_acc,omega_prec
   logical::king_check
   real(dp)::r_in,r_trap,r_isco,r_bending_wave,r_s,J_trap

   ! Compute the specific angular momentum at ISCO (Bardeen+1972) 
   L = Lam(a,grade)
   if(grade>=0)then
      L_ISCO = (factG * p%mBH(ipart) / (c_cgs/scale_v)) / L * ((L**2 - 2.0d0*a*sqrt(L) + a**2) / (L - 3.0d0 + 2.0d0*a/sqrt(L)))
   else
      L_ISCO = -(factG * p%mBH(ipart) / (c_cgs/scale_v)) / L * ((L**2 + 2.0d0*a*sqrt(L) + a**2) / (L - 3.0d0 - 2.0d0*a/sqrt(L)))
   end if

   ! Compute the specific angular momentum of the Blandford-Znajek jet
   phi_BH   = 52.6d0 - 20.2d0*a**3 - 14.9d0*a**2 + 34.0d0*a ! See Narayan+2022
   omega_BH = (a / 2.0d0) / (1.0d0 + sqrt(1.0d0 - a**2))    ! See Tchekhovskoy+2012
   L_BZ     = 1.0d0/(12.0d0*pi**2) * phi_BH**2 * omega_BH * (factG * p%mBH(ipart) / (c_cgs/scale_v)) ! See Talbot+2021

   ! Evolve the first step (ISM->Disc accretion + Disc->ISCO accretion + BH->Jet)
   J_BH_mag            = norm2(p%jp(ipart,1:ndim))
   J_BH_hat(1:ndim)    = p%jp(ipart,1:ndim) / (J_BH_mag + tiny(0.0_dp)) 

   J_D_temp(1:ndim)    = p%jD(ipart,1:ndim) + l_acc(1:ndim) - (s%r%jet_mass_loading/(1.0d0 + s%r%jet_mass_loading))*L_ISCO*m_acc_interior*J_BH_hat(1:ndim)
   J_BH_temp(1:ndim)   = p%jp(ipart,1:ndim) - (1.0d0/(1.0d0 + s%r%jet_mass_loading)) * L_BZ * m_acc_interior * J_BH_hat(1:ndim)

   ! Compute the (temporary) conserved angular momentum
   J_cons_temp(1:ndim) = J_BH_temp(1:ndim) + J_D_temp(1:ndim)

   ! Evolve the second step (ISCO accretion -> BH)
   J_BH_mag_new = norm2(J_BH_temp(1:ndim)) + (1.0d0/(1.0d0 + s%r%jet_mass_loading)) * L_ISCO * m_acc_interior

   ! Cap the new spin at 0.998 (Thorne+1974)
   J_BH_mag_new = min(J_BH_mag_new, factG * p%mBH(ipart)**2 * c_cgs/scale_v * s%r%BH_spin_max)

   ! Compute the warping mass of the Black Hole (Fiacconi+18;Talbot+21;Kao+25)
   eta_rad_thin = radiative_efficiency_thin(a,grade)
   m_bh_warp = (1.0D7/scale_m_msun) * (p%mD(ipart) * scale_m_msun / 1d4)**(35/82) * (s%r%disc_viscosity/0.1d0)**(-1/41) * (f_edd/(eta_rad_thin/0.1))**(-17/82) * a**(-25/82)

   ! Compute the critical accretion rate between the Diffusive (i.e. Bardeen-Peterson) and wave regimes (see Ingram & Motta 2019) following Kao+2025
   ! Here, we define it such that r_pt = r_warp. Also, we assume zeta = 1 (otherwise there's a factor of zeta^(-20/41))
   f_edd_crit = 16.0d0 * (s%r%disc_viscosity_zeta)**(-20/41) * (s%r%disc_viscosity/0.1d0)**(24/41) * (eta_rad/0.1) * (p%mBH(ipart) * scale_m_msun / 1d6)**(-4/41) * a**(20/41)
   if(s%r%disc_model_type==1)f_edd_crit=huge(0.0_dp) ! We should never enter this regime for pure-thin-disc models

   ! Find the temporary Black Hole and Disc unit vectors
   J_BH_temp_hat(1:ndim) = J_BH_temp(1:ndim) / (norm2(J_BH_temp(1:ndim)) + tiny(0.0_dp)) 
   J_D_temp_hat(1:ndim)  = J_D_temp(1:ndim)  / (norm2(J_D_temp(1:ndim))  + tiny(0.0_dp)) 

   ! Evolve the direction of the BH spin (Lense-Thirring etc.)
   if(f_edd>=f_edd_crit)then
      !!! We are in the wave regime, follow Ingram & Motta 2019 (see also Koudmani+24,Kao+25)
      ! Compute all needed radii (all normalised to r_s)
      r_isco = Lam(a,grade)/2.0d0
      r_trap = 15.0d0 * (f_edd / (eta_rad/0.1d0))
      r_bending_wave = 3.0d0 * a**(2/5) ! Note, we assume H/R = 1
      r_in = max(r_isco,r_bending_wave)
      r_s = 2.0d0*factG*p%mBH(ipart)/(c_cgs/scale_v)**2

      ! Compute the precession frequency
      omega_prec = (2.0d0*(c_cgs/scale_v)*a/r_s) * ((0.51d0+2.0d0)/(0.51d0-1.0d0)) * ((r_trap)**(0.51d0-1.0d0) - (r_in)**(0.51d0-1.0d0))/((r_trap)**(0.51d0+2.0d0) - (r_in)**(0.51d0+2.0d0))

      ! Compute the angular momentum contained within the photon trapping region
      !call compute_disc_angular_momentum()
      J_trap = norm2(p%jD(ipart,1:ndim))

      ! Compute the accretion timescale
      t_acc = (sigma_T * c_cgs / (4.0d0 * pi * factG_in_cgs * mH)) * (eta_rad/0.1) / scale_t / f_edd

      ! Update the Black Hole spin direction
      p%jp(ipart,1:ndim) = p%jp(ipart,1:ndim) - J_BH_mag_new*((J_trap/(J_BH_mag_new+tiny(0.0_dp)) )*omega_prec*cross(J_BH_temp_hat,J_D_temp_hat) + (2.0d0*pi/t_acc)*cross(J_BH_temp_hat,cross(J_BH_temp_hat,J_D_temp_hat)))
   else
      ! We are in the diffusive, first we check if the BH exceeds the warp mass
      if(p%mBH(ipart)>=m_bh_warp)then
         !!! We align instantly, following King+2005
         ! The Black Hole aligns along the total angular momentum
         p%jp(ipart,1:ndim) = J_BH_mag_new/(norm2(J_cons_temp)+tiny(0.0_dp))  * J_cons_temp(1:ndim)

         ! Check how the disc should align, following King+2005
         king_check = (dot_product(J_BH_temp_hat,J_D_temp_hat)>=-norm2(J_D_temp)/(2.0d0*(J_BH_mag_new + tiny(0.0_dp)) ))

         ! Align the disc as required
         if(king_check)then
            p%jD(ipart,1:ndim) = norm2(J_D_temp)/(norm2(J_cons_temp)+tiny(0.0_dp))  * J_cons_temp(1:ndim)
         else
            p%jD(ipart,1:ndim) = -norm2(J_D_temp)/(norm2(J_cons_temp)+tiny(0.0_dp))  * J_cons_temp(1:ndim)
         end if
      else
         !!! We are in the Bardeen-Peterson 1975 configuration
         ! Compute the gravito-magnetic timescale (Martin+2007;Perego+2009;Dotti+2013)
         t_GM = 1.7d5 * (p%mBH(ipart) * scale_m_msun / 1d6)**(-2/35) * (f_edd/(eta_rad_thin/0.1))**(-32/35) * (s%r%disc_viscosity/0.1d0)**(58/35) * a**(5/7) / scale_t

         ! Update the Black Hole spin direction
         p%jp(ipart,1:ndim) = p%jp(ipart,1:ndim) - (s%g%dtnew(ilevel)/t_GM) * J_BH_mag_new * (sin(pi/7.0d0)*cross(J_BH_temp_hat,J_D_temp_hat) + cos(pi/7.0d0)*cross(J_BH_temp_hat,cross(J_BH_temp_hat,J_D_temp_hat)))
      end if
   end if

   ! Evolve the disc angular momentum using the conserved angular momentum
   p%jD(ipart,1:ndim) = J_cons_temp(1:ndim) - p%jp(ipart,1:ndim)

end subroutine drive_angular_momentum_evolution

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine solve_for_internal_accretion_rate(s,p,ipart,a,grade,scale_m_msun,f_edd)
   use amr_parameters, only: dp
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ipart
   real(dp)::a,grade,scale_m_msun,f_edd
   !==================================================================
   ! This is the RAMSES routine to compute the internal accretion rate
   ! of the sub-grid disc based on the given (multi-)zone model.
   ! Written by Nicholas Choustikov (May 2025)
   !==================================================================

   ! Call the appropriate routine
   if(s%r%disc_model_type==1)then
      call thin_disk_one_zone_f_edd(s,p,ipart,a,grade,scale_m_msun,f_edd)
   else
      write(*,*)'Unknown internal accretion disk model used...'
   end if

   ! If desired, limit this accretion rate
   f_edd = min(f_edd, s%r%max_internal_accretion_rate)

end subroutine solve_for_internal_accretion_rate

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine thin_disk_one_zone_f_edd(s,p,ipart,a,grade,scale_m_msun,f_edd)
   use amr_parameters, only: dp,ndim
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ipart
   real(dp)::a,grade,scale_m_msun,f_edd
   !==================================================================
   ! Accretion rate computation assuming outer solution (region c) of 
   ! Shakura & Sunyaev 1973 solution.
   ! Written by Nicholas Choustikov (May 2025)
   !==================================================================
   real(dp)::J_BH_mag,J_D_mag

   J_BH_mag = norm2(p%jp(ipart,1:ndim))
   J_D_mag = norm2(p%jD(ipart,1:ndim))
   f_edd = 0.76d0 * (radiative_efficiency_thin(a,grade)/0.1d0) * (s%r%disc_viscosity/0.1d0)**(8/7) * (p%mD(ipart) * scale_m_msun / 1d4)**5 * (p%mBH(ipart) * scale_m_msun / 1d6)**(-47/7) * (a*J_D_mag*sign(1.0d0,dot_product(p%jp(ipart,1:ndim),p%jD(ipart,1:ndim)))/(3.0d0 * J_BH_mag))**(-25/7)

end subroutine thin_disk_one_zone_f_edd

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine compute_self_gravity_mass(s,p,ipart,a,scale_m_msun,f_edd,eta_rad,m_sg)
   use amr_parameters, only: dp
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ipart
   real(dp)::a,scale_m_msun,f_edd,eta_rad,m_sg
   !==================================================================
   ! This is the RAMSES routine to compute the self-gravity mass of the
   ! sub-grid disc based on the given (multi-)zone model.
   ! Written by Nicholas Choustikov (May 2025)
   !==================================================================

   if(f_edd>=0)then
      if(s%r%disc_model_type==1)then
         call thin_disk_one_zone_sg_mass(s,p,ipart,a,scale_m_msun,f_edd,eta_rad,m_sg)
      else
         write(*,*)'Unknown internal accretion disk model used...'
      end if
   else
      if(s%r%disc_model_type==1)then
         call thin_disk_one_zone_sg_mass_no_fedd(s,p,ipart,scale_m_msun,m_sg)
      else
         write(*,*)'Unknown internal accretion disk model used...'
      end if
   end if
end subroutine compute_self_gravity_mass

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine thin_disk_one_zone_sg_mass(s,p,ipart,a,scale_m_msun,f_edd,eta_rad,m_sg)
   use amr_parameters, only: dp
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ipart
   real(dp)::a,scale_m_msun,f_edd,eta_rad,m_sg
   !==================================================================
   ! Disc self-gravity mass computation assuming outer solution (region c) of 
   ! Shakura & Sunyaev 1973 solution.
   ! Written by Nicholas Choustikov (May 2025)
   !==================================================================
   
   m_sg = 2.0d4 * (s%r%disc_viscosity/0.1)**(-1/45) * (p%mBH(ipart) * scale_m_msun / 1d6)**(34/45) * (f_edd/(eta_rad/0.1))**(4/45)
   
end subroutine thin_disk_one_zone_sg_mass

subroutine thin_disk_one_zone_sg_mass_no_fedd(s,p,ipart,scale_m_msun,m_sg)
   use constants, only: c_cgs,factG_in_cgs,m_sun
   use amr_parameters, only: dp,ndim
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ipart
   real(dp)::scale_m_msun,m_sg
   !==================================================================
   ! Disc self-gravity mass computation assuming outer solution (region c) of 
   ! Shakura & Sunyaev 1973 solution. Here, we take the no-fedd approach.
   ! Written by Nicholas Choustikov (May 2025)
   !==================================================================
   
   m_sg = 19518.0d0 * (s%r%disc_viscosity/0.1)**(-5/63) * (p%mD(ipart) * scale_m_msun / 1d4)**(4/45) * (p%mBH(ipart) * scale_m_msun / 1d6)**(19/105) * (c_cgs * norm2(p%jD(ipart,1:ndim)) / (3 * factG_in_cgs * (p%mBH(ipart) * scale_m_msun*m_sun)**2))
   
end subroutine thin_disk_one_zone_sg_mass_no_fedd

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine radiative_efficiency(a,grade,f_edd,fedd_ADAF,f_edd_Edd,eta_rad)
   use amr_parameters, only:dp
   implicit none
   real(dp)::a,grade,f_edd,fedd_ADAF,f_edd_Edd
   real(dp)::eta_rad
   real(dp)::A_var,B_var,C_var
   real(dp)::transition_lower,transition_higher
   ! Compute the accretion-dependent radiative efficiency
   if(f_edd>=f_edd_Edd)then
      ! Madau+2014, Sadowski+2009
      if(grade>=0)then
         A_var   = (0.9663 - 0.9292*a)**(-0.5639)
         B_var   = (4.627  -  4.445*a)**(-0.5524)
         C_var   = (827.3  -  718.1*a)**(-0.706)
      else
         A_var   = (0.9663 + 0.9292*a)**(-0.5639)
         B_var   = (4.627  +  4.445*a)**(-0.5524)
         C_var   = (827.3  +  718.1*a)**(-0.706)
      end if
      eta_rad = 0.1d0/f_edd * A_var*(0.985d0/(1.6d0/f_edd + B_var) + 0.015d0/(1.6d0/f_edd + C_var))
   else
      transition_higher = (1.0d0 + (1.88d0/f_edd)**3)**(-1)
      transition_lower  = 1.0d0 - 1.0d0 * (1.0d0 + (fedd_ADAF/f_edd)**3)**(-1)
      eta_rad = radiative_efficiency_thin(a,grade) * (1.0d0 - transition_higher - transition_lower)
   end if
end subroutine radiative_efficiency

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine BZ_efficiency(a,f_edd,fedd_ADAF,eta_BZ)
   use constants, only: pi
   use amr_parameters, only:dp
   implicit none
   real(dp)::a,f_edd,fedd_ADAF
   real(dp)::eta_BZ
   real(dp)::omega_BH,phi_BH
   real(dp)::transition
   ! Compute the accretion-dependent BZ jet efficiency
   omega_BH = (a / 2.0d0) / (1.0d0 + sqrt(1.0d0 - a**2)) ! See Tchekhovskoy+2012
   phi_BH = 52.6d0 - 20.2d0*a**3 - 14.9d0*a**2 + 34.0d0*a ! See Narayan+2022
   transition = 1.0d0 - (1.0d0 + (fedd_ADAF/f_edd)**6)**(-1)
   eta_BZ = ((1.0d0 / (24.0d0 * pi**2)) * phi_BH**2 * omega_BH**2 * (1.0d0 + 1.38d0*omega_BH**2 - 9.2d0*omega_BH**4)) * (((f_edd/1.88d0)**(1.29d0) / (1.0d0 + (f_edd/1.88d0)**(1.29d0)))**2 + transition)
end subroutine BZ_efficiency

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine BZ_sink_B_spline_weights_CIC(s,x,xnei,ckey,vol,ilevel)
   use amr_parameters, only: ndim,dp,twotondim
   use ramses_commons, only: ramses_t
   implicit none
   type(ramses_t)::s
   real(dp),dimension(1:ndim)::x
   real(dp),dimension(1:ndim,1:twotondim)::xnei
   integer,dimension(1:ndim,1:twotondim)::ckey
   real(dp),dimension(1:twotondim)::vol
   integer::ilevel
   !==================================================================
   ! Simple routine to compute B-spline cells and weights for a given sink
   ! Nicholas Choustikov
   !==================================================================
   integer::idim,j
   integer,dimension(1:ndim)::ir,il
   real(dp),dimension(1:ndim)::dr,dl

   associate(r=>s%r,m=>s%m)

   ! CIC at level ilevel (dr: right cloud boundary; dl: left cloud boundary)
   do idim=1,ndim
      dr(idim)=x(idim)+0.5D0
      ir(idim)=int(dr(idim))
      dr(idim)=dr(idim)-ir(idim)
      dl(idim)=1.0D0-dr(idim)
      il(idim)=ir(idim)-1
   end do

   ! Periodic boundary conditions
   do idim=1,ndim
      if(il(idim)<0)il(idim)=m%ckey_max(ilevel+1)-1
      if(ir(idim)==m%ckey_max(ilevel+1))ir(idim)=0
   enddo

   ! Compute cloud volumes
   vol = cic_weight(dl,dr)

   ! Compute cartesian keys
   ckey = cic_index(il,ir)

   ! Compute neighbour positions
   do j = 1,twotondim
      do idim = 1,ndim
         xnei(idim,j) = dble(ckey(idim,j)) + 0.5d0
      end do
   end do

   end associate

end subroutine BZ_sink_B_spline_weights_CIC

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine BZ_psy_function(r,psy)
   use amr_commons
   implicit none

   real(dp)::r,psy
   logical::mode

   if(mode)then
      !!! Quasar mode
      psy = 1.0d0
   else
      !!! Radio mode
      psy = 1.0d0
   end if

end subroutine BZ_psy_function

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine BZ_dump_sink_data_fine_AGN(s,p,ipart,ilevel,scale_l,scale_t,scale_d,dMBH_overdt,dMEd_overdt,m_acc,f_edd,eta_rad,eta_BZ,edot_jet,pdot_jet,mdot_jet)
   use constants
   use amr_parameters, only: dp,ndim
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   use mdl_module, only: mdl_mkdir
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ilevel,ipart
   real(dp)::dMBH_overdt,dMEd_overdt,m_acc
   real(dp)::f_edd,eta_rad,eta_BZ
   real(dp)::edot_jet,pdot_jet,mdot_jet
   real(dp)::scale_l,scale_t,scale_d
   !==================================================================
   ! Simple routine to dump sink data to a CSV on every fine time step
   ! Nicholas Choustikov
   !==================================================================
   character(LEN=80)::filename
   integer::id_sink_loc,unit
   character(LEN=5)::nchar
   logical::file_exist
   real(dp)::a

   associate(r=>s%r, g=>s%g, mdl=>s%mdl)

   if(r%verbose_sink)write(*,*)'Entering output_sink_csv'

   ! Computing extra units and constants
   a = (c_cgs*scale_t/scale_l) * norm2(p%jp(ipart,1:ndim)) / (1.0d0 * p%mBH(ipart)**2)
   
   ! Check if the SINK file exists
   filename=TRIM('SINK')
   inquire(file=filename, exist=file_exist)
   if(.not.file_exist)call mdl_mkdir(mdl,filename)

   ! Get the filename for this sink
   call title(p%idp(ipart),nchar)
   filename=TRIM('SINK/sink_'//TRIM(nchar)//'.csv')

   ! If this is a new sink, then we need to make the file associated with that sink
   inquire(file=filename, exist=file_exist)
   unit = 10 !+g%myid
   if(.not.file_exist)then
      if(r%verbose_sink)write(*,*)'Creating file: ',filename
      open(unit=unit,file=filename,form='formatted')
      write(unit,*)'nstep,time,dt,mass,MBH,MD,dMBH,dMEd,m_acc,f_edd,eta_rad,eta_BZ,edot_jet,pdot_jet,mdot_jet,'
      close(unit)
   end if
   
   ! Open the sink file
   open(unit=unit,file=filename,form='formatted',status='unknown',position='append')

   ! Write data to the sink file
   write(unit,'(I10,21(A1,ES21.10),A1,I10)')g%nstep,',',g%t,',',g%dtnew(ilevel),',',p%mp(ipart),',',p%mBH(ipart),',',p%mD(ipart),',',dMBH_overdt,',',dMEd_overdt,',',m_acc,',',f_edd,',',eta_rad,',',eta_BZ,',',edot_jet,',',pdot_jet,',',mdot_jet
   ! Close the sink file
   close(unit)

   !42 format((8,1x,f23.15,1x),42(e23.15,1x))
   end associate
end subroutine BZ_dump_sink_data_fine_AGN

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

function Lam(a,grade)
   use amr_parameters, only:dp
   real(dp)::a,grade,Lam
   real(dp)::Z1,Z2
   ! Compute the ISCO radius function based on black hole spin
   Z1 = 1.0d0 + (1.0d0 - a**2)**(1.0d0/3.0d0) * ((1.0d0 + a)**(1.0d0/3.0d0) + (1.0d0 - a)**(1.0d0/3.0d0))
   Z2 = sqrt(3*a**2 + Z1**2)
   if(grade>=0)then
      Lam = 3.0d0 + Z2 - sqrt((3.0d0 - Z1)*(3.0d0 + Z1 + 2.0d0*Z2))
   else
      Lam = 3.0d0 + Z2 + sqrt((3.0d0 - Z1)*(3.0d0 + Z1 + 2.0d0*Z2))
   end if
end function Lam

function radiative_efficiency_thin(a,grade)
   use amr_parameters, only:dp
   real(dp)::a,grade,radiative_efficiency_thin
   ! Compute the radiative effiency around a spinning black hole
   radiative_efficiency_thin = 1.0d0 - sqrt(1.0d0 - 2.0d0/(3.0d0*Lam(a,grade)))
end function radiative_efficiency_thin
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module BZ_sink_module