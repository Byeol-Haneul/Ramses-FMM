module BZ_sink_module

contains
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine evolve_BH_disk_system(s,p,ipart,ilevel,factG,scale_d,scale_l,scale_t,m_acc)
   use constants
   use amr_parameters, only: ndim,dp
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t

   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ipart,ilevel
   real(dp)::factG,scale_d,scale_l,scale_t
   !==================================================================
   ! This is the RAMSES routine to drive the evolution of the internal
   ! Black hole/disk sub-grid model. This is designed to be modular,
   ! allowing easy coupling to any desired disk model, as long as the
   ! structure is followed.
   ! Written by Nicholas Choustikov (May 2025)
   !==================================================================
   real(dp)::a,J_BH_mag,J_Disc_mag
   real(dp)::f_edd,eta_rad,eta_BZ
   real(dp)::surplus_mass,scale_m_msun,m_acc

   associate(r=>s%r,g=>s%g,m=>s%m)

   if(r%verbose)write(*,*)'Entering evolve_BH_disk_system...'

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Compute necessary quantities from last timestep
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Get the black hole spin parameter
   ! TODO: Figure out exactly when this should be pro/retro-grade
   J_BH_mag = norm2(p%jp(ipart,1:ndim))
   
   a = (c_cgs*scale_l/scale_t) * J_BH_mag / (factG * p%mBH(ipart)**2)
   a = min(a,0.998d0) ! Limit spin parameter Thorne+1974
   a = a * sign(1.0d0, dot_product(p%jp(ipart,1:ndim),p%jD(ipart,1:ndim)))

   scale_m_msun = scale_d * scale_l**3 / m_sun

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Use prescribed disk structure to solve for accretion rate
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   call solve_for_internal_accretion_rate(s,p,ipart,a,scale_m_msun,f_edd)

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Compute all efficiencies
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   call radiative_efficiency(a,f_edd,r%fedd_ADAF,r%fedd_Edd,eta_rad)
   call BZ_efficiency(a,f_edd,r%fedd_ADAF,eta_BZ)

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Compute mass evolution of BH/Disc system
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   call drive_mass_evolution(s,p,ipart,ilevel,f_edd,a,m_acc,eta_rad,eta_BZ,scale_t,scale_m_msun,surplus_mass)

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Compute angular momentum evolution for BH/Disc system
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   !call drive_angular_momentum_evolution()

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! While we are here, compute some base feedback quantities?
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   !call 

   end associate
   
end subroutine evolve_BH_disk_system

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine solve_for_internal_accretion_rate(s,p,ipart,a,scale_m_msun,f_edd)
   use amr_parameters, only: dp
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ipart
   real(dp)::a,scale_m_msun,f_edd
   !==================================================================
   ! This is the RAMSES routine to compute the internal accretion rate
   ! of the sub-grid disc based on the given (multi-)zone model.
   ! Written by Nicholas Choustikov (May 2025)
   !==================================================================

   if(s%r%disc_model_type==1)then
      call thin_disk_one_zone_f_edd(s,p,ipart,a,scale_m_msun,f_edd)
   else
      write(*,*)'Unknown internal accretion disk model used...'
   end if

end subroutine solve_for_internal_accretion_rate

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine thin_disk_one_zone_f_edd(s,p,ipart,a,scale_m_msun,f_edd)
   use amr_parameters, only: dp,ndim
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ipart
   real(dp)::a,scale_m_msun,f_edd
   !==================================================================
   ! Accretion rate computation assuming outer solution (region c) of 
   ! Shakura & Sunyaev 1973 solution.
   ! Written by Nicholas Choustikov (May 2025)
   !==================================================================
   real(dp)::J_BH_mag,J_D_mag

   J_BH_mag = norm2(p%jp(ipart,1:ndim))
   J_D_mag = norm2(p%jD(ipart,1:ndim))
   f_edd = 0.76d0 * (radiative_efficiency_thin(a)/0.1d0) * (s%r%disc_viscosity/0.1d0)**(8/7) * (p%mD(ipart) * scale_m_msun / 1d4)**5 * (p%mBH(ipart) * scale_m_msun / 1d6)**(-47/7) * (a*J_D_mag*sign(1.0d0,dot_product(p%jp(ipart,1:ndim),p%jD(ipart,1:ndim)))/(3.0d0 * J_BH_mag))**(-25/7)

end subroutine thin_disk_one_zone_f_edd

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine drive_mass_evolution(s,p,ipart,ilevel,f_edd,a,m_acc,eta_rad,eta_BZ,scale_t,scale_m_msun,surplus_mass)
   use constants
   use amr_parameters, only: dp
   use ramses_commons, only: ramses_t
   use pm_commons, only: part_t
   implicit none
   type(ramses_t)::s
   type(part_t)::p
   integer::ipart,ilevel
   real(dp)::f_edd,a,m_acc,eta_rad,eta_BZ,surplus_mass
   real(dp)::scale_t,scale_m_msun
   !==================================================================
   ! Update the black hole and disc masses based on the internal evolution
   ! of the disc model. Here, we use a second-order scheme (Fiacconi+18;Talbot+21).
   ! Written by Nicholas Choustikov (May 2025)
   !==================================================================
   real(dp)::m_disc_temp,m_bh_temp
   real(dp)::t_salp,self_gravity_mass

   ! Compute the Salpeter time
   t_salp       = (sigma_T * c_cgs / (4.0d0 * pi * factG_in_cgs * mH)) * (eta_rad/0.1) / scale_t
   
   ! Do the predicter step of the integration
   m_bh_temp    = p%mBH(ipart) + ((1.0d0 - eta_rad - eta_BZ)/(1 + s%r%jet_mass_loading)) * (f_edd / t_salp) * p%mBH(ipart) * s%g%dtnew(ilevel)
   m_disc_temp  = p%mD(ipart)  - (f_edd/t_salp)*p%mBH(ipart)*s%g%dtnew(ilevel) + m_acc
   
   ! Do the corrector step of the integration
   p%mBH(ipart) = p%mBH(ipart) + ((1.0d0 - eta_rad - eta_BZ)/(1 + s%r%jet_mass_loading)) * (f_edd / t_salp) * 0.5d0 * (p%mBH(ipart) + m_bh_temp) * s%g%dtnew(ilevel)
   p%mD(ipart)  = p%mD(ipart)  - (f_edd/t_salp)*0.5d0*(p%mBH(ipart) + m_bh_temp)*s%g%dtnew(ilevel)

   ! Compute the self-gravity mass of the disc
   call compute_self_gravity_mass(s,p,ipart,a,scale_m_msun,f_edd,eta_rad,self_gravity_mass)

   ! Limit the disc to the self-gravity mass
   ! TODO: Should this happen during the accretion step to save from returning mass during feedback?
   surplus_mass  = max(p%mD(ipart)-self_gravity_mass, 0.0d0)
   p%mD(ipart)  = min(p%mD(ipart),self_gravity_mass)

   ! Update the dynamical mass of the particle (i.e. BH + Disc mass)
   p%mp(ipart)  = p%mBH(ipart) + p%mD(ipart)

end subroutine drive_mass_evolution

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

   if(s%r%disc_model_type==1)then
      call thin_disk_one_zone_sg_mass(s,p,ipart,a,scale_m_msun,f_edd,eta_rad,m_sg)
   else
      write(*,*)'Unknown internal accretion disk model used...'
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

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine radiative_efficiency(a,f_edd,fedd_ADAF,f_edd_Edd,eta_rad)
   use amr_parameters, only:dp
   implicit none
   real(dp)::a,f_edd,fedd_ADAF,f_edd_Edd
   real(dp)::eta_rad
   real(dp)::A_var,B_var,C_var
   real(dp)::transition_lower,transition_higher
   ! Compute the accretion-dependent radiative efficiency
   if(f_edd>=f_edd_Edd)then
      ! Madau+2014, Sadowski+2009
      A_var = (0.9663-0.9292*a)**(-0.5639)
      B_var = (4.627 - 4.445*a)**(-0.5524)
      C_var = (827.3 - 718.1*a)**(-0.706)
      eta_rad = 0.1d0/f_edd * A_var*(0.985d0/(1.6d0/f_edd + B_var) + 0.015d0/(1.6d0/f_edd + C_var))
   else
      transition_higher = (1.0d0 + (1.88d0/f_edd)**3)**(-1)
      transition_lower = 1.0d0 - 1.0d0 * (1.0d0 + (fedd_ADAF/f_edd)**3)**(-1)
      eta_rad = radiative_efficiency_thin(a) * (1.0d0 - transition_higher - transition_lower)
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

function Lam(a)
   use amr_parameters, only:dp
   real(dp)::a,Lam
   real(dp)::Z1,Z2
   ! Compute the ISCO radius function based on black hole spin
   Z1 = 1.0d0 + (1.0d0 - a**2)**(1.0d0/3.0d0) * ((1.0d0 + a)**(1.0d0/3.0d0) + (1.0d0 - a)**(1.0d0/3.0d0))
   Z2 = sqrt(3*a**2 + Z1**2)
   if(a>=0)then
      Lam = 3.0d0 + Z2 - sqrt((3.0d0 - Z1)*(3.0d0 + Z1 + 2.0d0*Z2))
   else
      Lam = 3.0d0 + Z2 + sqrt((3.0d0 - Z1)*(3.0d0 + Z1 + 2.0d0*Z2))
   end if
end function Lam

function radiative_efficiency_thin(a)
   use amr_parameters, only:dp
   real(dp)::a,radiative_efficiency_thin
   ! Compute the radiative effiency around a spinning black hole
   radiative_efficiency_thin = 1.0d0 - sqrt(1.0d0 - 2.0d0/(3.0d0*Lam(a)))
end function radiative_efficiency_thin
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module BZ_sink_module