subroutine newdt_fine(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  use poisson_commons, ONLY: gravity_type
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-----------------------------------------------------------
  ! This routine compute the time step using 3 constraints:
  ! 1- a Courant-type condition using particle velocity
  ! 2- the gravity free-fall time
  ! 3- 10% maximum variation for aexp 
  ! This routine also compute the particle kinetic energy.
  !-----------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,nx_loc
  integer::npart1,ip,info,ilev
  integer,dimension(1:nvector),save::ind_part
  real(kind=8)::dt_loc,dt_all,ekin_loc,ekin_all,dt_acc_min
  real(dp)::tff,fourpi,threepi2
  real(dp)::aton_time_step,dt_aton,dt_rt
  real(dp)::dx_min,dx,scale,dt_fact,limiting_dt_fact
  logical::highest_level

  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  threepi2=3.0d0*ACOS(-1.0d0)**2

  ! Save old time step
  dtold(ilevel)=dtnew(ilevel)

  ! Maximum time step
  dtnew(ilevel)=boxlen/smallc
  if(poisson.and.gravity_type<=0)then
     fourpi=4.0d0*ACOS(-1.0d0)
     if(cosmo)fourpi=1.5d0*omega_m*aexp
     tff=sqrt(threepi2/8./fourpi/rho_max(ilevel))
     dtnew(ilevel)=MIN(dtnew(ilevel),courant_factor*tff)
  end if
  if(cosmo)then
     dtnew(ilevel)=MIN(dtnew(ilevel),0.1/hexp)
  end if

  if(hydro)call courant_fine(ilevel)
  
111 format('   Entering newdt_fine for level ',I2)

end subroutine newdt_fine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
#ifdef TOTO
subroutine newdt_part(ilevel, dt_loc, ekin_loc)
  use amr_commons, only: dp, icoarse_min, icoarse_max, boxlen, ndim
  use pm_commons, only: vp, mp, part_level_offset
  use hydro_commons, only: courant_factor
  implicit none
  integer      :: ilevel
  real(kind=8) :: dt_loc, ekin_loc
  ! TODO: add description here
  !
  integer  :: ipart, idim, nx_loc, offset, nparts
  real(dp) :: dx, dx_loc, scale, dtpart, v2max

  offset = part_level_offset(ilevel)
  nparts = part_level_offset(ilevel + 1) - part_level_offset(ilevel)
  
  ! Compute cell spacing
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max - icoarse_min + 1)
  scale = boxlen / dble(nx_loc)
  dx_loc = dx * scale

  ! Compute minimum time step due to particle velocities
  v2max = 0.d0
  do idim = 1, ndim
     do ipart = offset + 1, offset + nparts
        v2max = MAX(v2max, vp(ipart, idim)**2)
     end do
  end do
  
  if(v2max > 0.0D0)then
     dtpart = courant_factor * dx_loc / sqrt(v2max)
     dt_loc = MIN(dt_loc, dtpart)
  end if

  ! Compute kinetic energy (WHY IS THIS ACTUALLY DONE HERE???)
  do idim = 1, ndim
     do ipart = offset + 1, offset + nparts
        ekin_loc = ekin_loc + 0.5D0 * mp(ipart) * vp(ipart, idim)**2
     end do
  end do
    
end subroutine newdt_part
#endif



