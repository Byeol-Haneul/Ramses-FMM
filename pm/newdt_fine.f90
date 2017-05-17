!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
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

  ! Save old time step
  dtold(ilevel)=dtnew(ilevel)

  ! Maximum time step
  dtnew(ilevel)=boxlen/smallc
  if(poisson.and.gravity_type<=0)then
     fourpi=4.0d0*ACOS(-1.0d0)
     if(cosmo)fourpi=1.5d0*omega_m*aexp
     threepi2=3.0d0*ACOS(-1.0d0)**2
     tff=sqrt(threepi2/8./fourpi/rho_max(ilevel))
     dtnew(ilevel)=MIN(dtnew(ilevel),courant_factor*tff)
  end if
  if(cosmo)then
     dtnew(ilevel)=MIN(dtnew(ilevel),0.1/hexp)
  end if

  ! Particle-based Courant condition
  if(pic)call newdt_part(ilevel)

  ! Hydro-based Courant condition
  if(hydro)call courant_fine(ilevel)
  
111 format('   Entering newdt_fine for level ',I2)

end subroutine newdt_fine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine newdt_part(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons, only: courant_factor
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer      :: ilevel
  !
  real(kind=8) :: dt_loc, ekin_loc, dt_all, ekin_all
  integer  :: ipart, idim, info
  real(dp) :: dx_loc, dtpart, v2max

  dt_all=dtnew(ilevel); dt_loc=dt_all
  ekin_all=0.0; ekin_loc=0.0

  ! Compute cell spacing
  dx_loc = boxlen/2**ilevel

  ! Compute minimum time step due to particle velocities
  v2max = 0.d0
  do idim = 1, ndim
     do ipart = headp(ilevel), tailp(ilevel)
        v2max = MAX(v2max, vp(ipart, idim)**2)
     end do
  end do
  
  if(v2max > 0.0D0)then
     dtpart = courant_factor * dx_loc / sqrt(v2max)
     dt_loc = MIN(dt_loc, dtpart)
  end if

  ! Compute kinetic energy
  do idim = 1, ndim
     do ipart = headp(ilevel), tailp(ilevel)
        ekin_loc = ekin_loc + 0.5D0 * mp(ipart) * vp(ipart, idim)**2
     end do
  end do
    
  ! Reduction for time step and kinetic energy
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ekin_loc,ekin_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  dt_all=dt_loc
  ekin_all=ekin_loc
#endif
  ekin_tot=ekin_tot+ekin_all
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

end subroutine newdt_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine newdt_fine_2(r,g,m,p,ilevel)
  use amr_parameters, only: dp,nvector
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
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

  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,111)ilevel

  ! Save old time step
  g%dtold(ilevel)=g%dtnew(ilevel)

  ! Maximum time step
  g%dtnew(ilevel)=r%boxlen/r%smallc
  if(r%poisson.and.r%gravity_type<=0)then
     fourpi=4.0d0*ACOS(-1.0d0)
     if(r%cosmo)fourpi=1.5d0*g%omega_m*g%aexp
     threepi2=3.0d0*ACOS(-1.0d0)**2
     tff=sqrt(threepi2/8./fourpi/g%rho_max(ilevel))
     g%dtnew(ilevel)=MIN(g%dtnew(ilevel),r%courant_factor*tff)
  end if
  if(r%cosmo)then
     g%dtnew(ilevel)=MIN(g%dtnew(ilevel),0.1/g%hexp)
  end if

  ! Particle-based Courant condition
  if(r%pic)call newdt_part_2(r,g,m,p,ilevel)

  ! Hydro-based Courant condition
  if(r%hydro)call courant_fine_2(r,g,m,ilevel)
  
111 format('   Entering newdt_fine for level ',I2)

end subroutine newdt_fine_2
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine newdt_part_2(r,g,m,p,ilevel)
  use amr_parameters, only: dp,nvector,ndim
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  integer::ilevel
  !
  real(kind=8) :: dt_loc, ekin_loc, dt_all, ekin_all
  integer  :: ipart, idim, info
  real(dp) :: dx_loc, dtpart, v2max

  dt_all=g%dtnew(ilevel); dt_loc=dt_all
  ekin_all=0.0; ekin_loc=0.0

  ! Compute cell spacing
  dx_loc = r%boxlen/2**ilevel

  ! Compute minimum time step due to particle velocities
  v2max = 0.d0
  do idim = 1, ndim
     do ipart = p%headp(ilevel), p%tailp(ilevel)
        v2max = MAX(v2max, p%vp(ipart, idim)**2)
     end do
  end do
  
  if(v2max > 0.0D0)then
     dtpart = r%courant_factor * dx_loc / sqrt(v2max)
     dt_loc = MIN(dt_loc, dtpart)
  end if

  ! Compute kinetic energy
  do idim = 1, ndim
     do ipart = p%headp(ilevel), p%tailp(ilevel)
        ekin_loc = ekin_loc + 0.5D0 * p%mp(ipart) * p%vp(ipart, idim)**2
     end do
  end do
    
  ! Reduction for time step and kinetic energy
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ekin_loc,ekin_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  dt_all=dt_loc
  ekin_all=ekin_loc
#endif
  g%ekin_tot=g%ekin_tot+ekin_all
  g%dtnew(ilevel)=MIN(g%dtnew(ilevel),dt_all)

end subroutine newdt_part_2
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################



