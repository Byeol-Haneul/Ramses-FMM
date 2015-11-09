!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine kick_drift(ilevel) ! FORMERLY KNOWN AS MOVE_FINE
  use pm_commons,      only: part_level_offset, xp, vp, ap, idp
  use amr_parameters,  only: dp, ndim, tracer, hydro, static, verbose
  use amr_commons,     only: ncpu, dtnew, myid
  implicit none


  integer, intent(in) :: ilevel
  integer :: offset, nparts, idim, ipart
  
  
  if(verbose)write(*,'("Test: " (I2))')ilevel 

  offset = part_level_offset(ilevel)
  nparts = part_level_offset(ilevel + 1) - part_level_offset(ilevel)

  do ipart = offset + 1, offset + nparts
     if (idp(ipart)==1023)print*,myid,xp(ipart,1),ipart
  end do
  call compute_particle_acceleration(ilevel, tracer .and. hydro)
  
  ! Accelerate and move all parts 
  do idim = 1, ndim

     ! Update velocity     
     if(static .or. tracer)then
        do ipart = offset + 1, offset + nparts 
           vp(ipart, idim) = ap(ipart, idim)
        end do
     else
        do ipart = offset + 1, offset + nparts 
           vp(ipart, idim) = vp(ipart, idim) &
                + ap(ipart, idim) * 0.5D0 * dtnew(ilevel)
           if (idp(ipart)==1023)print*,'newmovef',idim,dtnew(ilevel),xp(ipart,idim)
        end do
     end if

     ! Update position
     if(.not. static)then
        do ipart = offset + 1, offset + nparts 
           xp(ipart, idim) = xp(ipart, idim) &
                + vp(ipart, idim) * dtnew(ilevel)
        end do
     end if

  end do
  
end subroutine kick_drift
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine second_kick(ilevel) !FORMERLY KNOWN AS SYNCHRO FINE
  use pm_commons,      only: part_level_offset, xp, vp, ap, levelp, idp
  use amr_parameters,  only: dp, nvector, ndim, tracer, hydro, static, twotondim, poisson, verbose
  use hydro_commons,   only: uold
  use poisson_commons, only: f
  use amr_commons,     only: dtnew, dtold, myid, levelmin, nlevelmax
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h' 
  integer :: info
#endif



  integer, intent(in) :: ilevel



  real(dp), dimension(1:nvector),              save :: dteff
  integer :: offset, nparts, ioft, ind, idim, ip, np

  if(verbose)write(*,'("Entering second_kick, level " I2)')ilevel 
  
  offset = part_level_offset(ilevel)
  nparts = part_level_offset(ilevel + 1) - part_level_offset(ilevel)

  call compute_particle_acceleration(ilevel, .false.)

  ! Compute individual time steps
  do ioft = offset, offset + nparts - 1, nvector
     np = min(nvector, offset + nparts - ioft)

     do ip = 1, np
        if(levelp(ioft + ip) >= ilevel)then
           dteff(ip) = dtnew(levelp(ioft + ip))
        else
           dteff(ip) = dtold(levelp(ioft + ip))
        endif
     end do

     ! Update particles level
     do ip = 1, np
        levelp(ioft + ip) = ilevel
     end do
     
     do idim = 1, ndim
        do ip = 1, np
           if (.not. abs(ap(ioft+ip, idim)) < 1000.)print*,'ap',myid, ioft + ip, ap(ioft+ip, idim)
           if (.not. abs(vp(ioft+ip, idim)) < 1000.)print*,'ap',myid, ioft + ip, ap(ioft+ip, idim)
           if (.not. abs(dteff(ip)) < 1000.)print*,'ap',myid, ioft + ip, ap(ioft+ip, idim)
           vp(ioft + ip, idim) = vp(ioft + ip, idim) &
                + ap(ioft + ip, idim) * 0.5D0 * dteff(ip)
        end do
     end do
  end do

  call MPI_BARRIER(MPI_COMM_WORLD,info) 
  if(verbose)write(*,'("Leaving second_kick, level " I2)')ilevel 
end subroutine second_kick
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine compute_particle_acceleration(ilevel, read_gas_velocity)
  use pm_commons,      only: part_level_offset, xp, ap, idp, &
                             part_hkey
  use amr_parameters,  only: dp, nvector, ndim, twotondim, poisson, verbose
  use hydro_commons,   only: uold
  use poisson_commons, only: f
  use amr_commons,     only: dtnew, ncpu, myid
  use particle_communication, only: build_communicator, part_data_to_domain_dp, domain_data_to_part_dp
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h' 
  integer :: info
#endif




  integer, intent(in) :: ilevel
  logical, intent(in) :: read_gas_velocity
  
  real(dp), allocatable, dimension(:,:) :: xp_remote, ap_remote
  integer,  dimension(1:ncpu, 1:4)       :: communicator
  integer,  dimension(1:nvector, 1:twotondim), save :: cell_index
  real(dp), dimension(1:nvector, 1:twotondim), save :: vol

  integer :: offset, nparts, ioft, np, ip, ind, idim, ipart, local_oft, npart_recv, nparts_local
  
  offset = part_level_offset(ilevel)
  nparts = part_level_offset(ilevel + 1) - part_level_offset(ilevel)

  if(verbose)write(*,'("Entering compute_particle_acceleration, level " I2)')ilevel 

  !  call compute_particle_histogram(offset, nparts)
  call check_sorted(offset,nparts)
  call build_communicator(communicator, npart_recv, &
       nparts, nparts_local, local_oft, &
       part_hkey(offset + 1 : offset + nparts, 2), &
       part_hkey(offset + 1 : offset + nparts, 1), &
       part_hkey(offset + 1 : offset + nparts, 0), & 
       ilevel)

  allocate(xp_remote(1:npart_recv, 1:3), ap_remote(1:npart_recv, 1:3))
  call part_data_to_domain_dp(communicator, xp(offset + 1 : offset + nparts, 1), xp_remote(:, 1))
  call part_data_to_domain_dp(communicator, xp(offset + 1 : offset + nparts, 2), xp_remote(:, 2))
  call part_data_to_domain_dp(communicator, xp(offset + 1 : offset + nparts, 3), xp_remote(:, 3))

  print*,'xp_remote',myid,xp_remote(1,1:3)
 ! Deal with remote particles
  do ioft = 0, npart_recv - 1, nvector
     np = min(nvector, npart_recv - ioft)

     call cic(xp_remote(ioft + 1: ioft + np, 1: ndim), cell_index, vol, np, ilevel, 2)
     
     ap_remote(ioft + 1: ioft + np, 1: ndim) = 0.0D0
     if(read_gas_velocity)then
        do idim = 1, ndim
           do ind = 1, twotondim              
              do ip = 1, np
                 ap_remote(ioft + ip, idim) = ap_remote(ioft + ip, idim) + uold(cell_index(ip,ind),idim+1) * vol(ip,ind)
              end do
           end do
        end do
     endif
     
     if(poisson)then
        do idim = 1,ndim
           do ind = 1,twotondim
              do ip = 1,np
                 ap_remote(ioft + ip, idim) = ap_remote(ioft + ip, idim) + f(cell_index(ip,ind),idim) * vol(ip,ind)
                 if (abs(xp_remote(ioft + ip,1) - 42.075531959999999) < 1.d-6 .and. idim==1)print*,'newforce',cell_index(ip,ind), vol(ip,ind),ip,ap_remote(ioft + ip,1), ioft + ip, myid
              end do
           end do
        end do
     endif
  end do
  if (myid==2)print*,'ap_remote',ap_remote(207,1:3)
  call domain_data_to_part_dp(communicator, ap_remote(:,3), ap(offset + 1 : offset + nparts, 3))
  call domain_data_to_part_dp(communicator, ap_remote(:,2), ap(offset + 1 : offset + nparts, 2))
  call domain_data_to_part_dp(communicator, ap_remote(:,1), ap(offset + 1 : offset + nparts, 1))
  if (myid==1)print*,'ap',myid,ap(260, 1:3), idp(260) 
  deallocate(xp_remote, ap_remote)


  ! Deal with local particles
  do ioft = offset + local_oft, offset + local_oft + nparts_local - 1, nvector
     np = min(nvector, offset + local_oft + nparts_local - ioft)

     call cic(xp(ioft + 1: ioft + np, 1: ndim), cell_index, vol, np, ilevel, 2)

     ap(ioft + 1: ioft + np, 1: ndim) = 0.0D0
     if(read_gas_velocity)then
        do idim = 1, ndim
           do ind = 1, twotondim              
              do ip = 1, np
                  ap(ioft + ip, idim) = ap(ioft + ip, idim) + uold(cell_index(ip,ind),idim+1) * vol(ip,ind)
              end do
           end do
        end do
     endif
     
     if(poisson)then
        do idim = 1,ndim
           do ind = 1,twotondim
              do ip = 1,np
                  ap(ioft + ip, idim) =  ap(ioft + ip, idim) + f(cell_index(ip,ind),idim) * vol(ip,ind)
                 if (idp(ioft + ip)==1023 .and. idim==1)print*,'newforce',cell_index(ip,ind), vol(ip,ind),ip
              end do
           end do
        end do
     endif
  end do

#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Just a reminder that this option is not built in yet
  call clean_stop
#endif
  call MPI_BARRIER(MPI_COMM_WORLD,info) 
  if(verbose)write(*,'("Leaving compute_particle_acceleration, level " I2)')ilevel 
end subroutine compute_particle_acceleration
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
