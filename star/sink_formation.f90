module sink_formation_module

  type :: out_sink_formation_t
     real(kind=8)::mass
  end type out_sink_formation_t

contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine m_sink_formation(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_module, only: mdl_wtime
  use clump_merger_module, only: r_deallocate_clump
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  type(out_sink_formation_t)::output_sink
  double precision::ttend, ttstart

  write(*,*)'Entering sink formation'
  ttstart = mdl_wtime(pst%s%mdl)

  !----------------------------
  ! Find sink formation sites
  !----------------------------
  call m_formation_site(pst)

  !----------------------------
  ! Create sink particles
  !----------------------------
  if(pst%s%c%npeak_tot>0)then
     call r_sink_formation(pst,pst%s%r%levelmin,1,output_sink,2)
     if(output_sink%mass>0)then
        pst%s%g%mass_sink_tot=pst%s%g%mass_sink_tot+output_sink%mass
     endif
  endif

  !------------------------------
  ! Deallocate all peak arrays
  !------------------------------
  call r_deallocate_clump(pst,pst%s%r%levelmin,1)

  ttend = mdl_wtime(pst%s%mdl)
  print '(A,F14.7)',' Time elapsed in creating sinks:',ttend-ttstart

end subroutine m_sink_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_sink_formation(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(out_sink_formation_t)::output,next_output

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SINK_FORMATION,pst%iUpper+1,input_size,output_size,ilevel)
     call r_sink_formation(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
  else
     call sink_formation(pst%s%r,pst%s%g,pst%s%m,pst%s%sink,pst%s%c,output%mass)
  endif

end subroutine r_sink_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine sink_formation(r,g,m,p,c,msink_loc)
  use rng
  use constants
  use hydro_parameters, only:nvar
  use amr_parameters, only:dp,ndim,twotondim
  use amr_commons, only:run_t,global_t,mesh_t
  use pm_commons, only:part_t
  use clfind_commons, only:clump_t
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(clump_t)::c
  real(kind=8)::msink_loc
  !-------------------------------------------------------------------
  ! Spawn sink particles from clumps using various formation criteria.
  ! We use the RAMSES clump finder PHEW for the clumps detection.
  !-------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:g%ncpu)::nsite_cpu_tot,nsink_cpu_tot
#endif
  integer(kind=8),dimension(0:g%ncpu)::nsite_cum,nsink_cum
  integer,dimension(1:g%ncpu)::nsite_cpu,nsink_cpu
  integer::i,j,icpu,nsite,nsink,nsink_loc
  integer::peak_nr
  logical::ok

#if NDIM>2
  !---------------------------
  ! Count sink formation sites
  !---------------------------
  nsite=0
  c%form_sink=0
  ! Loop over peaks
  do j=1,c%npeak
     ok=.true.
     !-------------------------------------
     ! Add here all sink formation criteria
     !-------------------------------------
     if(c%relevance(j)<=c%relevance_threshold)ok=.false.
     if(c%clump_mass(j)<=c%mass_threshold*g%mp_min)ok=.false.
     if(c%occupied(j)==1)ok=.false.
!     if(r%ivar_refine>0.and.c%var_refine(j)<=r%var_cut_refine)ok=.false.
     ! Set sink formation flag
     if(ok)c%form_sink(j)=1
     if(ok)nsite=nsite+1
  end do

  !---------------------------------------------------------
  ! Compute number of sink formation sites across all CPUs.
  !---------------------------------------------------------
  nsite_cpu=0
  nsite_cpu(g%myid)=nsite
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nsite_cpu,nsite_cpu_tot,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nsite_cpu=nsite_cpu_tot
#endif
  nsite_cum=0
  do icpu=1,g%ncpu
     nsite_cum(icpu)=nsite_cum(icpu-1)+int(nsite_cpu(icpu),kind=8)
  end do

  !--------------------------
  ! Create new sink particles
  !--------------------------
  nsink_loc=0
  msink_loc=0.0d0
  ! Loop over peaks
  do j=1,c%npeak
     if(c%form_sink(j).eq.1)then
        nsink_loc=nsink_loc+1
        p%npart=p%npart+1
        if(p%npart>r%nsinkmax)then
           write(*,*)'Not enough memory for sink particle'
           write(*,*)'Increase nsinkmax in the namelist'
           stop
        endif
        ! Compute sink particle position from peak position
        p%xp(p%npart,1)=c%peak_pos(j,1)
        p%xp(p%npart,2)=c%peak_pos(j,2)
        p%xp(p%npart,3)=c%peak_pos(j,3)
        ! Compute sink particle velocity from peak velocity
        p%vp(p%npart,1)=c%peak_vel(j,1)-c%peak_acc(j,1)*0.5d0*g%dtnew(c%peak_level(j))
        p%vp(p%npart,2)=c%peak_vel(j,2)-c%peak_acc(j,2)*0.5d0*g%dtnew(c%peak_level(j))
        p%vp(p%npart,3)=c%peak_vel(j,3)-c%peak_acc(j,3)*0.5d0*g%dtnew(c%peak_level(j))
        ! Compute sink particle old force from peak acceleration
        p%fp(p%npart,1)=c%peak_acc(j,1)
        p%fp(p%npart,2)=c%peak_acc(j,2)
        p%fp(p%npart,3)=c%peak_acc(j,3)
        ! Compute sink particle mass
        p%mp(p%npart)=0
        ! Compute sink particle birth time using proper time
        p%tp(p%npart)=g%texp
        ! Compute level
        p%levelp(p%npart)=c%peak_level(j)
     endif
  end do
  p%tailp(r%nlevelmax)=p%tailp(r%nlevelmax)+nsink_loc

  !---------------------------------------------------------
  ! Compute number of new sinks across all CPUs.
  !---------------------------------------------------------
  nsink_cpu=0
  nsink_cpu(g%myid)=nsink_loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nsink_cpu,nsink_cpu_tot,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nsink_cpu=nsink_cpu_tot
#endif
  nsink_cum=0
  do icpu=1,g%ncpu
     nsink_cum(icpu)=nsink_cum(icpu-1)+int(nsink_cpu(icpu),kind=8)
  end do

  !--------------------------------------
  ! Compute new sink particle index
  !--------------------------------------
  do i=p%npart-nsink_loc+1,p%npart
     p%idp(i)=p%npart_tot+nsink_cum(g%myid-1)+i
  end do
  p%npart_tot=p%npart_tot+nsink_cum(g%ncpu)

#endif

end subroutine sink_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine m_formation_site(pst)
  use amr_parameters, only: flen
  use mdl_module, only: mdl_wtime
  use ramses_commons, only: pst_t
  use clump_merger_module, only: r_deallocate_clump
#ifdef GRAV
  use rho_fine_module, only: m_rho_fine
#endif
  implicit none
  type(pst_t)::pst
  logical::keep_alive
  !----------------------------------------------------------------------
  ! This is the master routine for the RAMSES sink formation sites finder
  !----------------------------------------------------------------------  

#if NDIM==3 && defined(GRAV)

  associate(r=>pst%s%r,g=>pst%s%g,mdl=>pst%s%mdl,p=>pst%s%p,star=>pst%s%star)

  !--------------------------------------------------------------
  ! Compute rho from gas density or dark matter or star particles
  !--------------------------------------------------------------
  call m_rho_fine(pst,r%levelmin,r%rho_type_sink)

  !----------------------------------------------
  ! Find relevant peak patches as formation sites
  !----------------------------------------------
  call r_sink_clump(pst,r%levelmin,1)

  end associate
#endif

end subroutine m_formation_site
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_sink_clump(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  
  integer::ilevel
  integer::rID
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SINK_CLUMP,pst%iUpper+1,input_size,0,ilevel)
     call r_sink_clump(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call sink_clump(pst%s)
  endif
  
end subroutine r_sink_clump
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine sink_clump(s)
  use ramses_commons, only: ramses_t
  use clump_finder_module
  use clump_merger_module
  implicit none
  type(ramses_t)::s
  
#if NDIM==3 && defined(GRAV)

  !-----------------------------------------------
  ! Store clump finder parameters in clump object.
  !-----------------------------------------------
  s%c%relevance_threshold = 3
  s%c%density_threshold = 80
  s%c%saddle_threshold = -1
  s%c%mass_threshold = 100
  !----------------------------------------------------------------------
  ! Count and collect all cells above the prescribed density threshold.
  ! We call these cell test particles for the watershed algorithm.
  !----------------------------------------------------------------------
  call collect_test(s)
  if(s%c%ntest_tot==0)return
  !----------------------------------------------------------------------
  ! Count and collect all density peaks.
  ! We also compute for each test particle the coordinates of its
  ! densest neighbor.
  !----------------------------------------------------------------------
  call collect_peak(s)
  if(s%c%npeak_tot==0)return
  !----------------------------------------------------------------------
  ! Perform a segmentation of the density field using the watershed
  ! algorithm. We get well defined peak patches around each peak.
  ! As a result, each pair of neighboring peak patches are separated
  ! by their saddle surface.
  !----------------------------------------------------------------------
  call collect_patch(s)
  !----------------------------------------------------------------------
  ! Allocate all peak patch based arrays
  !----------------------------------------------------------------------
  call allocate_peak_patch_arrays(s)
  !----------------------------------------------------------------------
  ! Update the MPI communicator for peaks
  !----------------------------------------------------------------------
  call build_peak_communicator(s)
  !----------------------------------------------------------------------
  ! We build the saddle density matrix.
  ! Each pair of peaks is connected by a unique saddle point.
  ! The saddle point is the densest point on the saddle surface,
  !----------------------------------------------------------------------
  call collect_saddle(s)
  !----------------------------------------------------------------------
  ! Update the MPI communicator for peaks
  !----------------------------------------------------------------------
  call build_peak_communicator(s)
  !----------------------------------------------------------------------
  ! Merge peaks based on a relevance criterion.
  ! Peaks that are due to random noise fluctuations or peaks that
  ! have similar peak density values are merged into relevant peaks
  !----------------------------------------------------------------------
  call merge_clumps(s,'relevance')
  !----------------------------------------------------------------------
  ! Compute relevant peak properties such as mass and number of cells
  !----------------------------------------------------------------------
  call compute_clump_properties(s,s%r%rho_type_sink)
  !----------------------------------------------------------------------
  ! Compute additional halo or particle-based clump properties.
  !----------------------------------------------------------------------
  if(s%r%rho_type_sink.eq.1)then
     call particle_clump_properties(s,s%p)
  endif
  if(s%r%rho_type_sink.eq.2)then  
     call particle_clump_properties(s,s%star)
  endif
  !---------------------------------------------
  ! Determine which peaks are occupied by a sink
  !---------------------------------------------
  call occupied_peak(s)

#endif
end subroutine sink_clump
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine occupied_peak(s)
  use ramses_commons, only: ramses_t
  use clump_merger_module
  implicit none
  type(ramses_t)::s
  integer::i,no_peak,global_peak_id,peak_nr

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,p=>s%sink)

  !--------------------------------
  ! Reads peak id of sink particles
  !--------------------------------
  call particle_peak_id(s,p,no_peak)
  !---------------------------
  ! Flag occupied peaks with 1
  !---------------------------
  c%occupied=0
  do i=1,p%npart
     global_peak_id=p%workp(i)
     if (global_peak_id /=0 ) then
        call get_local_peak_id(s,global_peak_id,peak_nr)
        c%occupied(peak_nr)=1
     end if
  end do
#ifndef WITHOUTMPI
  ! Update peak communicator
  call build_peak_communicator(s)
  ! Collect results from all MPI domains
  call virtual_peak_int(s,c%occupied,'max')
#endif

  end associate

end subroutine occupied_peak
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module sink_formation_module
