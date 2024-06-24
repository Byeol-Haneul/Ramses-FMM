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
  ! Call the clump finder
  !----------------------------
  call m_sink_finder(pst,.true.) ! Create no output and need to keep alive

  !----------------------------
  ! Create sink particles
  !----------------------------
  call r_sink_formation(pst,pst%s%r%levelmin,1,output_sink,2)
  if(output_sink%mass>0)then
     pst%s%g%mass_sink_tot=pst%s%g%mass_sink_tot+output_sink%mass
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
     call sink_formation(pst%s%r,pst%s%g,pst%s%m,pst%s%sink,pst%s%c,ilevel,output%mass)
  endif

end subroutine r_sink_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine sink_formation(r,g,m,p,c,ilevel,msink_loc)
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
  integer::ilevel
  real(kind=8)::msink_loc
  !-------------------------------------------------------------------
  ! Spawn star particles according to various star formation models.
  ! We use a random Poisson process.
  !-------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:g%ncpu)::nsite_cpu_tot,nsink_cpu_tot
#endif
  integer(kind=8),dimension(0:g%ncpu)::nsite_cum,nsink_cum
  integer,dimension(1:g%ncpu)::nsite_cpu,nsink_cpu
  integer::i,ind,igrid,idim,icpu,ngrid,nleaf,nsite,nsink,nsink_loc
  integer::peak_nr
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dx,vol,factG,f,d,mask
#if NENER>0
  integer::irad
#endif
  logical::ok

#if NDIM>2
  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Set some constants
  dx=r%boxlen/2**ilevel
  vol=dx**ndim
  factG=1d0
  if(r%cosmo)factG=3d0/4d0/twopi*g%omega_m*g%aexp

  !---------------------------------------------------------
  ! Count potential sink formation sites.
  !---------------------------------------------------------
  nsite=0
  ! Loop over octs with vector sweeps
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
        ! Select leaf cells
        ok = .not. m%grid(igrid)%refined(ind)
        ! Select flagged cells
        f = m%grid(igrid)%flag2(ind)
        ok = ok .and. f > 0
#ifdef HYDRO
        ! Select cells in zoom region
        if(r%ivar_refine>0)then
           mask = m%grid(igrid)%uold(ind,r%ivar_refine)/d
           ok = ok .and. mask > r%var_cut_refine
        endif
#endif
        ! Count number of formation sites
        if(ok)then
           nsite=nsite+1
        endif
     end do
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

  nsink_loc=0
  msink_loc=0.0d0
  ! Loop over octs
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
        ! Select leaf cells
        ok = .not. m%grid(igrid)%refined(ind)
        ! Select dense enough cells
        f = m%grid(igrid)%flag2(ind)
        ok = ok .and. f > 0
#ifdef HYDRO
        ! Select cells in zoom region
        if(r%ivar_refine>0)then
           mask = m%grid(igrid)%uold(ind,r%ivar_refine)/d
           ok = ok .and. mask > r%var_cut_refine
        endif
#endif
        ! Create new sink particle
        if(ok)then
           nsink_loc=nsink_loc+1
           p%npart=p%npart+1
           if(p%npart>r%nsinkmax)then
              write(*,*)'Not enough memory for sink particle'
              write(*,*)'Increase nsinkmax in the namelist'
              stop
           endif
           ! Compute sink particle coordinate from cell centers
           p%xp(p%npart,1)=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx-m%skip(1)
           p%xp(p%npart,2)=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx-m%skip(2)
           p%xp(p%npart,3)=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx-m%skip(3)
           ! Compute sink particle velocity from clump velocity
           p%vp(p%npart,1)=c%peak_vel(peak_nr,1)
           p%vp(p%npart,2)=c%peak_vel(peak_nr,1)
           p%vp(p%npart,3)=c%peak_vel(peak_nr,1)
#ifdef GRAV
           ! Remove half a kick (will be added later)
           p%vp(p%npart,1)=p%vp(p%npart,1)-m%grid(igrid)%f(ind,1)*0.5d0*g%dtnew(ilevel)
           p%vp(p%npart,2)=p%vp(p%npart,2)-m%grid(igrid)%f(ind,2)*0.5d0*g%dtnew(ilevel)
           p%vp(p%npart,3)=p%vp(p%npart,3)-m%grid(igrid)%f(ind,3)*0.5d0*g%dtnew(ilevel)
#endif
           ! Compute sink particle mass
           p%mp(p%npart)=0
           ! Compute sink particle birth time using proper time
           p%tp(p%npart)=g%texp
           ! Compute level
           p%levelp(p%npart)=ilevel

        endif
     end do
  end do
  p%tailp(r%nlevelmax)=p%tailp(r%nlevelmax)+nsink_loc

  !---------------------------------------------------------
  ! Compute number of new stars across all CPUs.
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

subroutine m_sink_finder(pst,keep_alive)
    use amr_parameters, only: flen
    use mdl_module, only: mdl_wtime
    use ramses_commons, only: pst_t
#ifdef GRAV
    use rho_fine_module, only: m_rho_fine
#endif
    implicit none
    type(pst_t)::pst
    logical::keep_alive
    !-----------------------------------------------------------------------
    ! This is the master routine for the RAMSES clump finder.
    !-----------------------------------------------------------------------  
    character(LEN=5)::nchar
    character(LEN=flen)::filename,filedir
    integer,dimension(1:flen/4)::input_array
    double precision::ttend, ttstart=0.0
    integer::dummy(1)
  
  
#if NDIM==3 && defined(GRAV)
  
    associate(r=>pst%s%r,g=>pst%s%g,mdl=>pst%s%mdl,p=>pst%s%p,star=>pst%s%star)
  
    write(*,*)'Entering clump finder'
    ttstart = mdl_wtime(mdl)
  
    !-----------------------------------------------------------------------
    ! Compute rho from gas density and/or dark matter and/or star particles
    !-----------------------------------------------------------------------
    call m_rho_fine(pst,r%levelmin,r%rtype_sink)
    
    !------------------------------------------
    ! Find relevant peak patches and halos
    !------------------------------------------
    call r_sink_finder(pst,r%levelmin,1)
  
    !------------------------------------------
    ! Deallocate all peak arrays if needed
    !------------------------------------------
    if(.not. keep_alive)then
       call r_deallocate_clump(pst,r%levelmin,1)
    endif
  
    ttend = mdl_wtime(mdl)
    print '(A,F14.7)',' Time elapsed in finding clumps:',ttend-ttstart
  
    end associate
#endif
  
  end subroutine m_sink_finder
  !################################################################
  !################################################################
  !################################################################
  !################################################################
  recursive subroutine r_sink_finder(pst,ilevel,input_size)
    use mdl_module
    use ramses_commons, only: pst_t
    use mdl_parameters
    implicit none
    type(pst_t)::pst
    integer,VALUE::input_size
  
    integer::ilevel
    integer::rID
    if(pst%nLower>0)then
       rID = mdl_send_request(pst%s%mdl,MDL_CLUMP_FINDER,pst%iUpper+1,input_size,0,ilevel)
       call r_sink_finder(pst%pLower,ilevel,input_size)
       call mdl_get_reply(pst%s%mdl,rID,0)
    else
       call sink_finder(pst%s)
    endif
    
  end subroutine r_sink_finder
  !###########################################################
  !###########################################################
  !###########################################################
  !###########################################################
  subroutine sink_finder(s)
    use ramses_commons, only: ramses_t
    use clump_finder_module
    implicit none
    type(ramses_t)::s
  
#if NDIM==3 && defined(GRAV)
  
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
    call compute_clump_properties(s)
    !----------------------------------------------------------------------
    ! Compute additional halo or particle-based clump properties.
    !----------------------------------------------------------------------
    call particle_clump_properties(s)

#endif
  end subroutine sink_finder

end module sink_formation_module
