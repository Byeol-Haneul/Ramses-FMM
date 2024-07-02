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
  if(pst%s%c%npeak_tot>0)then
    !----------------------------
    ! Create sink particles
    !----------------------------
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
     if(c%relevance(j)<=r%sink_relevance_threshold)ok=.false.
     if(c%clump_mass(j)<=r%sink_mass_threshold*g%mp_min)ok=.false.
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
        p%vp(p%npart,1)=c%peak_vel(j,1)
        p%vp(p%npart,2)=c%peak_vel(j,2)
        p%vp(p%npart,3)=c%peak_vel(j,3)
        ! Compute sink particle old force from peak acceleration
        p%fp(p%npart,1)=c%peak_acc(j,1)
        p%fp(p%npart,2)=c%peak_acc(j,2)
        p%fp(p%npart,3)=c%peak_acc(j,3)
        ! Compute sink particle mass
        p%mp(p%npart)=0
        ! Compute sink particle birth time using proper time
        p%tp(p%npart)=g%texp
        ! Compute level
        p%levelp(p%npart)=r%nlevelmax
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
  
  !----------------------------------------------------------------------
  ! Count and collect all cells above the prescribed density threshold.
  ! We call these cell test particles for the watershed algorithm.
  !----------------------------------------------------------------------
  call collect_test(s,r%sink_density_threshold)
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
  call merge_clumps(s,'relevance',r%sink_mass_threshold,r%sink_relevance_threshold,r%sink_density_threshold,r%sink_saddle_threshold)
  !----------------------------------------------------------------------
  ! Compute relevant peak properties such as mass and number of cells
  !----------------------------------------------------------------------
  call compute_clump_properties(s,s%r%rho_type_sink)
  !----------------------------------------------------------------------
  ! Compute additional halo or particle-based clump properties.
  !----------------------------------------------------------------------
  if(s%r%rho_type_sink.eq.1)then
     call particle_clump_properties(s,s%p,r%sink_saddle_threshold,r%sink_mass_threshold,r%sink_relevance_threshold)
  endif
  if(s%r%rho_type_sink.eq.2)then  
     call particle_clump_properties(s,s%star,r%sink_saddle_threshold,r%sink_mass_threshold,r%sink_relevance_threshold)
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
subroutine move_sink(s)
    use ramses_commons, only: ramses_t
    use amr_parameters, only: dp,ndim,twotondim,twopi
    use move_fine_module, only: m_kick_drift_part
    use pm_parameters, only: action_kick_only,action_kick_drift
  implicit none
  type(ramses_t)::s
  real(dp)::dx_min
  real(dp),dimension(1:r%nsinkmax,1:ndim)    ::xsinkold,vsinkold,fsinkold
  
  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,p=>s%sink)
#if NDIM==3
#endif
  end associate
  dx_min=r%boxlen/2**r%nlevelmax/r%aexp
  xsinkold(1:p%npart,1:ndim) = p%xp(1:p%npart,1:ndim)
  vsinkold(1:p%npart,1:ndim) = p%vp(1:p%npart,1:ndim)
  fsinkold(1:p%npart,1:ndim) = p%fp(1:p%npart,1:ndim)
  call m_kick_drift_part(pst,ilevel,action_kick_only,3)
  call m_kick_drift_part(pst,ilevel,action_kick_drift,3)
  if (sink_descent) then
    xsink_graddescent(1:p%npart,1:ndim)=0.0
    fsink_norm=NORM2(p(isink,1:ndim))
    gamma_grad_descent = 0.0d0
    graddescent_over_dt = 0.0d0
    do idim=1,ndim
        gamma_grad_descent = gamma_grad_descent + (p%xp(isink,idim)-xsinkold(isink,idim))*(p%fp(isink,idim)-fsinkold(isink,idim))
    enddo
    if(gamma_grad_descent>0.0)then
        gamma_grad_descent = fudge_graddescent*g%dtnew(ilevel)*SQRT(ABS(gamma_grad_descent)/(NORM2(p%fp(isink,1:ndim)-fsinkold(isink,1:ndim)))**2)
        ! Require thatthe sink cannot move more than half a grid
        if(gamma_grad_descent*fsink_norm>dx_min/2.0) then
           xsink_graddescent(isink,1:ndim) = p%fp(isink,1:ndim) * dx_min/2.0/fsink_norm
        else
           xsink_graddescent(isink,1:ndim) = p%fp(isink,1:ndim) * gamma_grad_descent
        endif
        ! Uopdate the sink position
        p%xp(isink,1:ndim)=p%xp(isink,1:ndim)+ xsink_graddescent(isink,1:ndim)
        ! Store the descent velocity for the time-stepping
        graddescent_over_dt(isink) = NORM2(xsink_graddescent(isink,1:ndim))/g%dtnew(ilevel)
     endif
  endif


end subroutine






!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine update_sink(s,ilevel)
    use ramses_commons, only: ramses_t
    use amr_parameters, only: dp,ndim,twotondim,twopi
    use constant, only: M_sun,yr2sec
    use move_fine_module, only: m_kick_drift_part
    use pm_parameters, only: action_kick_only
  implicit none
  type(ramses_t)::s
  integer::ilevel

  real(dp),dimension(1:ndim)::r_rel
  real(dp)::rr,factG,rmax,rmax2,dx_min
  real(dp),allocatable,dimension(:)    ::msum_overlap
  logical::iyoung,jyoung,overlap,merge_flag
  integer::lev,isink,jsink
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  
  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,p=>s%sink)
#if NDIM==3

  dx_min=r%boxlen/2**r%nlevelmax/r%aexp
  rmax=dble(r%ir_cloud)*dx_min ! Linking length in physical units
  rmax2=rmax*rmax

  ! Lifetime of first larson core in code units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  t_larson1=r%merging_timescale*yr2sec/scale_t
  ! Gravitational constant
  factG=1d0
  if(r%cosmo)factG=3d0/4d0/twopi*g%omega_m*g%aexp

  ! Set overlap mass to sink mass
  allocate(msum_overlap(1:p%npart))
  msum_overlap=p%mp(p%npart)

  ! Check for overlapping sinks
  do isink=1,p%npart
    if (p%mp(isink)>0.)then
        do jsink=i+1,p%npart
            ! Compute relative distance
            r_rel(1:ndim)=p%xp(jsink,1:ndim)-p%xp(isink,1:ndim)
            do idim=1,ndim
                if (r%periodic(idim) .and. r_rel(idim)>r%boxlen/2.0d0) r_rel(idim)=r_rel(idim)-r%boxlen
                if (r%periodic(idim) .and. r_rel(idim)<-r%boxlen/2.0d0) r_rel(idim)=r_rel(idim)+r%boxlen
            end do
            rr=sum(r_rel**2)
            ! Check for overlap
            overlap=rr<4*rmax2
            if(overlap)then
                msum_overlap(isink)=msum_overlap(isink)+p%mp(jsink)
                msum_overlap(jsink)=msum_overlap(jsink)+p%mp(isink)

                ! Merging based on relative distance
                merge_flag=rr<4*dx_min**2 ! Sinks are within two cells from each other

                ! Merging based on relative velocity
                if(r%mass_merger_vel_check>0 .and. (p%mp(isink)+p%mp(jsink)).ge.r%mass_merger_vel_check*M_sun/(scale_d*scale_l**ndim)) then
                    v1_v2=(p%vp(isink,1)-p%vp(jsink,1))**2+(p%vp(isink,2)-p%vp(jsink,2))**2+(p%vp(isink,3)-p%vp(jsink,3))**2
                    merge_flag=merge_flag .and. 2*factG*(p%mp(isink)+p%mp(jsink))/sqrt(rr)>v1_v2
                end if

                ! Merging based on sink age
                if (r%merging_timescale>0d0)then
                    iyoung=(g%texp-p%tp(isink)<t_larson1)
                    jyoung=(g%texp-p%tp(jsink)<t_larson1)
                    merge_flag=merge_flag .and. (iyoung .or. jyoung)
                endif
                
                if (merge_flag.)then
                    if(g%myid==1)then
                        write(*,*)'> Merging sink ',p%idp(jsink),' into sink ',p%idp(isink)
                        if(verbose_AGN)then
                        write(*,*)'>> Sink #1: ',p%idp(isink)
                        write(*,*)p%mp(isink)/M_sun*(scale_d*scale_l**ndim)
                        write(*,*)p%xp(isink,1:ndim)
                        write(*,*)'>> Sink #2: ',p%idp(jsink)
                        write(*,*)p%mp(jsink)/M_sun*(scale_d*scale_l**ndim)
                        write(*,*)p%xp(jsink,1:ndim)
                        endif
                    endif
                    ! Set new values of remaining sink (keep one with larger index)
                    ! Compute centre of mass quantities
                    mcom     =(p%mp(isink)+p%mp(jsink))
                    xcom(1:ndim)=p%mp(isink,1:ndim)+p%mp(jsink)*r_rel(1:ndim)/mcom
                    vcom(1:ndim)=(p%mp(isink)*p%vp(isink,1:ndim)+p%mp(jsink)*p%vp(jsink,1:ndim))/mcom
                    !lcom(1:ndim)=p%mp(isink)*cross((p%xp(isink,1:ndim)-xcom(1:ndim)),p%vp(isink,1:ndim)-vcom(1:ndim))+ &
                    !      &    p%mp(jsink)*cross((p%xp(jsink,1:ndim)-xcom(1:ndim)),p%vp(jsink,1:ndim)-vcom(1:ndim))

                    ! Compute merged quantities
                    p%mp(isink)        = mcom
                    p%xp(isink,1:ndim)        = xcom(1:ndim)
                    p%vp(isink,1:ndim)        = vcom(1:ndim)
                    p%tp(isink)        = min(p%tp(isink),p%tp(jsink))
                    p%idp(isink)       = min(p%idp(isink),p%idp(jsink))

                    ! Zero mass of the sink that was merged in
                    p%mp(jsink)=0
                endif
            endif
        enddo
    endif
  enddo        

#endif
  end associate
end subroutine update_sink


!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module sink_formation_module
