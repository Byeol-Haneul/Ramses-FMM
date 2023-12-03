module clump_finder_module

contains

#if NDIM==3
subroutine m_clump_finder(pst,create_output,keep_alive)
!subroutine clump_finder(r,g,m,create_output,keep_alive)
    use amr_parameters, only:dp,ndim,twotondim
    use amr_commons, only:run_t,global_t,mesh_t
    use hydro_parameters, only:nvar
    use rho_fine_module, only: m_rho_fine
    use clfind_commons
    use ramses_commons, only: pst_t,ramses_t
    use mdl_module
    use mdl_parameters
#ifdef GRAV
  use rho_fine_module, only: m_rho_fine
#endif
  use hilbert
#ifndef WITHOUTMPI
    use mpi
    integer(i8b)::nmove_all,nzero_all
#endif
    implicit none
    type(pst_t)::pst
    type(run_t)::r
    type(global_t)::g
    type(mesh_t)::m
    type(peak_t)::p
    type(ramses_t)::s
    logical::create_output,keep_alive
#ifndef WITHOUTMPI
    integer::info
    integer,dimension(1:g%ncpu)::nsite_cpu_tot
#endif
    integer(kind=8),dimension(0:g%ncpu)::nsite_cum,ntest_cum,npeak_cum
    integer,dimension(1:g%ncpu)::nsite_cpu,nsite_cpu_all
    integer::ind,igrid,idim,icpu,ngrid,nleaf,nsite,now_level,next_level 

    integer::istep,nskip,ilevel,nmove,nzero,ipart,jpart,ip
    integer::i,levelmin_part
    integer(kind=8)::ntest_all,nmove_tot,nzero_tot
    integer(kind=8),dimension(1:g%ncpu)::ntest_cpu,ntest_cpu_all,npeak_cpu,npeak_cpu_all
    logical::verbose_all=.false.
    real(kind=8)::d,d0
    integer::action,ivar_clump
    logical::ok
    real(kind=8)::dx,vol
    integer::npeaks,npeaks_tot,npeaks_max 


    if(r%verbose)write(*,*)' Entering clump_finder'

    !-----------------------------------------------------------------------
    ! Compute rho from gas density and/or dark matter and/or star particles
    !-----------------------------------------------------------------------
    call m_rho_fine(pst,r%levelmin)


    ! Set some constants
    !dx=r%boxlen/2**ilevel
    !vol=dx**ndim

    !if (create_output) then
    !    if(g%nstep_coarse==g%nstep_coarse_old.and.g%nstep_coarse>0)return
    !    if(g%nstep_coarse==0.and.r%nrestart>0)return
    !endif

    ! if(r%verbose.and.g%myid==1)write(*,*)' Entering clump_finder'

    !---------------------------------------------------------------
    ! Compute rho from gas density or dark matter particles
    !---------------------------------------------------------------
    if(ivar_clump==0 .or. ivar_clump==-1)then
        call m_rho_fine(pst,r%levelmin)
    endif

    ntest=0
    do ilevel=r%levelmin,r%nlevelmax
        ! Loop over octs
        do igrid=m%head(ilevel),m%tail(ilevel)
            ! Loop over cells
            do ind=1,twotondim
                ! Select leaf cells
                ok = .not. m%grid(igrid)%refined(ind)
                ! Select dense enough cells
                if(r%hydro)then
                    d = m%grid(igrid)%uold(ind,ivar_clump)
                else
                    d = m%grid(igrid)%rho(ind)
                endif
                ok = ok .and. d > d0
                !count and flag !create 'testparticles'
                ! Compute test particle map
                !m%grid(igrid)%flag2(ind) = 0
                if(ok)then
                    ntest=ntest+1
                    !m%grid(igrid)%flag2(ind) = 1   
                endif
            end do
        end do
    end do
!---------------------------------------------------------
! Compute number of peaks across all CPUs.
!---------------------------------------------------------
    ntest_cpu=0
    ntest_cpu(g%myid)=ntest
#ifndef WITHOUTMPI
#ifndef LONGINT
    call MPI_ALLREDUCE(ntest_cpu,ntest_cpu_all,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
    call MPI_ALLREDUCE(ntest_cpu,ntest_cpu_all,g%ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
    ntest_cpu=ntest_cpu_all
#endif
    ntest_cum=0
    do icpu=1,g%ncpu
        ntest_cum(icpu)=ntest_cum(icpu-1)+int(ntest_cpu(icpu),kind=8)
    end do
    ntest_all = ntest_cum(g%ncpu)


    p%npart=ntest   
    p%npart_tot=ntest_cum(g%ncpu)
    itest = 0
    do ilevel=r%levelmin,r%nlevelmax
        ! Loop over octs
        do igrid=m%head(ilevel),m%tail(ilevel)
            ! Loop over cells
            do ind=1,twotondim
                ! Select leaf cells
                ok = .not. m%grid(igrid)%refined(ind)
                ! Select dense enough cells
                if(r%hydro)then
                    d = m%grid(igrid)%uold(ind,ivar_clump)
                else
                    d = m%grid(igrid)%rho(ind)
                endif
                ok = ok .and. d > d0
                if(ok)then
                    if(ntest>0)then
                        itest=itest+1
                        p%npart=p%npart+1                   ! Local 'test particle' index
                        p%levelp(p%npart)=ilevel        ! Level
                        ! Compute peak coordinate from cell centers
                        p%xp(p%npart,1)=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx-m%skip(1) 
                        p%xp(p%npart,2)=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx-m%skip(2)
                        p%xp(p%npart,3)=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx-m%skip(3)
                        p%maxxp(p%npart,1)=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx-m%skip(1) 
                        p%maxxp(p%npart,2)=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx-m%skip(2)
                        p%maxxp(p%npart,3)=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx-m%skip(3)
                        p%denp(p%npart)=d ! Save density values here!
                        !m%grid(igrid)%flag2(ind)=itest+ntest_cum(g%myid)  ! Initialize flag2 to GLOBAL test particle index
                    endif
                endif
            end do
        end do
    end do

    do i=1,ntest
        p%idp(i)=ntest_cum(g%myid-1)+i
    end do
    p%npart_tot=ntest_cum(g%ncpu)
    
    !-----------------------------------------------------------------------
    ! Sort cells above threshold according to their density
    !-----------------------------------------------------------------------
    if (ntest>0) then
        allocate(testp_sort(ntest))
        do i=1,ntest
            denp(i)=-p%denp(i)
            testp_sort(i)=i
        end do
        call quick_sort_dp(denp(1),testp_sort(1),ntest)
        deallocate(denp)
    endif
    
    do i=1,ntest
        p%sortp(i) = testp_sort(i)
    end do
    deallocate(testp_sort)

    !-----------------------------------------------------------------------
    ! Count number of density peaks and share info across processors
    !-----------------------------------------------------------------------
    npeaks=0
    if(ntest>0)then
        if(ivar_clump==0 .or. ivar_clump==-1)then  ! case 1: count peaks
            !----------------------------------------------------------------------
            ! Count the peaks (cell with no denser neighbor)
            ! Store the index of the densest neighbor (can be the cell itself)
            ! for later usage.
            !----------------------------------------------------------------------
            ! Group chunks of nvector cells (of the same level) and send them
            ! to the routine that constructs the neighboring cells
            ip=0
            do ipart=1,ntest
                ip=ip+1
                next_level = p%levelp(p%sortp(ipart+1))
                now_level = p%levelp(p%sortp(ipart))

                if(next_level /= now_level)then
                    call neighborsearch(r,p,ip,npeaks,now_level)
                    ip=0
                endif
            end do
            if (ip>0)then
                call neighborsearch(r,p,ip,npeaks,now_level)
            endif
        endif
    endif

!---------------------------------------------------------
! Compute number of peaks across all CPUs.
!---------------------------------------------------------
    npeak_cpu=0
    npeak_cpu(g%myid)=npeaks
#ifndef WITHOUTMPI
#ifndef LONGINT
    call MPI_ALLREDUCE(npeak_cpu,npeak_cpu_all,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
    call MPI_ALLREDUCE(npeak_cpu,npeak_cpu_all,g%ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
    npeak_cpu=npeak_cpu_all
#endif
    npeak_cum=0
    do icpu=1,g%ncpu
        npeak_cum(icpu)=npeak_cum(icpu-1)+int(npeak_cpu(icpu),kind=8)
    end do
#ifdef WITHOUTMPI
    npeaks_tot=npeaks
#else
    npeaks_tot=npeak_cum(g%ncpu)
#endif
    if (myid==1.and.npeaks_tot>0) &
        & write(*,'(" Total number of density peaks found=",I10)')npeaks_tot


!----------------------------------------------------------------------
! Flag peaks with global peak id using flag2 array
! Compute peak density using max_dens array
!----------------------------------------------------------------------
ipeak = 0
if(ntest>0)then
    if(ivar_clump==0 .or. ivar_clump==-1)then
        do i=1,ntest
            ipart=p%sortp(i)
            if(p%peak(ipart))then
               ipeak=ipeak+1
               ckey(1:ndim)=int(p%xp(ipart,1:ndim)/dx_loc)
               ! Get parent cell at level ilevel using cache
                hash_cell(0)=ilevel+1
                hash_cell(1:ndim)=ckey(1:ndim)
                call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icellp,flush_cache=.true.,fetch_cache=.true.)
                gridp%flag2(icellp)=ipeak+npeak_cum(g%myid)
            endif
        end do
    else
        if(hydro)then
            call flag_peaks(uold(1,ivar_clump),nskip)
        endif
    endif
endif



!---------------------------------------------------------------------
! Determine peak-patches around each peak
! Main step:
! - order cells in descending density
! - get peak id from densest neighbor
! - nmove is number of peak id's passed along
! - done when nmove_tot=0 (for single core, only one sweep is necessary)
!---------------------------------------------------------------------
if (myid==1.and.ntest_all>0)write(*,*)'Finding peak patches'






end subroutine m_clump_finder

#endif

subroutine neighborsearch(r,p,np,count,ilevel)
    use mdl_module
    use amr_parameters, only:dp,ndim,nvector
    use amr_commons, only:run_t,mesh_t,oct
    use clfind_commons
    use ramses_commons, only: pst_t,ramses_t
    use multigrid_fine_coarse, only:pack_fetch_phi,unpack_fetch_phi
    use rho_fine_module, only:init_flush_rho,pack_flush_rho,unpack_flush_rho
    use godunov_fine_module, only: init_flush_godunov,pack_flush_godunov,unpack_flush_godunov
    use marshal, only: pack_fetch_refine,unpack_fetch_refine
    use boundaries, only: init_bound_refine
    use cache_commons
    use cache
    use nbors_utils
    implicit none
    type(ramses_t)::s
    type(run_t)::r
    type(peak_t)::p
    type(mesh_t)::m
    type(msg_twin_realdp)::dummy_twin_realdp
    !------------------------------------------------------------
    ! This routine constructs all neighboring leaf cells at levels
    ! ilevel-1, ilevel, ilevel+1.
    ! Depending on the action case value, fuctions performing
    ! further checks for the neighbor cells are called.
    ! xx is on input the array containing the density field
    !------------------------------------------------------------
    integer::np,count,ilevel,curr_level,x,y,z
    integer::i,j,k,ind,nx_loc,ipart,icellp,icelln
    real(dp)::dx,dx_loc,scale,vol_loc
    integer ,dimension(1:nvector)::ind_grid,grid
    integer ,dimension(1:99)::cell_index,cell_levl,test_levl
    real(dp),dimension(1:99,1:ndim)::xtest,xrel
    real(dp),dimension(1:nvector)::density_max
    real(dp),dimension(1:3)::skip_loc
    logical ,dimension(1:nvector)::okpeak
    integer::ntestpos,ntp,idim,ipos
    integer,dimension(1:ndim)::ckey,ckey_nbor
    real(dp),dimension(1:ndim)::xcen,xnei
    ! Number of neighboring cells to deposit mass/momentum/energy
    integer, parameter::nPnei=48
    real(dp),dimension(1:3,1:nPnei)::xPnei
    type(oct),pointer::gridp,gridn
    logical::ok_level,ok_leaf,ok
    integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor

    
    ! Arrays to define neighbors (center=[0,0,0])
    ! normalized to dx = 1 = size of the central leaf cell in which a SN particle sits
    ! from -0.75 to 0.75
    ind=0
    do k=1,4
        do j=1,4
            do i=1,4
                ok=.true.
                if((i==1.or.i==4).and.(j==1.or.j==4).and.(k==1.or.k==4)) ok=.false. ! edge
                if((i==2.or.i==3).and.(j==2.or.j==3).and.(k==2.or.k==3)) ok=.false. ! centre
                if(ok)then
                    ind = ind+1
                    x = (i-1)+0.5d0 - 2
                    y = (j-1)+0.5d0 - 2
                    z = (k-1)+0.5d0 - 2
                    xPnei(1,ind) = x/2d0
                    xPnei(2,ind) = y/2d0
                    xPnei(3,ind) = z/2d0
                endif
            enddo
        enddo
    enddo

#if NDIM==3
    ! Mesh spacing in that level
    dx_loc=r%boxlen/2**ilevel
    vol_loc=dx_loc**ndim

    if(r%hydro)then
        call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_large_realdp)/32,&
                pack=pack_fetch_refine,unpack=unpack_fetch_refine,&
                init=init_flush_godunov, flush=pack_flush_godunov,&
                combine=unpack_flush_godunov, bound=init_bound_refine)
    else
        call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_twin_realdp)/32,&
                pack=pack_fetch_phi,unpack=unpack_fetch_phi,&
                init=init_flush_rho, flush=pack_flush_rho, combine=unpack_flush_rho)
    endif
    okpeak = .true.
    do i=1,np
        ipart = p%sortp(i)
        density_max(i)=p%denp(ipart) ! get cell density (1.0001 probably not necessary)
        p%peak(ipart)=okpeak(i)
        ! Set pointers to null
        icellp=0; icelln=0
        !nullify(gridp)
        !nullify(gridn)
        !do j=1,nPnei
        !    nullify(grid_nbor(j)%p)
        !end do

        ! Find parent cell at level ilevel
        ckey(1:ndim)=int(p%xp(ipart,1:ndim)/dx_loc)
        xcen(1:ndim)=ckey(1:ndim)+0.5
        ! Get parent cell at level ilevel using cache
        hash_cell(0)=ilevel+1
        hash_cell(1:ndim)=ckey(1:ndim)
        call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icellp,flush_cache=.true.,fetch_cache=.true.,lock=.true.)

        !ind_max(j)=ind_cell(j) !save cell index
        ok_level = associated(gridp)
        if(.not. ok_level)then
            write(*,*)"Something went wrong in neighbor searching in peak finding"
            write(*,*)"Current level grid should exist..."
            stop
        endif

        ! Collect all neighboring cell from hash table
        do j=1,nPnei
            curr_level = ilevel+1
            ! Compute neighboring cell coordinates
            xnei(1:ndim)=xcen(1:ndim)+xPnei(1:ndim,j)
            ! Periodic boundary conditions
            do idim=1,ndim
                if(xnei(idim)<0.0d0)xnei(idim)=xnei(idim)+m%ckey_max(ilevel+1)
                if(xnei(idim)>=m%ckey_max(ilevel+1))xnei(idim)=xnei(idim)-m%ckey_max(ilevel+1)
            end do
            ! Get neighboring cell at ilevel
            ckey_nbor(1:ndim)=int(xnei(1:ndim))
            hash_nbor(0)=ilevel+1
            hash_nbor(1:ndim)=ckey_nbor(1:ndim)
            call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.,lock=.true.)
            ! If missing, get neighboring cell at ilevel-1
            if(.not.associated(gridn))then
                call unlock_cache(s,gridn)
                ckey_nbor(1:ndim)=int(xnei(1:ndim)/2.0)
                hash_nbor(0)=ilevel
                hash_nbor(1:ndim)=ckey_nbor(1:ndim)
                curr_level = ilevel
                call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.,lock=.true.)
            ! If refined, get neighboring cell at ilevel+1
            else if (gridn%refined(icelln))then
                call unlock_cache(s,gridn)
                ckey_nbor(1:ndim)=int(xnei(1:ndim)*2.0)
                hash_nbor(0)=ilevel+2
                hash_nbor(1:ndim)=ckey_nbor(1:ndim)
                curr_level = ilevel+2
                call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.,lock=.true.)
            endif
            !TO DO:
            ! if redfined, needs to loop for redfinding results
            
            if (.not. gridn%refined(icelln))then
                ! if hydro, use uold
                if(r%hydro)then
                    if(gridn%uold(icelln,ivar_clump)>density_max(i))then
                        okpeak(i)=.false.                 ! cell is no peak
                        density_max(i)=gridn%uold(icelln,ivar_clump) ! change densest neighbor dens
                    endif
                else
                    if(gridn%rho(icelln)>density_max(i))then
                        okpeak(i)=.false.                 ! cell is no peak
                        density_max(i)=gridn%rho(icelln) ! change densest neighbor dens
                    endif
                endif
                if(.not. okpeak(i))then
                    !change test particle properties
                    p%levelp(ipart)=curr_level-1       ! Level
                    ! Compute peak coordinate from cell centers
                    p%maxxp(ipart,1)=(2*gridn%ckey(1)+MOD((icelln-1)  ,2)+0.5)*dx-m%skip(1) 
                    p%maxxp(ipart,2)=(2*gridn%ckey(2)+MOD((icelln-1)/2,2)+0.5)*dx-m%skip(2)
                    p%maxxp(ipart,3)=(2*gridn%ckey(3)+MOD((icelln-1)/4,2)+0.5)*dx-m%skip(3)
                    p%denp(ipart)=gridn%rho(icelln)
                    p%peak(ipart)=okpeak(i)
                endif
            endif
        enddo
        ! Unlock all octs
        call unlock_cache(s,gridp)
    end do
    call close_cache(s,m%grid_dict)
    do j=1,np
        if(okpeak(j))then
            count=count+1
        endif
    end do
#endif
end subroutine neighborsearch
end module clump_finder_module

