module clump_finder_module

contains

subroutine clump_finder(r,g,m,create_output,keep_alive)
    use amr_parameters, only:dp,ndim,twotondim
    use amr_commons, only:run_t,global_t,mesh_t
    use clfind_commons
    use hydro_parameters, only:nvar
    use rho_fine_module, only: m_rho_fine
    use ramses_commons, only: pst_t
    use mdl_module
    use mdl_parameters
#ifndef WITHOUTMPI
    use mpi
#endif
    implicit none
    type(pst_t)::pst
    type(run_t)::r
    type(global_t)::g
    type(mesh_t)::m
    type(peak_t)::p
    !type(part_t)::s
    logical::create_output,keep_alive
#ifndef WITHOUTMPI
    integer::info
    integer,dimension(1:g%ncpu)::nsite_cpu_tot
#endif
    integer(kind=8),dimension(0:g%ncpu)::nsite_cum,ntest_cum
    integer,dimension(1:g%ncpu)::nsite_cpu,nsite_cpu_all
    integer::ind,igrid,idim,icpu,ngrid,nleaf,nsite

    integer::istep,nskip,ilevel,nmove,nzero
    integer::i,levelmin_part
    integer(kind=8)::ntest_all,nmove_tot,nzero_tot
    integer(kind=8),dimension(1:g%ncpu)::ntest_cpu,ntest_cpu_all
    integer,dimension(1:g%ncpu)::npeaks_per_cpu_tot
    logical::verbose_all=.false.
    real(kind=8)::d,d0
    integer::action,ivar_clump
    logical::ok
    real(kind=8)::dx

    ! Set some constants
    dx=r%boxlen/2**ilevel
    vol=dx**ndim

    if (create_output) then
        if(g%nstep_coarse==g%nstep_coarse_old.and.g%nstep_coarse>0)return
        if(g%nstep_coarse==0.and.r%nrestart>0)return
    endif

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


    p%npart=ntest   
    p%npart_tot=ntest_cum(g%cpu)
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
                        p%npart=p%npart+1                   ! Local 'test particle' index
                        p%levelp(p%npart)=ilevel        ! Level
                        ! Compute peak coordinate from cell centers
                        p%xp(p%npart,1)=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx-m%skip(1) 
                        p%xp(p%npart,2)=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx-m%skip(2)
                        p%xp(p%npart,3)=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx-m%skip(3)
                        p%denp(p%npart)=d ! Save density values here!
                        !m%grid(igrid)%flag2(ind)=ntest  ! Initialize flag2 to GLOBAL test particle index
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
                now_dens = p%denp(p%sortp(ipart))
                now_xp = p%xp(p%sortp(ipart))

                if(next_level /= now_level)then
                    call neighborsearch(now_dens,now_xp,ind_max,ip,npeaks,now_level,1)
                    do jpart=1,ip
                        imaxp(ind_part(jpart))=ind_max(jpart)
                    end do
                    ip=0
                endif
            end do
        endif
    endif

end subroutine clump_finder

subroutine neighborsearch(xx,ind_cell,ind_max,np,count,ilevel,action)
    use amr_commons
    implicit none
    integer::np,count,ilevel,action
    integer,dimension(1:nvector)::ind_max,ind_cell
    real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
    !------------------------------------------------------------
    ! This routine constructs all neighboring leaf cells at levels
    ! ilevel-1, ilevel, ilevel+1.
    ! Depending on the action case value, fuctions performing
    ! further checks for the neighbor cells are called.
    ! xx is on input the array containing the density field
    !------------------------------------------------------------
    integer::j,ind,nx_loc,i1,j1,k1,i2,j2,k2,i3,j3,k3,ix,iy,iz
    integer::i1min,i1max,j1min,j1max,k1min,k1max
    integer::i2min,i2max,j2min,j2max,k2min,k2max
    integer::i3min,i3max,j3min,j3max,k3min,k3max
    real(dp)::dx,dx_loc,scale,vol_loc
    integer ,dimension(1:nvector)::clump_nr,indv,ind_grid,grid,ind_cell_coarse

    real(dp),dimension(1:twotondim,1:3)::xc
    integer ,dimension(1:99)::cell_index,cell_levl,test_levl
    real(dp),dimension(1:99,1:ndim)::xtest,xrel
    logical ,dimension(1:99)::ok
    real(dp),dimension(1:nvector)::density_max
    real(dp),dimension(1:3)::skip_loc
    logical ,dimension(1:nvector)::okpeak
    integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
    integer ,dimension(1:threetondim)::nbors_father_cells_pass
    integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
    integer::ntestpos,ntp,idim,ipos


    integer(kind=8),dimension(0:ndim)::hash_key,hash_father,hash_nbor
    integer,dimension(1:3,1:8)::shift_oct=reshape(&
       & (/-1,-1,-1,+1,-1,-1,-1,+1,-1,+1,+1,-1,&
       &   -1,-1,+1,+1,-1,+1,-1,+1,+1,+1,+1,+1/),(/3,8/))
    integer,dimension(1:3,1:8)::start_oct=reshape(&
       & (/ 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1,&
       &    1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0/),(/3,8/))
#if NDIM==3
    ! Mesh spacing in that level
    dx_loc=r%boxlen/2**ilevel
    vol_loc=dx_loc**ndim
    ! Integer constants
    i1min=0; i1max=1; j1min=0; j1max=1; k1min=0; k1max=1
    i2min=0; i2max=2; j2min=0; j2max=2; k2min=0; k2max=2
    i3min=0; i3max=3; j3min=0; j3max=3; k3min=0; k3max=3
    ! Set position of cell centers relative to grid center
    do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
        if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
        if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
    end do

    call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_twin_realdp)/32,&
                pack=pack_fetch_phi,unpack=unpack_fetch_phi,&
                init=init_flush_rho, flush=pack_flush_rho, combine=unpack_flush_rho)
    ! some preliminary action...
    do j=1,np
        ckey(1:ndim)=int(p%xp(p%sortp(j),1:ndim)/dx_loc)
        !indv(j)=(ind_cell(j)-ncoarse-1)/ngridmax+1 ! cell position in grid
        ind_grid(j)=ind_cell(j)-ncoarse-(indv(j)-1)*ngridmax ! grid index
        density_max(j)=p%denp(p%sortp(j))*1.0001d0 ! get cell density (1.0001 probably not necessary)
        ind_max(j)=ind_cell(j) !save cell index
        !if (action.ge.4)clump_nr(j)=flag2(ind_cell(j)) ! save clump number
    end do

#endif
end subroutine neighborsearch
end module clump_finder_module

