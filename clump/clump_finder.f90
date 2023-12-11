module clump_finder_module

contains

#if NDIM==3
subroutine m_clump_finder(pst,create_output,keep_alive)
!subroutine clump_finder(r,g,m,create_output,keep_alive)
    use amr_parameters, only:dp,ndim,twotondim
    use amr_commons, only:mesh_t,oct
    use hydro_parameters, only:nvar
    use rho_fine_module, only: m_rho_fine
    use clfind_commons
    use ramses_commons, only: pst_t,ramses_t
    use mdl_module
    use mdl_parameters
    use nbors_utils
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
    !type(run_t)::r
    type(peak_t)::p
    !type(global_t)::g
    !type(mesh_t)::m
    type(ramses_t)::s
    logical::create_output,keep_alive
    integer::ind,igrid,idim,icpu,ngrid,nleaf,nsite,now_level,next_level 
    !integer(kind=8),dimension(1:g%ncpu)::ntest_cpu
    !integer(kind=8),dimension(1:g%ncpu)::npeak_cpu
    integer::istep,nskip,ilevel,nmove,nzero,ipart,jpart,ip,ipeak
    integer::i,levelmin_part
    integer(kind=8)::ntest_tot,nmove_tot,nzero_tot
    
    logical::verbose_all=.false.
    real(kind=8)::d
    integer::action,ivar_clump
    logical::ok
    real(kind=8)::dx,vol
    integer::npeaks,npeaks_tot,icellp,icelln,ntest,ntest_all
    integer,dimension(1:ndim)::ckey,ckey_nbor
    integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
    type(oct),pointer::gridp,gridn
    real(dp)::dx_loc,scale,vol_loc
    associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,mdl=>pst%s%mdl)
#ifndef WITHOUTMPI
    integer::info
    !integer,dimension(1:g%ncpu)::nsite_cpu_tot,ntest_cpu_all,npeak_cpu_all,nsite_cpu_all
#endif
    !integer,dimension(0:g%ncpu)::nsite_cum,ntest_cum,npeak_cum
    !integer,dimension(1:g%ncpu)::nsite_cpu
    
    

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

!---------------------------------------------------------
! Compute number of test particle
!---------------------------------------------------------
! Loop over all finer levels from fine to coarse
    do i=r%nlevelmax,r%levelmin,-1
        if(m%noct_tot(i)>0)then
            ! Collect local ntest from all CPU
            call r_collect_test(pst,p,i,1)
        endif
    end do

    call r_collect_peak(s,pst,r%levelmin,p,1)

    if(npeak_cum(g%ncpu)>0)then

        ! initialize the clump type
        call r_init_peak(pst,c,p,r%levelmin,1)

        call r_collect_saddle(s,pst,p,r%levelmin,1)

        !------------------------------------------
        ! Merge irrelevant peaks
        !------------------------------------------
        if(g%myid==1.and.clinfo)write(*,*)"Now merging irrelevant peaks."

        call merge_clumps('relevance')


        !------------------------------------------
        ! Compute clumps properties
        !------------------------------------------
        if(g%myid==1.and.clinfo)write(*,*)"Computing relevant clump properties."
        !call compute_clump_properties()

        !------------------------------------------
        ! Merge clumps into haloes
        !------------------------------------------
        if(saddle_threshold>0)then
            if(g%myid==1.and.clinfo)write(*,*)"Now merging peaks into halos."
            call merge_clumps('saddleden')
        endif

        !------------------------------------------
        ! Output clumps properties to file
        !------------------------------------------
        if(r%verbose)then
            write(*,*)"Output status of peak memory."
        endif
        
        !if(clinfo.and.saddle_threshold.LE.0)call write_clump_properties(.false.)
        !if(create_output.and..not.unbind)then
        !    if(g%myid==1)write(*,*)"Outputing clump properties to disc."
            !call write_clump_properties(.true.)
            !if(r%pic)call output_part_clump_id()
            ! output the clump field
            !if (output_clump_field)then
            !    if(g%myid==1)write(*,*)"Outputing clump field to disc"
            !    !call write_clump_field
            !end if
        !endif
    endif

end associate

end subroutine m_clump_finder

#endif

subroutine neighborsearch(s,r,p,np,count,ilevel,peak)
    use mdl_module
    use amr_parameters, only:dp,ndim,nvector
    use amr_commons, only:run_t,mesh_t,oct,global_t
    use ramses_commons, only: ramses_t
    use clfind_commons
    use clump_merger, only: get_local_peak_id
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
    type(global_t)::g
    type(peak_t)::p
    type(mesh_t)::m
    type(msg_twin_realdp)::dummy_twin_realdp
    type(msg_large_realdp)::dummy_large_realdp
    !------------------------------------------------------------
    ! This routine constructs all neighboring leaf cells at levels
    ! ilevel-1, ilevel, ilevel+1.
    ! Depending on the action case value, fuctions performing
    ! further checks for the neighbor cells are called.
    ! xx is on input the array containing the density field
    !------------------------------------------------------------
    integer::np,count,ilevel,curr_level,x,y,z
    integer::i,j,k,ind,nx_loc,ipart,icellp,icelln,ipeak,jpeak
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
    real(dp),dimension(1:nPnei)::av_dens
    type(oct),pointer::gridp,gridn
    logical::ok_level,ok_leaf,ok,peak
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
        dx_loc=r%boxlen/2**p%levelp(ipart)
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
                if(peak)then
                    ! if hydro, use uold
                    if(r%hydro)then
                        if(gridn%uold(icelln,1)>density_max(i))then
                            okpeak(i)=.false.                 ! cell is no peak
                            density_max(i)=gridn%uold(icelln,1) ! change densest neighbor dens
                        endif
                    else
                        if(gridn%rho(icelln)>density_max(i))then
                            okpeak(i)=.false.                 ! cell is no peak
                            density_max(i)=gridn%rho(icelln) ! change densest neighbor dens
                        endif
                    endif
                    if(.not. okpeak(i))then
                        !change test particle properties
                        p%levelpm(ipart)=curr_level-1       ! Level
                        ! Compute peak coordinate from cell centers
                        p%maxxp(ipart,1)=(2*gridn%ckey(1)+MOD((icelln-1)  ,2)+0.5)*dx_loc-m%skip(1) 
                        p%maxxp(ipart,2)=(2*gridn%ckey(2)+MOD((icelln-1)/2,2)+0.5)*dx_loc-m%skip(2)
                        p%maxxp(ipart,3)=(2*gridn%ckey(3)+MOD((icelln-1)/4,2)+0.5)*dx_loc-m%skip(3)
                        p%denpm(ipart)=gridn%rho(icelln)
                        p%peak(ipart)=okpeak(i)
                    endif
                else
                    !neighboring cell is in a clump and is in another clump
                    if(gridn%flag2(icelln)/=0 .and. gridn%flag2(icelln)/=p%pid(ipart))then
                        if(r%hydro)then
                            av_dens(j)=(gridn%uold(icelln,1)+p%denp(ipart))/2 !average density of cell and neighbor cell
                        else
                            av_dens(j)=(gridn%uold(icelln,1)+p%denp(ipart))/2 !average density of cell and neighbor cell
                        endif
                        call get_local_peak_id(g,p%pid(ipart),ipeak)
                        call get_local_peak_id(g,gridn%flag2(icelln),jpeak)
                        if (get_value(ipeak,jpeak,sparse_saddle_dens) < av_dens(j))then
                            call set_value(ipeak,jpeak,av_dens(j),sparse_saddle_dens)
                         endif
                         if (get_value(jpeak,ipeak,sparse_saddle_dens) < av_dens(j))then
                            call set_value(jpeak,ipeak,av_dens(j),sparse_saddle_dens)
                         endif
                    endif
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


!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_collect_test(pst,p,ilevel,input_size)!ntest,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  use clfind_commons
  implicit none
  type(pst_t)::pst
  type(peak_t)::p
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  real(kind=8)::ntest,next_ntest

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RESET_RHO,pst%iUpper+1,input_size,0,ilevel)!,output_size,ilevel)
     call r_collect_test(pst%pLower,p,ilevel,input_size)!,ntest,output_size)
     call mdl_get_reply(pst%s%mdl,rID,0)!,output_size,next_ntest)
     !ntest=ntest+next_ntest
  else
     call collect_test(pst%s%r,pst%s%g,pst%s%m,p,ilevel)!,ntest)
  endif

end subroutine r_collect_test
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine collect_test(r,g,m,p,ilevel)!,ntest)
    use amr_parameters, only: twotondim
    use amr_commons, only: run_t,global_t,mesh_t
    use clfind_commons
    implicit none
    type(run_t)::r
    type(global_t)::g
    type(mesh_t)::m
    type(peak_t)::p
    integer,  intent(in)  :: ilevel
    !integer, intent(out)::ntest
    integer::ntest,ntest_all
    integer(kind=8),dimension(0:g%ncpu)::nsite_cum,ntest_cum
    integer,dimension(1:g%ncpu)::nsite_cpu,nsite_cpu_all
    integer::ind,igrid,idim,icpu,ngrid,nleaf,nsite,now_level,next_level 

    integer::istep,nskip,nmove,nzero,ipart,jpart,ip
    integer::i,levelmin_part
    integer(kind=8)::ntest_tot,nmove_tot,nzero_tot
    integer(kind=8),dimension(1:g%ncpu)::ntest_cpu,ntest_cpu_all,npeak_cpu,npeak_cpu_all
    logical::verbose_all=.false.
    real(kind=8)::d,dx_loc
    integer::action,ivar_clump
    logical::ok
    real(kind=8)::dx,vol
    ntest=0
    dx_loc=r%boxlen/2**ilevel
#ifdef GRAV
  ! Initialize density field to zero
    do igrid=m%head(ilevel),m%tail(ilevel)
        p%headp(ilevel)=ntest+1
        ! Loop over cells
        do ind=1,twotondim
            ! Select leaf cells
            ok = .not. m%grid(igrid)%refined(ind)
            ! Select dense enough cells
            !if(r%hydro)then
            !    d = m%grid(igrid)%uold(ind,ivar_clump)
            !else

            d = m%grid(igrid)%rho(ind)
            !endif
            ok = ok .and. d > density_threshold
            !count and flag !create 'testparticles'
            ! Compute test particle map
            m%grid(igrid)%flag2(ind) = 0
            if(ok)then
                ntest=ntest+1
                p%npart=ntest                  ! Local 'test particle' index
                p%levelp(p%npart)=ilevel        ! Level
                p%levelpm(p%npart)=ilevel        ! Level
                ! Compute peak coordinate from cell centers
                p%xp(p%npart,1)=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx_loc-m%skip(1) 
                p%xp(p%npart,2)=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx_loc-m%skip(2)
                p%xp(p%npart,3)=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx_loc-m%skip(3)
                p%maxxp(p%npart,1)=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx_loc-m%skip(1) 
                p%maxxp(p%npart,2)=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx_loc-m%skip(2)
                p%maxxp(p%npart,3)=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx_loc-m%skip(3)
                p%denp(p%npart)=d
                p%denpm(ipart)=d
                m%grid(igrid)%flag2(ind) = 1   
            endif
        end do
        p%tailp(ilevel)=ntest
    end do

    !---------------------------------------------------------
    ! Compute number of test particles across all CPUs.
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
        call quick_sort_dp(denp,testp_sort,ntest)
        deallocate(denp)
    endif
    
    do i=1,ntest
        p%sortp(i) = testp_sort(i)
    end do
    deallocate(testp_sort)

#endif  

end subroutine collect_test

!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_collect_peak(s,pst,p,ilevel,input_size)!ntest,output_size)
  use mdl_module
  use ramses_commons, only: pst_t,ramses_t
  use mdl_parameters
  use clfind_commons
  implicit none
  type(pst_t)::pst
  type(peak_t)::p
  type(ramses_t)::s
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RESET_RHO,pst%iUpper+1,input_size,0,ilevel)!,output_size,ilevel)
     call r_collect_peak(s,pst%pLower,p,ilevel,input_size)!,ntest,output_size)
     call mdl_get_reply(pst%s%mdl,rID,0)!,output_size,next_ntest)
  else
     call collect_peak(s,pst%s%r,pst%s%g,pst%s%m,p,ilevel)!,ntest)
  endif

end subroutine r_collect_peak
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine collect_peak(s,r,g,m,p,ilevel)
    use amr_parameters, only: twotondim
    use amr_commons, only: run_t,global_t,mesh_t,oct
    use ramses_commons, only: ramses_t
    use clfind_commons
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
    type(global_t)::g
    type(mesh_t)::m
    type(peak_t)::p
    type(msg_twin_realdp)::dummy_twin_realdp
    type(msg_large_realdp)::dummy_large_realdp
    type(oct),pointer::gridp,gridn,gridpm
    integer,  intent(in)  :: ilevel
    integer::npeaks,npeaks_tot,icpu,next_level,now_level
    integer::ipart,jpart,ip,i,icellp,icellpm,ipeak
    integer(kind=8),dimension(1:g%ncpu)::npeak_cpu,npeak_cpu_all
    integer,dimension(1:ndim)::ckey,ckey_nbor
    integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
    real(dp)::dx_loc
    
    !-----------------------------------------------------------------------
    ! Count number of density peaks and share info across processors
    !-----------------------------------------------------------------------
    !----------------------------------------------------------------------
    ! Count the peaks (cell with no denser neighbor)
    ! Store the index of the densest neighbor (can be the cell itself)
    ! for later usage.
    !----------------------------------------------------------------------
    ! Group chunks of nvector cells (of the same level) and send them
    ! to the routine that constructs the neighboring cells

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

    npeaks=0
    ip=0
    do ipart=1,p%npart
        ip=ip+1
        next_level = p%levelp(p%sortp(ipart+1))
        now_level = p%levelp(p%sortp(ipart))
        if(next_level /= now_level)then
            call neighborsearch(s,r,p,ip,npeaks,now_level,.true.)
            ip=0
        endif
    end do
    if (ip>0)then
        call neighborsearch(s,r,p,ip,npeaks,now_level,.true.)
    endif

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
    if (g%myid==1.and.npeaks_tot>0) &
        & write(*,'(" Total number of density peaks found=",I10)')npeaks_tot
    p%npeak_tot=npeaks_tot
!----------------------------------------------------------------------
! Flag peaks with global peak id using flag2 array
! Compute peak density using max_dens array
!----------------------------------------------------------------------
    ipeak = 0
    if(p%npart>0)then
        do i=1,p%npart
            ipart=p%sortp(i)
            if(p%peak(ipart))then
                ipeak=ipeak+1
                dx_loc=r%boxlen/2**p%levelp(ipart)
                ckey(1:ndim)=int(p%xp(ipart,1:ndim)/dx_loc)
                ! Get parent cell at level ilevel using cache
                hash_cell(0)=p%levelp(ipart)+1
                hash_cell(1:ndim)=ckey(1:ndim)
                call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icellp,flush_cache=.true.,fetch_cache=.true.)
                gridp%flag2(icellp)=ipeak+npeak_cum(g%myid-1)
                p%pid(ipart) = ipeak+npeak_cum(g%myid-1)
            endif
        end do
    endif

!---------------------------------------------------------------------
! Determine peak-patches around each peak
! Main step:
! - order cells in descending density
! - get peak id from densest neighbor
! - nmove is number of peak id's passed along
! - done when nmove_tot=0 (for single core, only one sweep is necessary)
!---------------------------------------------------------------------
    if (g%myid==1)write(*,*)'Finding peak patches'
    
    if(p%npart>0)then
        do ipart=1,p%npart
            jpart=p%sortp(ipart)
            if(.not. p%peak(jpart))then
                dx_loc=r%boxlen/2**p%levelp(jpart)
                ckey(1:ndim)=int(p%maxxp(jpart,1:ndim)/dx_loc)
                hash_cell(0)=p%levelpm(jpart)+1
                hash_cell(1:ndim)=ckey(1:ndim)
                call get_parent_cell(s,hash_cell,m%grid_dict,gridpm,icellpm,flush_cache=.true.,fetch_cache=.true.)
                ckey(1:ndim)=int(p%xp(jpart,1:ndim)/dx_loc)
                hash_cell(0)=p%levelp(jpart)+1
                hash_cell(1:ndim)=ckey(1:ndim)
                call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icellp,flush_cache=.true.,fetch_cache=.true.)
                gridp%flag2(icellp) = gridpm%flag2(icellp)
                p%pid(jpart) = gridpm%flag2(icellp)
            endif
        end do
    endif

    call close_cache(s,m%grid_dict)

end subroutine collect_peak

!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_collect_saddle(s,pst,p,ilevel,input_size)!ntest,output_size)
  use mdl_module
  use ramses_commons, only: pst_t,ramses_t
  use mdl_parameters
  use clfind_commons
  implicit none
  type(pst_t)::pst
  type(peak_t)::p
  type(ramses_t)::s
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RESET_RHO,pst%iUpper+1,input_size,0,ilevel)!,output_size,ilevel)
     call r_collect_saddle(s,pst%pLower,p,ilevel,input_size)!,ntest,output_size)
     call mdl_get_reply(pst%s%mdl,rID,0)!,output_size,next_ntest)
  else
     call collect_saddle(s,pst%s%r,pst%s%g,pst%s%m,p,ilevel)!,ntest)
  endif

end subroutine r_collect_saddle

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine collect_saddle(s,r,g,m,p,ilevel)
    use amr_commons
    use clfind_commons
    use amr_commons, only: run_t,global_t,mesh_t
    use ramses_commons, only: ramses_t
    implicit none
    integer::ipart,jpart,ip,now_level,next_level,ilevel
    integer::dummyint
    type(run_t)::r
    type(peak_t)::p
    type(global_t)::g
    type(mesh_t)::m
    type(ramses_t)::s
    !---------------------------------------------------------------------------
    ! subroutine which creates a npeaks**2 sized array of saddlepoint densities
    ! by looping over all testparticles and passing them to neighborcheck with
    ! case 4, which means that saddlecheck will be called for each neighboring
    ! leaf cell. There it is checked, whether the two cells (original cell and
    ! neighboring cell) are connected by a new densest saddle.
    !---------------------------------------------------------------------------
    
    ip=0
    do ipart=1,p%npart
        ip=ip+1
        jpart = p%sortp(ipart)
        now_level= p%levelp(jpart)! level
        next_level=0 !level of next particle
        if(ipart<p%npart)next_level=p%levelp(jpart+1)
        if(next_level /= now_level)then
            call neighborsearch(s,r,p,ip,dummyint,now_level,.false.)
            ip=0
        endif
    end do
    if (ip>0)then
        call neighborsearch(s,r,p,ip,dummyint,now_level,.false.)
    endif

end subroutine collect_saddle

!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_init_peak(pst,c,p,ilevel,input_size)!ntest,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  use clfind_commons
  implicit none
  type(pst_t)::pst
  type(peak_t)::p
  type(clump_t)::c
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RESET_RHO,pst%iUpper+1,input_size,0,ilevel)!,output_size,ilevel)
     call r_init_peak(pst%pLower,c,p,ilevel,input_size)!,ntest,output_size)
     call mdl_get_reply(pst%s%mdl,rID,0)!,output_size,next_ntest)
  else
     call init_peak(c,p,ilevel)!,ntest)
  endif

end subroutine r_init_peak

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_peak(c,p,ilevel)
    !use amr_commons
    use clfind_commons
    !use amr_commons, only: run_t,global_t,mesh_t
    !use ramses_commons, only: ramses_t
    implicit none
    integer::ipart,jpart,ip,now_level,next_level,ilevel
    integer::dummyint
    type(peak_t)::p
    type(clump_t)::c
    !---------------------------------------------------------------------------
    ! subroutine which initializes the clump used to merge 
    ! the tpye property idp can be used to link clump with peak_t
    !---------------------------------------------------------------------------
    ip = 0
    do i=1,p%npart
        if(p%peak(i))then
            ip=ip+1
            c%npart=ip                  ! Local 'test particle' index
            c%levelp(ip) = p%levelp(i)        ! Level
            ! Compute peak coordinate from cell centers
            c%xp(ip,1) = p%xp(i,1) 
            c%xp(ip,2) = p%xp(i,2) 
            c%xp(ip,3) = p%xp(i,3) 
            c%denp(ip) = p%denp(i)   
            c%pid = p%idp(i)
            c%idp = i
            c%idc = i
        endif
    enddo


end subroutine init_peak




end module clump_finder_module

