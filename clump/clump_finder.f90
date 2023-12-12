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
    type(clump_t)::c
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

    call r_collect_peak(s,pst,p,r%levelmin,1)

    if(npeak_cum(g%ncpu)>0)then

        ! initialize the clump type
        call r_init_peak(pst,c,p,r%levelmin,1)

        call r_collect_saddle(s,pst,p,r%levelmin,1)

        !------------------------------------------
        ! Merge irrelevant peaks
        !------------------------------------------
        if(g%myid==1.and.clinfo)write(*,*)"Now merging irrelevant peaks."

        call r_merge_clumps(s,pst,c,p,'relevance',r%levelmin,1)


        !------------------------------------------
        ! Compute clumps properties
        !------------------------------------------
        if(g%myid==1.and.clinfo)write(*,*)"Computing relevant clump properties."
        call r_compute_clump_properties(s,pst,c,p,r%levelmin,1)

        !------------------------------------------
        ! Merge clumps into haloes
        !------------------------------------------
        if(saddle_threshold>0)then
            if(g%myid==1.and.clinfo)write(*,*)"Now merging peaks into halos."
            call r_merge_clumps(s,pst,c,p,'saddleden',r%levelmin,1)
        endif

        !------------------------------------------
        ! Output clumps properties to file
        !------------------------------------------
        if(r%verbose)then
            write(*,*)"Output status of peak memory."
        endif
        
        !if(clinfo.and.saddle_threshold.LE.0)call write_clump_properties(r,g,mdl,c,.false.)
        !if(create_output.and..not.unbind)then
        !    ! if unbind, output will be written in unbinding() routine
        !    if(g%myid==1)write(*,*)"Outputing clump properties to disc."
        !    call write_clump_properties(r,g,mdl,c,.true.)
            !if(r%pic)call output_part_clump_id()
            ! output the clump field
            !if (output_clump_field)then
            !    if(g%myid==1)write(*,*)"Outputing clump field to disc"
            !    call write_clump_field
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
    !use clump_merger, only: get_local_peak_id
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
                if(r%hydro)then
                    p%vel(p%npart,1)=m%grid(igrid)%uold(ind,2)
                    p%vel(p%npart,2)=m%grid(igrid)%uold(ind,3)
                    p%vel(p%npart,3)=m%grid(igrid)%uold(ind,4)
                endif
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
    integer::ipart,jpart,ip,now_level,next_level,ilevel,i
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
            c%idc = p%idp(i)
        endif
    enddo


end subroutine init_peak

!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_local_peak_id(g,global_peak_id,local_peak_id)
    use amr_commons
    use clfind_commons
    implicit none
    type(global_t)::g
    integer::local_peak_id
    integer::ihash,ikey,jkey
    integer::global_peak_id

    if(    global_peak_id> npeak_cum(g%myid-1) .and. &
         & global_peak_id<=npeak_cum(g%myid))then
       local_peak_id=global_peak_id-npeak_cum(g%myid-1)
    else
       ihash=MOD(global_peak_id,nhash)+1 ! compute the simple prime hash key
       if(hkey(ihash)==0)then ! hash table is empty
          hkey(ihash)=hfree
          local_peak_id=hfree
          gkey(hfree)=global_peak_id
          hfree=hfree+1
          if(hfree.eq.npeaks_max)then
             write(*,*)'Too many peaks'
             write(*,*)'Increase npeaks_max'
             stop
          endif
       else
          ikey=hkey(ihash) ! collision in the hash table
          do while(ikey>0)
             jkey=ikey
             if(gkey(ikey)==global_peak_id)exit
             ikey=nkey(ikey)
          end do
          if(ikey==0)then ! peak doesn't already exist
             nkey(jkey)=hfree
             local_peak_id=hfree
             gkey(hfree)=global_peak_id
             hfree=hfree+1
             hcollision=hcollision+1
             if(hfree.eq.npeaks_max)then
                write(*,*)'Too many peaks'
                write(*,*)'Increase npeaks_max'
                stop
             endif
          else            ! peak already exists
             local_peak_id=ikey
          end if
       end if
    end if
  
end subroutine get_local_peak_id




recursive subroutine r_compute_clump_properties(s,pst,c,p,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t,ramses_t
  use mdl_parameters
  use clfind_commons
  implicit none
  type(pst_t)::pst
  type(peak_t)::p
  type(clump_t)::c
  type(ramses_t)::s
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel


  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RESET_RHO,pst%iUpper+1,input_size,0,ilevel)!,output_size,ilevel)
     call r_compute_clump_properties(s,pst%pLower,c,p,ilevel,input_size)!,ntest,output_size)
     call mdl_get_reply(pst%s%mdl,rID,0)!,output_size,next_ntest)
  else
     call compute_clump_properties(pst%s%r,pst%s%g,c,p,ilevel)!,ntest)
  endif

end subroutine r_compute_clump_properties
!################################################################
!################################################################
!################################################################
!################################################################



subroutine compute_clump_properties(r,g,c,p,ilevel)
    use amr_commons
    use clfind_commons
    !use multigrid_fine_coarse, only:pack_fetch_phi,unpack_fetch_phi
    !use rho_fine_module, only:init_flush_rho,pack_flush_rho,unpack_flush_rho
    !use godunov_fine_module, only: init_flush_godunov,pack_flush_godunov,unpack_flush_godunov
    !use marshal, only: pack_fetch_refine,unpack_fetch_refine
    !use boundaries, only: init_bound_refine
    !use cache_commons
    !use cache
    !use nbors_utils
    implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
    type(run_t)::r
    type(global_t)::g
    type(peak_t)::p
    type(clump_t)::c
    !type(msg_twin_realdp)::dummy_twin_realdp
    !type(msg_large_realdp)::dummy_large_realdp
    !type(oct),pointer::gridp,gridn,gridpm
    !integer,dimension(1:ndim)::ckey,ckey_nbor
    !integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
    real(dp)::dx_loc,vol_loc
    integer::ipart,grid,peak_nr,ilevel,global_peak_id,ipeak,plevel
    real(dp)::zero=0
    !variables needed temporarily store cell properties
    real(dp)::d=0, vol=0
    ! variables related to the size of a cell on a given level
    real(dp)::dx
    real(dp),dimension(1:r%nlevelmax)::volume
    real(dp),dimension(1:3)::skip_loc,xcell
    real(dp),dimension(1:twotondim,1:3)::xc
    integer::nx_loc,ind,ix,iy,iz,idim,npeaks,ntest
    logical,dimension(1:ndim)::period
    logical::periodic
    

#ifndef WITHOUTMPI
    integer::i
    real(dp)::tot_mass_tot
#endif

    periodic=r%periodic(1)
    periodic=periodic.or.r%periodic(2)
    periodic=periodic.or.r%periodic(3)

    !peak-patch related arrays before sharing information with other cpus

    npeaks = c%npart
    ntest = p%npart

    !min_dens=huge(zero)
    !n_cells=0; n_cells_halo=0
    !halo_mass=0d0; clump_mass=0d0; clump_vol=0d0
    !center_of_mass=0d0; clump_velocity=0d0
    !peak_pos=0d0

    if(r%verbose)write(*,*)'Entering compute clump properties'


    !if(r%hydro)then
    !    call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
    !            hilbert=m%domain,pack_size=storage_size(dummy_large_realdp)/32,&
    !            pack=pack_fetch_refine,unpack=unpack_fetch_refine,&
    !            init=init_flush_godunov, flush=pack_flush_godunov,&
    !            combine=unpack_flush_godunov, bound=init_bound_refine)
    !else
    !    call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
    !            hilbert=m%domain,pack_size=storage_size(dummy_twin_realdp)/32,&
    !            pack=pack_fetch_phi,unpack=unpack_fetch_phi,&
    !            init=init_flush_rho, flush=pack_flush_rho, combine=unpack_flush_rho)
    !endif

    !--------------------------------------------------------------------------
    ! loop over all cells above the threshold
    !--------------------------------------------------------------------------
    do ipart=1,ntest
        global_peak_id=p%pid(ipart)
        if (global_peak_id /=0 ) then
            call get_local_peak_id(g,global_peak_id,peak_nr)
                
            dx_loc=r%boxlen/2**p%levelp(ipart)
            vol_loc=dx_loc**ndim
            !ckey(1:ndim)=int(c%xp(peak_nr,1:ndim)/dx_loc)
            !hash_cell(0)=c%levelp(peak_nr)+1
            !hash_cell(1:ndim)=ckey(1:ndim)
            !call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icellp,flush_cache=.true.,fetch_cache=.true.)
             
            d = p%denp(ipart)!c%denp(peak_nr)
    
            ! Cell volume
            vol=vol_loc
    
            ! Number of leaf cells per clump
            c%n_cells(peak_nr)=c%n_cells(peak_nr)+1
    
            ! Min density
            c%min_dens(peak_nr)=min(d,c%min_dens(peak_nr))
    
            ! Clump mass
            c%clump_mass(peak_nr)=c%clump_mass(peak_nr)+vol*d
    
            ! Clump volume
            c%clump_vol(peak_nr)=c%clump_vol(peak_nr)+vol
    
            ! Clump center of mass location
            c%center_of_mass(peak_nr,1:3)=c%center_of_mass(peak_nr,1:3)+vol*d*c%xp(peak_nr,1:3)
    
            ! Clump center of mass velocity

            if (r%hydro)then
                c%clump_velocity(peak_nr,1:3)=c%clump_velocity(peak_nr,1:3)+vol*p%vel(ipart,1:3)
            endif
        end if
    end do
    
    !--------------------------------------------------------------------------
    ! Loop over local peaks and identify true peak positions
    !--------------------------------------------------------------------------
    do ipeak=1,npeaks
        ! Peak cell coordinates
        plevel=c%levelp(ipeak)
        dx=0.5D0**plevel
        xcell(1:ndim)=c%xp(ipeak,1:3)
        !call true_max(xcell(1),xcell(2),xcell(3),plevel)
        c%peak_pos(ipeak,1:3)=xcell(1:3)
    end do

    ! Compute specific quantities
    do ipeak=1,npeaks
        if (c%relevance(ipeak)>0..and.c%n_cells(ipeak)>0)then
        c%center_of_mass(ipeak,1:3)=c%center_of_mass(ipeak,1:3)/c%clump_mass(ipeak)
        c%clump_velocity(ipeak,1:3)=c%clump_velocity(ipeak,1:3)/c%clump_mass(ipeak)
        end if
        ! Initialize halo mass to clump mass
        c%halo_mass(ipeak)=c%clump_mass(ipeak)
    end do

    ! Calculate total mass above threshold
    tot_mass=sum(c%clump_mass(1:npeaks))

    ! Compute further properties of the clumps
    !do ipeak=1,npeaks
    !    if (c%relevance(ipeak)>0..and.c%n_cells(ipeak)>0)then
    !    av_dens(ipeak)=c%clump_mass(ipeak)/c%clump_vol(ipeak)
    !    end if
    !end do

    ! For periodic boxes, recompute center of mass relative to peak position
    if(periodic)then
        !center_of_mass=0d0;
        do ipart=1,ntest
            global_peak_id=p%pid(ipart)!flag2(icellp(ipart))
            if (global_peak_id /=0 ) then
                call get_local_peak_id(g,global_peak_id,peak_nr)

                dx_loc=r%boxlen/2**p%levelp(ipart)
                vol_loc=dx_loc**ndim
                xcell(1:ndim)=p%xp(ipart,1:ndim)

                do idim=1,ndim
                    if (period(idim) .and. (xcell(idim)-c%peak_pos(peak_nr,idim))>r%boxlen*0.5)xcell(idim)=xcell(idim)-r%boxlen
                    if (period(idim) .and. (xcell(idim)-c%peak_pos(peak_nr,idim))<r%boxlen*(-0.5))xcell(idim)=xcell(idim)+r%boxlen
                end do

                d = p%denp(ipart)

                ! Cell volume
                vol=vol_loc

                ! Clump center of mass location
                c%center_of_mass(peak_nr,1:3)=c%center_of_mass(peak_nr,1:3)+vol*d*xcell(1:3)

            end if
        end do

        ! Compute specific quantity
        do ipeak=1,npeaks
            if (c%relevance(ipeak)>0..and.c%n_cells(ipeak)>0)then
            c%center_of_mass(ipeak,1:3)=c%center_of_mass(ipeak,1:3)/c%clump_mass(ipeak)
            end if
        end do
    endif



end subroutine compute_clump_properties
!################################################################
!################################################################
!################################################################
!################################################################



recursive subroutine r_merge_clumps(s,pst,c,p,action,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t,ramses_t
  use mdl_parameters
  use clfind_commons
  implicit none
  character(len=9)::action
  type(pst_t)::pst
  type(peak_t)::p
  type(clump_t)::c
  type(ramses_t)::s
  type(mesh_t)::m
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RESET_RHO,pst%iUpper+1,input_size,0,ilevel)!,output_size,ilevel)
     call r_merge_clumps(s,pst%pLower,c,p,action,ilevel,input_size)!,ntest,output_size)
     call mdl_get_reply(pst%s%mdl,rID,0)!,output_size,next_ntest)
  else
     call merge_clumps(s,pst%s%r,pst%s%g,pst%s%m,c,p,action,ilevel)!,ntest)
  endif

end subroutine r_merge_clumps


subroutine merge_clumps(s,r,g,m,c,p,action,ilevel)
    use amr_commons
    use clfind_commons
    use ramses_commons, only:ramses_t
    use multigrid_fine_coarse, only:pack_fetch_phi,unpack_fetch_phi
    use rho_fine_module, only:init_flush_rho,pack_flush_rho,unpack_flush_rho
    use godunov_fine_module, only: init_flush_godunov,pack_flush_godunov,unpack_flush_godunov
    use marshal, only: pack_fetch_refine,unpack_fetch_refine
    use boundaries, only: init_bound_refine
    use cache_commons
    use cache
    use nbors_utils
    implicit none
    character(len=9)::action
    type(ramses_t)::s
    type(run_t)::r
    type(global_t)::g
    type(mesh_t)::m
    type(peak_t)::p
    type(clump_t)::c
    type(msg_twin_realdp)::dummy_twin_realdp
    type(msg_large_realdp)::dummy_large_realdp
    type(oct),pointer::gridp,gridn,gridpm
    integer,dimension(1:ndim)::ckey,ckey_nbor
    integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
    integer::ilevel
    real(dp)::dx_loc
#ifndef WITHOUTMPI
    integer::info
#endif
    !---------------------------------------------------------------------------
    ! This routine merges the irrelevant clumps
    ! -clumps are sorted by ascending max density
    ! -irrelevent clumps are merged to most relevant neighbor
    !---------------------------------------------------------------------------
    integer::j,i,merge_to,ipart,icellp,mergelevel_max,npeaks,npeaks_tot
    integer::current,nmove,ipeak,jpeak,iter
    integer::nsurvive,nzero,idepth
    integer::ilev,global_peak_id,n_cells_halo
    real(dp)::value_iij,zero=0,relevance_peak
    integer,dimension(1:npeaks_max)::alive,ind_sort
    real(dp),dimension(1:npeaks_max)::peakd
    logical::do_merge=.false.
    real(dp)::halo_mass

#ifndef WITHOUTMPI
    integer::mergelevel_max_global
    integer::nmove_all,nsurvive_all,nzero_all
#endif
    if(r%verbose)then
        if(action.EQ.'relevance')then
        write(*,*)'Now merging irrelevant clumps'
        endif
        if(action.EQ.'saddleden')then
        write(*,*)'Now merging clumps into halos'
        endif
    endif
    !npeaks = npeak_cpu(g%myid)-npeak_cpu(g%myid-1)
    ! Initialize new_peak array to global peak id
    ! All peaks are alive at the start
    npeaks = c%npart 
    do i=1,npeaks
        if(action.EQ.'relevance')then
            c%alive(i) = 1
        endif
        if(action.EQ.'saddleden')then
            if(c%relevance(i)>relevance_threshold)then
                c%alive(i) = 1
            else
                c%alive(i) = 0
            endif
        endif
    end do
    ! Sort peaks by maximum peak density in ascending order
    
    if (npeaks>0) then
        allocate(testp_sort(npeaks))
        do i=1,npeaks
            denp(i)=c%denp(i)
            testp_sort(i)=i
        end do
        call quick_sort_dp(denp,testp_sort,npeaks)
        deallocate(denp)
    endif

    do i=1,npeaks
        c%sortp(i) = testp_sort(i)
    end do
    deallocate(testp_sort)

    ! Loop over peak levels
    nzero=npeaks_tot
    idepth=0
    do while(nzero>0)

        ! Compute maximum saddle density for each clump
        do i=1,hfree-1
            call get_max(g,i,sparse_saddle_dens)
        end do
        
        ! Merge peaks
        nmove=npeak_cum(g%ncpu)
        iter=0
        do while(nmove>0)
            nmove=0
            do i=npeaks,1,-1
                ipeak=c%sortp(i)
                merge_to=c%idc(ipeak)
                if(c%alive(ipeak)>0)then
                    if(action.EQ.'relevance')then
                        if(sparse_saddle_dens%maxval(ipeak)>0)then
                            relevance_peak=c%denp(ipeak)/sparse_saddle_dens%maxval(ipeak)
                        else
                            relevance_peak=c%denp(ipeak)/density_threshold
                         end if
                        do_merge=(relevance_peak<relevance_threshold)
                    endif
                    if(action.EQ.'saddleden')then
                        do_merge=(sparse_saddle_dens%maxval(ipeak)>saddle_threshold)
                    endif
                    if(do_merge)then
                        if(sparse_saddle_dens%maxloc(ipeak)>0)then
                           call get_local_peak_id(g,sparse_saddle_dens%maxloc(ipeak),jpeak)
                           if(c%denp(jpeak)>c%denp(ipeak))then
                                merge_to=c%idc(jpeak)
                           else if(c%denp(jpeak)==c%denp(ipeak))then
                                merge_to=MIN(c%idc(ipeak),c%idc(jpeak))
                           endif
                        endif
                    endif
                endif
                if(c%idc(ipeak).NE.merge_to)then
                    nmove=nmove+1
                    c%idc(ipeak)=merge_to
                endif
            enddo
            iter=iter+1
        enddo
        ! Transfer matrix elements of merged peaks to surviving peaks
        ! Create new local duplicated peaks and update communicator
        do ipeak=1,hfree-1
            if(c%alive(ipeak)>0)then
                merge_to=c%idc(ipeak)
                if(ipeak.LE.npeaks)then
                    global_peak_id=npeak_cum(g%myid-1)+ipeak
                else
                    global_peak_id=gkey(ipeak)
                endif
                if(merge_to.NE.global_peak_id)then
                    call get_local_peak_id(g,merge_to,jpeak)
                    current=sparse_saddle_dens%first(ipeak) ! first element of line ipeak
                    do while(current>0) ! walk the line
                        j=sparse_saddle_dens%col(current)
                        value_iij=sparse_saddle_dens%val(current) ! value of the matrix
                        ! Copy the value of density only if larger
                        if(value_iij>get_value(jpeak,j,sparse_saddle_dens))then
                            call set_value(jpeak,j,value_iij,sparse_saddle_dens)
                            call set_value(j,jpeak,value_iij,sparse_saddle_dens)
                        end if
                        current=sparse_saddle_dens%next(current)
                    end do
                    call set_value(jpeak,jpeak,zero,sparse_saddle_dens)
                endif
            endif
        enddo
        ! Set alive to zero for newly merged peaks
        nzero=0
        nsurvive=0
        do ipeak=1,npeaks
            if(c%alive(ipeak)>0)then
                merge_to=c%idc(ipeak)
                if(merge_to.NE.(npeak_cum(g%myid-1)+ipeak))then
                  c%alive(ipeak)=0
                  c%lev_peak(ipeak)=idepth
                  nzero=nzero+1
                else
                    nsurvive=nsurvive+1
                end if
            endif
        end do

        ! Remove all matrix elements corresponding to merged peaks
        do ipeak=1,hfree-1
            current=sparse_saddle_dens%first(ipeak) ! first element of line ipeak
            do while(current>0) ! walk the line
            j=sparse_saddle_dens%col(current)
            if(c%alive(ipeak)==0 .or. c%alive(j)==0)then
                call set_value(ipeak,j,zero,sparse_saddle_dens)
            endif
            current=sparse_saddle_dens%next(current)
            end do
        end do
        if(r%verbose)write(*,*)'level=',idepth,'nmove=',nzero,'survived=',nsurvive
        idepth=idepth+1
    end do
    ! End loop over peak levels
    mergelevel_max=idepth-2 ! last level has no more clumps, also idepth=idepth+1 still happens on last level.
    
    ! Compute maximum saddle density for each surviving clump
    ! Create new local duplicated peaks and update communicator
    do i=1,hfree-1
        call get_max(g,i,sparse_saddle_dens)
    end do

    ! Count surviving peaks
    nsurvive=0
    do ipeak=1,npeaks
        if(c%alive(ipeak)>0)then
            nsurvive=nsurvive+1
        endif
    end do

    if(g%myid==1)then
        if(action.EQ.'relevance')then
           write(*,*)'Found',nsurvive,' relevant peaks'
        endif
        if(action.EQ.'saddleden')then
           write(*,*)'Found',nsurvive,' halos'
        endif
     endif
    if(action.EQ.'relevance')then
        ! Compute relevance
        do ipeak=1,npeaks
            if(c%alive(ipeak)>0)then
            if (sparse_saddle_dens%maxval(ipeak)>0)then
                relevance_peak=c%denp(ipeak)/sparse_saddle_dens%maxval(ipeak)
            else
                relevance_peak=c%denp(ipeak)/density_threshold
            end if
            c%relevance(ipeak)=relevance_peak
            else
            c%relevance(ipeak)=0
            endif
        end do
        ! Merge all peaks to deepest level
        do ilev=idepth-2,0,-1
            do ipeak=1,npeaks
                if(c%lev_peak(ipeak)==ilev)then
                    merge_to=c%idc(ipeak)
                    call get_local_peak_id(g,merge_to,jpeak)
                    c%idc(ipeak)=c%idc(jpeak)
                endif
            end do
        end do

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

        ! Update flag2 field
        do ipart=1,p%npart
            dx_loc=r%boxlen/2**p%levelp(ipart)
            ckey(1:ndim)=int(p%xp(ipart,1:ndim)/dx_loc)
            hash_cell(0)=p%levelp(ipart)+1
            hash_cell(1:ndim)=ckey(1:ndim)
            call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icellp,flush_cache=.true.,fetch_cache=.true.)
            if (gridp%flag2(icellp)>0)then
                call get_local_peak_id(g,gridp%flag2(icellp),ipeak)
                merge_to=c%idc(ipeak)
                !call get_local_peak_id(g,merge_to,jpeak)
                gridp%flag2(icellp)=merge_to
                !update test particles peak
                p%pid(ipart) = merge_to
            end if
        end do



        call close_cache(s,m%grid_dict)
        
    endif

    if(action.EQ.'saddleden')then
        ! Compute peak index for the halo
        do ipeak=1,npeaks
            c%ind_halo(ipeak)=c%idc(ipeak)
        end do
        do ilev=idepth-2,0,-1
            do ipeak=1,npeaks
               if(c%lev_peak(ipeak)==ilev)then
                  merge_to=c%ind_halo(ipeak)
                  call get_local_peak_id(g,merge_to,jpeak)
                  c%ind_halo(ipeak)=c%ind_halo(jpeak)
               endif
            end do
        end do
        ! Compute halo masses
        halo_mass=0
        n_cells_halo=0
        do ipeak=1,npeaks
            merge_to=c%ind_halo(ipeak)
            call get_local_peak_id(g,merge_to,jpeak)
            c%halo_mass(jpeak)=c%halo_mass(jpeak)+c%clump_mass(ipeak)
            c%n_cell_halo(jpeak)=c%n_cell_halo(jpeak)+c%n_cells(ipeak)
        end do

        ! Assign back halo mass to peak
        do ipeak=1,npeaks
            merge_to=c%ind_halo(ipeak)
            call get_local_peak_id(g,merge_to,jpeak)
            c%halo_mass(ipeak)=c%halo_mass(jpeak)
        end do

    endif

    

end subroutine merge_clumps

!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_max(g,i,mat)
    use amr_commons, only: global_t
    use sparse_matrix
    use clfind_commons,ONLY: gkey,npeak_cum
    type(sparse_mat)::mat
    type(global_t)::g
    integer::i,npeaks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! get maximum in i-th line by walking the linked list
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer::current,icol

    npeaks = npeak_cum(g%myid)-npeak_cum(g%myid-1)
  
    mat%maxval(i)=0
    mat%maxloc(i)=0
  
    ! walk the line...
    current=mat%first(i)
    do  while( current /= 0 )
       if(mat%maxval(i)<mat%val(current))then
          mat%maxval(i)=mat%val(current)
          icol=mat%col(current)
          if(icol<=npeaks)then
             mat%maxloc(i)=npeak_cum(g%myid-1)+icol
          else
             mat%maxloc(i)=gkey(icol)
          endif
       end if
       current=mat%next(current)
    end do
  end subroutine get_max

  

!################################################################
!################################################################
!################################################################
!################################################################
subroutine write_clump_properties(r,g,c,to_file)
    use amr_commons
    !use pm_commons,ONLY:mp
    !use hydro_commons,ONLY:mass_sph
    use clfind_commons
    use amr_commons, only:run_t,global_t
    use mdl_module, only: mdl_mkdir, mdl_wtime
    implicit none
#ifndef WITHOUTMPI
    integer,parameter::tag=1101
    integer::dummy_io,info,info2
#endif
    logical::to_file
    type(run_t)::r
    type(global_t)::g
    type(clump_t)::c
    !---------------------------------------------------------------------------
    ! this routine writes the clump properties to screen and to file
    !---------------------------------------------------------------------------
  
    integer::i,j,jj,ilun,ilun2,n_rel,n_rel_tot,nx_loc,npeaks,IOGROUPSIZEREP,ifout
    real(dp)::rel_mass,rel_mass_tot,scale,particle_mass=0
    character(LEN=80)::fileloc,filedir
    character(LEN=5)::nchar,ncharcpu
    real(dp),dimension(1:c%npart)::peakd
    integer,dimension(1:c%npart)::ind_sort
  
#ifndef WITHOUTMPI
    real(dp)::particle_mass_tot
#endif
  
    if (.not. to_file)return
  
    !if(r%hydro)then
    !    particle_mass=MINVAL(mp, MASK=(mp > 0))
    !else
    !    particle_mass=mass_sph
    !endif
    particle_mass = 0.0
    npeaks = c%npart
    IOGROUPSIZEREP = 1
    ifout = 1
  
    if(to_file)then
       ilun=20
       ilun2=22
    else
       ilun=6
       ilun2=6
    end if
  
    ! print results in descending order to screen/file
    rel_mass=0
    n_rel=0
  
    if (to_file .eqv. .true.) then
       ! first create directories
       call title(ifout,nchar)
       filedir='output_'//TRIM(nchar)
       !call mdl_mkdir(mdl,filedir)
       !call create_output_dirs(filedir)
       ! Wait for the token
#ifndef WITHOUTMPI
       if(IOGROUPSIZE>0) then
          if (mod(g%myid-1,IOGROUPSIZE)/=0) then
             call MPI_RECV(dummy_io,1,MPI_INTEGER,g%myid-1-1,tag,&
                  & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
          end if
       endif
#endif
  
       if(IOGROUPSIZEREP>0)then
          call title(((g%myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
          fileloc=TRIM(filedir)//'/group_'//TRIM(ncharcpu)//'/clump_'//TRIM(nchar)//'.txt'
       else
          fileloc=TRIM(filedir)//'/clump_'//TRIM(nchar)//'.txt'
       endif
       call title(g%myid,nchar)
       fileloc=TRIM(fileloc)//TRIM(nchar)
       open(unit=ilun,file=fileloc,form='formatted')
  
       if(saddle_threshold>0)then
          call title(ifout,nchar)
          if(IOGROUPSIZEREP>0)then
             call title(((g%myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
             fileloc=TRIM(filedir)//'/group_'//TRIM(ncharcpu)//'/halo_'//TRIM(nchar)//'.txt'
          else
             fileloc=TRIM(filedir)//'/halo_'//TRIM(nchar)//'.txt'
          endif
          call title(g%myid,nchar)
          fileloc=TRIM(fileloc)//TRIM(nchar)
          open(unit=ilun2,file=fileloc,form='formatted')
       endif
    end if
  
    if (to_file .or. g%myid==1)then
       write(ilun,'(144A)')'   index  halo   lev   parent      ncell    peak_x             peak_y             peak_z     '//&
            '        rho-               rho+               rho_av             mass_cl            relevance   '
       if(saddle_threshold>0)then
          write(ilun2,'(135A)')'     index      ncell    peak_x             peak_y             peak_z     '//&
               '        rho+               mass      '
       endif
    end if
  
  
    do j=npeaks,1,-1
        jj=c%sortp(j)
        if (c%relevance(jj) > relevance_threshold .and. c%halo_mass(jj) > mass_threshold*particle_mass)then
        write(ilun,'(I8,X,I8,1X,I2,X,I10,X,I10,8(X,1PE18.9E2))')&
                jj+npeak_cum(g%myid-1)&
                ,c%ind_halo(jj)&
                ,c%lev_peak(jj)&
                ,c%idc(jj)&
                ,c%n_cells(jj)&
                ,c%peak_pos(jj,1)&
                ,c%peak_pos(jj,2)&
                ,c%peak_pos(jj,3)&
                ,c%min_dens(jj)&
                ,c%max_dens(jj)&
                ,c%clump_mass(jj)/c%clump_vol(jj)&
                ,c%clump_mass(jj)&
                ,c%relevance(jj)
        rel_mass=rel_mass+c%clump_mass(jj)
        n_rel=n_rel+1
        end if
        if(saddle_threshold>0)then
        if(c%ind_halo(jj).EQ.jj+npeak_cum(g%myid-1).AND.c%halo_mass(jj) > mass_threshold*particle_mass)then
            write(ilun2,'(I10,X,I10,5(X,1PE18.9E2))')&
            jj+npeak_cum(g%myid-1)&
                ,c%n_cell_halo(jj)&
                ,c%peak_pos(jj,1)&
                ,c%peak_pos(jj,2)&
                ,c%peak_pos(jj,3)&
                ,c%max_dens(jj)&
                ,c%halo_mass(jj)
        endif
        endif
    end do
  
    if (to_file)then
       close(ilun)
       if(saddle_threshold>0)then
          close(ilun2)
       endif
    end if
  
    ! Send the token
#ifndef WITHOUTMPI
    if(IOGROUPSIZE>0) then
        if(mod(g%myid,IOGROUPSIZE)/=0 .and.(g%myid.lt.g%ncpu))then
            dummy_io=1
            call MPI_SEND(dummy_io,1,MPI_INTEGER,g%myid-1+1,tag, &
                & MPI_COMM_WORLD,info2)
        end if
    endif
#endif
  
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(n_rel,n_rel_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    n_rel=n_rel_tot
    call MPI_ALLREDUCE(rel_mass,rel_mass_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    rel_mass=rel_mass_tot
#else
    n_rel_tot = n_rel
    rel_mass_tot = rel_mass
#endif
    if(g%myid==1)then
       if(clinfo)write(*,'(A,1PE12.5)')' Total mass [code units] above threshold =',tot_mass
       if(clinfo)write(*,'(A,I10,A,1PE12.5)')' Total mass [code units] in',n_rel_tot,' listed clumps =',rel_mass_tot
    endif
  
end subroutine write_clump_properties




end module clump_finder_module

