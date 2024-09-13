module clump_merger_module
contains
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_peak(s,global_peak_id,local_peak_id,flush_cache,fetch_cache,lock)
  use amr_commons
  use ramses_commons, only: ramses_t
  use clfind_commons
  use cache_commons
  use hash
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  integer(kind=8)::global_peak_id
  integer::local_peak_id
  logical::flush_cache
  logical::fetch_cache
  logical,optional::lock

  logical::failed_request
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::send_request_status_clump
  integer::info
#endif
  integer::i,iskip,ntile_response,icounter,peak_cpu
  integer::send_request_id_clump,response_id_clump
  integer(kind=8)::gpid
  integer::lpid
  
  associate(g=>s%g,c=>s%c,m=>s%m,mdl=>s%mdl)

#ifndef WITHOUTMPI
  ! If counter is good, check on incoming messages and perform actions
  if(mdl%mail_counter==32)then
     call check_mail(s,MPI_REQUEST_NULL,m%grid_dict)
     mdl%mail_counter=0
  endif
  mdl%mail_counter=mdl%mail_counter+1
#endif

  ! Peak is in the local processor memory
  if(    global_peak_id > c%npeak_cum(g%myid-1) .and. &
       & global_peak_id <= c%npeak_cum(g%myid) )then

     local_peak_id=global_peak_id-c%npeak_cum(g%myid-1)
     return

  else

#ifndef WITHOUTMPI
     ! Get position in cache memory from hash table
     local_peak_id = hash_getp_simple(c%peak_dict,global_peak_id)

     if(present(lock).and.local_peak_id>0)then
        if(lock)call lock_cache_clump(s,local_peak_id)
     endif

     ! If peak exists, all good
     if(local_peak_id > 0)return

     ! Get remote cpu to which peak belongs
     call get_global_peak_cpu(s,global_peak_id,peak_cpu)

     if(fetch_cache)then

        ! Send a request to the relevant cpu
        mdl%send_request_array_clump(1:2)=transfer(global_peak_id,local_peak_id,2)

        ! Post RECV for the expected response
        call MPI_IRECV(mdl%recv_fetch_array_clump,mdl%size_fetch_array_clump,MPI_INTEGER,peak_cpu-1,msg_tag_clump,MPI_COMM_WORLD,response_id_clump,info)

        ! Post SEND for the request
        call MPI_ISEND(mdl%send_request_array_clump,mdl%size_request_array_clump,MPI_INTEGER,peak_cpu-1,request_tag_clump,MPI_COMM_WORLD,send_request_id_clump,info)

        ! While waiting for reply, check on incoming messages and perform actions
        call check_mail(s,response_id_clump,m%grid_dict)

        ! Wait for ISEND completion to free memory in corresponding MPI buffer
        call MPI_WAIT(send_request_id_clump,send_request_status_clump,info)

        ! Check header for type of response
        iskip=1
        failed_request=(mdl%recv_fetch_array_clump(iskip).NE.1)
        iskip=iskip+1

        ! Number of tiles in the response buffer
        ntile_response=mdl%recv_fetch_array_clump(iskip)
        iskip=iskip+1

        ! Loop over tiles
        do i=1,ntile_response

           ! Find next available cache line
           if(c%locked(c%free_cache))then
              icounter=0
              do while(c%locked(c%free_cache))
                 c%free_cache=c%free_cache+1
                 icounter=icounter+1
                 if(c%free_cache>c%ncachemax)c%free_cache=1
                 if(icounter>c%ncachemax)then
                    write(*,*)'PE ',g%myid,'clump cache entirely locked'
                    stop
                 endif
              end do
           end if

           ! Next available peak in memory
           lpid=c%npeak+c%free_cache

           ! If cache line is occupied, free it.
           if(c%occupied(c%free_cache))call destage_clump(s,lpid,m%grid_dict)

           ! Get global peak id from message header
           gpid=transfer(mdl%recv_fetch_array_clump(iskip:iskip+1),gpid)
           iskip=iskip+2

           ! Insert local peak id in hash table
           call hash_setp_simple(c%peak_dict,gpid,lpid)

           c%gid(c%free_cache)=gpid
           c%parent_cpu(c%free_cache)=peak_cpu
           c%occupied(c%free_cache)=.true.
           c%dirty(c%free_cache)=.false.

           ! Set the request peak local peak id
           if(global_peak_id==gpid)local_peak_id=lpid

           ! Unpack response to fetch request
           call unpack_fetch_clump%proc(c,lpid,mdl%size_msg_array_clump,mdl%recv_fetch_array_clump(iskip:iskip+mdl%size_msg_array_clump-1))

           ! If we also have also a flush cache...
           ! This is for combined read-write cache operations
           if(flush_cache)then

              c%dirty(c%free_cache)=.true.

              ! Set initialisation rule for combiner operations
              call init_flush_clump%proc(c,lpid)

           endif

           ! Go to next free cache line
           c%free_cache=c%free_cache+1
           c%ncache=c%ncache+1
           if(c%free_cache.GT.c%ncachemax)then
              c%free_cache=1
           endif
           if(c%ncache.GT.c%ncachemax)c%ncache=c%ncachemax

           ! Go to next tile
           iskip=iskip+mdl%size_msg_array_clump

        end do

     else if(flush_cache)then

	! Find next available cache line
        if(c%locked(c%free_cache))then
           icounter=0
           do while(c%locked(c%free_cache))
              c%free_cache=c%free_cache+1
              icounter=icounter+1
              if(c%free_cache>c%ncachemax)c%free_cache=1
              if(icounter>c%ncachemax)then
                 write(*,*)'PE ',g%myid,'clump cache entirely locked'
                 stop
              endif
           end do
        end if

        ! Next available peak in memory
        local_peak_id=c%npeak+c%free_cache

        ! If cache line is occupied, free it.
        if(c%occupied(c%free_cache))call destage_clump(s,local_peak_id,m%grid_dict)

        ! Insert local peak id in hash table
        call hash_setp_simple(c%peak_dict,global_peak_id,local_peak_id)

        c%gid(c%free_cache)=global_peak_id
        c%parent_cpu(c%free_cache)=peak_cpu
        c%occupied(c%free_cache)=.true.
        c%dirty(c%free_cache)=.true.

        ! Set initialisation rule for combiner operations
        call init_flush_clump%proc(c,local_peak_id)

        ! Go to next free cache line
        c%free_cache=c%free_cache+1
        c%ncache=c%ncache+1
        if(c%free_cache.GT.c%ncachemax)then
           c%free_cache=1
        endif
        if(c%ncache.GT.c%ncachemax)c%ncache=c%ncachemax

     endif

     if(present(lock).and.local_peak_id>0)then
        if(lock)call lock_cache_clump(s,local_peak_id)
     endif

#endif
  endif

  end associate

end subroutine get_peak
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine lock_cache_clump(s,local_peak_id)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  integer::local_peak_id
  !
  ! This routine locks a cache line because
  ! it will be updated later.
  !
  integer::icache
  if(local_peak_id<=s%c%npeak)return
  icache=local_peak_id-s%c%npeak
  s%c%locked(icache)=.true.
end subroutine lock_cache_clump
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine unlock_cache_clump(s,local_peak_id)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  integer::local_peak_id
  !
  ! This routine unlocks a cache line because
  ! it has been updated and can be flushed.
  !
  integer::icache
  if(local_peak_id<=s%c%npeak)return
  icache=local_peak_id-s%c%npeak
  s%c%locked(icache)=.false.
end subroutine unlock_cache_clump
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_global_peak_cpu(s,global_peak_id,peak_cpu)
  use amr_commons
  use clfind_commons
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  integer(kind=8)::global_peak_id
  integer::peak_cpu
  ! get the mpi-domain a peak with input global id
  integer::icpu

  associate(g=>s%g,c=>s%c)

    peak_cpu = g%ncpu
    do icpu = 1,g%ncpu
       if(    global_peak_id > c%npeak_cum(icpu-1) .and. &
            & global_peak_id <= c%npeak_cum(icpu))then
          peak_cpu = icpu
       endif
    end do

  end associate

end subroutine get_global_peak_cpu
!################################################################
!################################################################
!################################################################
!################################################################
subroutine allocate_peak_patch_arrays(s)
  use amr_parameters, ONLY: ndim, dp, nbin
  use clfind_commons
  use ramses_commons, ONLY: ramses_t
  implicit none
  type(ramses_t)::s

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c)

  !------------------------------------------
  ! Compute the max size of peak-based arrays
  !------------------------------------------
  c%npeak_max=c%npeak+c%ncachemax

  !-------------------------------
  ! Allocate peak-patch properties
  !-------------------------------
  allocate(c%saddle_dens(1:c%npeak_max))
  allocate(c%saddle_nbor(1:c%npeak_max))

  allocate(c%lev_peak(1:c%npeak_max))
  allocate(c%new_peak(c%npeak_max))
  allocate(c%relevance(1:c%npeak_max))

  allocate(c%peak_pos(1:c%npeak_max,1:ndim))
  allocate(c%peak_vel(1:c%npeak_max,1:ndim))
  allocate(c%peak_acc(1:c%npeak_max,1:ndim))

  allocate(c%min_dens(1:c%npeak_max))
  allocate(c%n_cells(1:c%npeak_max))
  allocate(c%clump_mass(1:c%npeak_max))
  allocate(c%clump_vol(1:c%npeak_max))

  allocate(c%particle_mass(1:c%npeak_max))

  !-------------------------------
  ! Allocate halo-patch properties
  !-------------------------------
  allocate(c%halo_mass(1:c%npeak_max))
  allocate(c%ind_halo(1:c%npeak_max))
  allocate(c%n_cells_halo(1:c%npeak_max))

  allocate(c%max_peak_mass(1:c%npeak_max))
  allocate(c%ind_max_mass(1:c%npeak_max))
  allocate(c%ind_halo_1(1:c%npeak_max))
  allocate(c%ind_halo_2(1:c%npeak_max))
  allocate(c%ind_halo_3(1:c%npeak_max))
  allocate(c%ind_central(1:c%npeak_max))

  allocate(c%mass_bin(1:c%npeak_max,1:nbin))

  !-----------------------------------
  ! Allocate sink particles properties
  !-----------------------------------
  if(r%sink)then
     allocate(c%occupied_sink(1:c%npeak_max))
     allocate(c%form_sink(1:c%npeak_max))
  endif

  !--------------------
  ! Allocate hash table
  !--------------------
  call init_empty_hash_simple(c%peak_dict,c%ncachemax)

  !----------------------
  ! Allocate cache memory
  !----------------------
  allocate(c%dirty(1:c%ncachemax))
  allocate(c%locked(1:c%ncachemax))
  allocate(c%occupied(1:c%ncachemax))
  allocate(c%parent_cpu(1:c%ncachemax))
  allocate(c%gid(1:c%ncachemax))
  c%dirty=.false.
  c%locked=.false.
  c%occupied=.false.
  c%free_cache=1
  c%ncache=0

  !---------------------
  ! Allocate cache comms
  !---------------------
  call init_cache_clump(s%mdl)

  !------------------------------------------------
  ! Initialize all peak based arrays for clump finder
  !------------------------------------------------
  c%lev_peak=0; c%new_peak=0; c%ind_halo=0; c%relevance=1

  end associate

end subroutine allocate_peak_patch_arrays
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_deallocate_clump(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CLUMP_DEALLOC,pst%iUpper+1,input_size,0,ilevel)
     call r_deallocate_clump(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call deallocate_peak_patch_arrays(pst%s)
  endif

end subroutine r_deallocate_clump
!################################################################
!################################################################
!################################################################
!################################################################
subroutine deallocate_peak_patch_arrays(s)
  use clfind_commons
  use ramses_commons, only: ramses_t
  use hash, only: reset_entire_hash_simple
  implicit none
  type(ramses_t)::s

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c)

  ! Deallocate test particle arrays
  if(c%ntest>0)then
     deallocate(c%cell)
     deallocate(c%grid)
     deallocate(c%level)
     deallocate(c%hash)
  endif
  if(c%ntest_tot==0)return

  ! Deallocate cumulative peak count CPU array
  deallocate(c%npeak_cum)
  if(c%npeak_tot==0)return

  ! Deallocate peak-patch arrays
  deallocate(c%peak_cell)
  deallocate(c%peak_grid)
  deallocate(c%peak_level)
  deallocate(c%max_dens)

  deallocate(c%saddle_dens)
  deallocate(c%saddle_nbor)

  deallocate(c%lev_peak)
  deallocate(c%new_peak)
  deallocate(c%relevance)

  deallocate(c%peak_pos)
  deallocate(c%peak_vel)
  deallocate(c%peak_acc)

  deallocate(c%min_dens)
  deallocate(c%n_cells)
  deallocate(c%clump_mass)
  deallocate(c%clump_vol)
  deallocate(c%particle_mass)

  ! Deallocate halo-patch arrays
  deallocate(c%halo_mass)
  deallocate(c%ind_halo)
  deallocate(c%n_cells_halo)

  deallocate(c%mass_bin)

  deallocate(c%max_peak_mass)
  deallocate(c%ind_max_mass)
  deallocate(c%ind_halo_1)
  deallocate(c%ind_halo_2)
  deallocate(c%ind_halo_3)
  deallocate(c%ind_central)

  if(r%sink)then
     deallocate(c%occupied_sink)
     deallocate(c%form_sink)
  endif

  ! Deallocate hash table
  call reset_entire_hash_simple(c%peak_dict)

  ! Dellocate cache memory
  deallocate(c%dirty)
  deallocate(c%locked)
  deallocate(c%occupied)
  deallocate(c%parent_cpu)
  deallocate(c%gid)

  ! Deallocate cache comms
  call kill_cache_clump(s%mdl)

  end associate

end subroutine deallocate_peak_patch_arrays
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine collect_saddle(s)
  use amr_parameters, only: twotondim,ndim
  use amr_commons, only:oct,nbor
  use ramses_commons, only: ramses_t
  use cache_commons
  use cache
  use nbors_utils
  implicit none
  type(ramses_t)::s
  !===================================================================
  ! This is the clump finder routine for collecting the densest saddle
  ! points between each patch and its corresponding neighbor.
  ! Written by Ziyong Wu (mini-ramses version December 2023).
  !==================================================================
  type(msg_int4_small_realdp)::dummy_int4_small_realdp
  type(msg_saddle_clump)::dummy_saddle_clump
  type(oct),pointer::gridn
  integer:: ilevel
  integer::icpu,next_level,now_level,icelln,idim,j,jpeak,k
  integer::ipart,jpart,ip,i,icellp,icellpm,ipeak,itest,igrid,ind,peak_cen,peak_nbor
  integer(kind=8),dimension(1:s%g%ncpu)::npeak_cpu,npeak_cpu_all
  integer,dimension(1:ndim)::ckey,ckey_nbor
  integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
  real(dp)::dens_cen,dens_ave,dens_nbor,x,y,z
  real(dp),dimension(1:ndim)::xcen,xnei
  integer, parameter::nSnei=56
  real(dp),dimension(1:ndim,1:nSnei)::xSnei
  type(nbor),dimension(1:nSnei) :: grid_nbor
  integer(kind=8),dimension(1:nSnei)::icell_nbor
  integer(kind=8)::global_peak_id
  logical::ok

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c)

  !--------------------------------------------------------
  ! Arrays to define neighbors (center=[0,0,0])
  ! normalized to dx = 1 = size of the central leaf cell
  ! from -0.75 to 0.75
  !--------------------------------------------------------
  ind=0
  do k=1,4
     do j=1,4
        do i=1,4
           ok=.true.
           if((i==2.or.i==3).and.(j==2.or.j==3).and.(k==2.or.k==3)) ok=.false. ! centre
           if(ok)then
              ind = ind+1
              x = (i-1)+0.5d0 - 2
              y = (j-1)+0.5d0 - 2
              z = (k-1)+0.5d0 - 2
              xSnei(1,ind) = x/2d0
#if NDIM>1
              xSnei(2,ind) = y/2d0
#endif
#if NDIM>2
              xSnei(3,ind) = z/2d0
#endif
           endif
        enddo
     enddo
  enddo

  !----------------------------------------------------
  ! Compute densest saddle point and associated peak id
  !----------------------------------------------------
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
       hilbert=m%domain,pack_size=storage_size(dummy_int4_small_realdp)/32,&
       pack=pack_fetch_saddle, unpack=unpack_fetch_saddle)

  call open_cache_clump(s,pack_size=storage_size(dummy_saddle_clump)/32,&
       init=init_flush_saddle,flush=pack_flush_saddle,combine=unpack_flush_saddle)

  c%saddle_dens = 0
  c%saddle_nbor = 0

  do itest=1,c%ntest
     ilevel=c%level(itest)
     igrid=c%grid(itest)
     ind=c%cell(itest)

     peak_cen = m%grid(igrid)%flag1(ind)
#ifdef GRAV
     dens_cen = m%grid(igrid)%rho(ind)
#endif
     ! Set pointers to null
     icelln=0
     nullify(gridn)
     do j=1,nSnei
        nullify(grid_nbor(j)%p)
     end do

     xcen(1)=2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5
#if NDIM>1
     xcen(2)=2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5
#endif
#if NDIM>2
     xcen(3)=2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5
#endif
     ! Collect all neighboring cell from hash table
     do j=1,nSnei

        ! Compute neighboring cell coordinates
        xnei(1:ndim)=xcen(1:ndim)+xSnei(1:ndim,j)
        ! Periodic boundary conditions
        do idim=1,ndim
           if(xnei(idim)<                0.0d0)xnei(idim)=xnei(idim)+m%ckey_max(ilevel+1)
           if(xnei(idim)>=m%ckey_max(ilevel+1))xnei(idim)=xnei(idim)-m%ckey_max(ilevel+1)
        end do

        ! Get neighboring cell at ilevel
        ckey_nbor(1:ndim)=int(xnei(1:ndim))
        hash_nbor(0)=ilevel+1
        hash_nbor(1:ndim)=ckey_nbor(1:ndim)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.false.,fetch_cache=.true.)

        ! If missing, get neighboring cell at ilevel-1
        if(.not.associated(gridn))then
           ckey_nbor(1:ndim)=int(xnei(1:ndim)/2.0)
           hash_nbor(0)=ilevel
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.false.,fetch_cache=.true.)

           ! If refined, get neighboring cell at ilevel+1
        else if (gridn%refined(icelln))then
           ckey_nbor(1:ndim)=int(xnei(1:ndim)*2.0)
           hash_nbor(0)=ilevel+2
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.false.,fetch_cache=.true.)
        endif

        ! Lock grid in cache
        call lock_cache(s,gridn)

        grid_nbor(j)%p => gridn
        icell_nbor(j) = icelln

     end do

     do j=1,nSnei
        gridn => grid_nbor(j)%p ! Gather neighboring grid
        icelln = icell_nbor(j)

        peak_nbor = gridn%flag1(icelln)
#ifdef GRAV
        dens_nbor = gridn%rho(icelln)
#endif
        ok = peak_cen/=0
        ok = ok .and. peak_nbor/=0
        ok = ok .and. peak_cen/=peak_nbor

        ! If saddle density is larger, set new densest saddle point
        if(ok)then
           dens_ave = 0.5*(dens_cen+dens_nbor)
           global_peak_id = peak_cen
           call get_peak(s,global_peak_id,ipeak,flush_cache=.true.,fetch_cache=.false.)
           if(dens_ave>c%saddle_dens(ipeak))then
              c%saddle_dens(ipeak)=dens_ave
              c%saddle_nbor(ipeak)=peak_nbor
           end if
        endif
     end do

     ! Unlock neighboring grids
     do j=1,nSnei
        gridn => grid_nbor(j)%p
        call unlock_cache(s,gridn)
     end do

  end do

  call close_cache(s,m%grid_dict)

  end associate

end subroutine collect_saddle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_saddle(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_int4_small_realdp)::msg

  do ind=1,twotondim
     msg%flg(ind)=grid%flag1(ind)
  end do
  do ind=1,twotondim
     if(grid%refined(ind))then
        msg%ref(ind)=1
     else
        msg%ref(ind)=0
     endif
  end do
#ifdef GRAV
  do ind=1,twotondim
     msg%realdp(ind)=grid%rho(ind)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_saddle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_saddle(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_int4_small_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     grid%flag1(ind)=msg%flg(ind)
  end do
  do ind=1,twotondim
     if(msg%ref(ind)==1)then
        grid%refined(ind)=.true.
     else
        grid%refined(ind)=.false.
     endif
  end do
#ifdef GRAV
  do ind=1,twotondim
     grid%rho(ind)=msg%realdp(ind)
  end do
#endif

end subroutine unpack_fetch_saddle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_saddle(c,local_peak_id)
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  c%saddle_dens(local_peak_id)=0d0
  c%saddle_nbor(local_peak_id)=0

end subroutine init_flush_saddle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_saddle(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_saddle_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_saddle_clump)::msg

  msg%nbor=c%saddle_nbor(local_peak_id)
  msg%dens=c%saddle_dens(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_saddle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_saddle(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_saddle_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_saddle_clump)::msg

  msg=transfer(msg_array,msg)

  if(msg%dens>c%saddle_dens(local_peak_id))then
     c%saddle_dens(local_peak_id)=msg%dens
     c%saddle_nbor(local_peak_id)=msg%nbor
  endif

end subroutine unpack_flush_saddle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine merge_clumps(s,action)
  use amr_commons, only: dp, ndim
  use ramses_commons, only: ramses_t
  use cache_commons, only: msg_merge_clump, msg_halo_clump
  use cache
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  character(len=9)::action
  !---------------------------------------------------------
  ! This routine merges the irrelevant clumps
  ! - clumps are sorted by ascending max density
  ! - irrelevent clumps are merged to most relevant neighbor
  !---------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::j,i,ipart,igrid,ind,itest
  integer::current,nmove,ipeak,jpeak,iter
  integer::nsurvive,nzero,idepth
  integer::ilev,ilevel,mergelevel_max
  integer(kind=8)::global_peak_id,merge_to
  real(dp)::value_iij,zero=0,relevance_peak
  real(dp)::d,dx_loc,vol
  integer,dimension(1:s%c%npeak_max)::alive,ind_sort
  real(dp),dimension(1:s%c%npeak_max)::peakd
  logical::do_merge=.false.
  type(msg_merge_clump)::dummy_merge_clump
  type(msg_halo_clump)::dummy_halo_clump
#ifndef WITHOUTMPI
  integer::mergelevel_max_global
  integer::nmove_all,nsurvive_all,nzero_all
#endif

  associate(g=>s%g,r=>s%r,m=>s%m,c=>s%c)

  if (r%verbose.and.g%myid==1)then
     if(action.EQ.'relevance')then
        write(*,*)'Now merging irrelevant clumps'
     endif
     if(action.EQ.'saddleden')then
        write(*,*)'Now merging clumps into halos'
     endif
  endif

  ! Initialize new_peak array to global peak id
  ! All peaks are alive at the start
  do i=1,c%npeak
     c%new_peak(i)=c%npeak_cum(g%myid-1)+i
     if(action.EQ.'relevance')then
        alive(i)=1
     endif
     if(action.EQ.'saddleden')then
        if(c%relevance(i)>c%relevance_threshold)then
           alive(i)=1
        else
           alive(i)=0
        endif
     endif
  end do

  ! Sort peaks by maximum peak density in ascending order
  do i=1,c%npeak
     peakd(i)=c%max_dens(i)
     ind_sort(i)=i
  end do
  call quick_sort_dp(peakd,ind_sort,c%npeak)

  ! Loop over peak levels
  nzero=c%npeak_tot
  idepth=0
  do while(nzero>0)

     ! Merge peaks
     nmove=c%npeak_tot
     iter=0
     do while(nmove>0)

        call open_cache_clump(s,pack_size=storage_size(dummy_merge_clump)/32,&
             pack=pack_fetch_merge,unpack=unpack_fetch_merge)

        nmove=0
        do i=c%npeak,1,-1
           ipeak=ind_sort(i)
           merge_to=c%new_peak(ipeak)
           if(alive(ipeak)>0)then
              if(action.EQ.'relevance')then
                 if(c%saddle_dens(ipeak)>0)then
                    relevance_peak=c%max_dens(ipeak)/c%saddle_dens(ipeak)
                 else
                    relevance_peak=c%max_dens(ipeak)/c%density_threshold
                 end if
                 do_merge=relevance_peak<c%relevance_threshold
              endif
              if(action.EQ.'saddleden')then
                 do_merge=(c%saddle_dens(ipeak)>c%saddle_threshold)
              endif
              if(do_merge)then
                 if(c%saddle_nbor(ipeak)>0)then
                    global_peak_id=c%saddle_nbor(ipeak)
                    call get_peak(s,global_peak_id,jpeak,flush_cache=.false.,fetch_cache=.true.)
                    if(c%max_dens(jpeak)>c%max_dens(ipeak))then
                       merge_to=c%new_peak(jpeak)
                    else if(c%max_dens(jpeak)==c%max_dens(ipeak))then
                       merge_to=MIN(c%new_peak(ipeak),c%new_peak(jpeak))
                    endif
                 endif
              endif
           endif
           if(c%new_peak(ipeak).NE.merge_to)then
              nmove=nmove+1
              c%new_peak(ipeak)=merge_to
           endif
        end do

        call close_cache(s,m%grid_dict)

        iter=iter+1
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(nmove,nmove_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        nmove=nmove_all
#endif
        if(g%myid==1.and.r%verbose)write(*,*)'niter=',iter,'nmove=',nmove
     end do

     ! Update flag1 field
     call open_cache_clump(s,pack_size=storage_size(dummy_merge_clump)/32,&
          pack=pack_fetch_merge,unpack=unpack_fetch_merge)
     do itest=1,c%ntest
        igrid=c%grid(itest)
        ind=c%cell(itest)
        global_peak_id=m%grid(igrid)%flag1(ind)
        if (global_peak_id>0)then
           call get_peak(s,global_peak_id,ipeak,flush_cache=.false.,fetch_cache=.true.)
           m%grid(igrid)%flag1(ind)=c%new_peak(ipeak)
        end if
     end do
     call close_cache(s,m%grid_dict)

     ! Compute new saddle points and corresponding nboring peaks.
     call collect_saddle(s)

     ! Set alive to zero for newly merged peaks
     nzero=0
     nsurvive=0
     do ipeak=1,c%npeak
        if(alive(ipeak)>0)then
           merge_to=c%new_peak(ipeak)
           if(merge_to.NE.(c%npeak_cum(g%myid-1)+ipeak))then
              alive(ipeak)=0
              c%lev_peak(ipeak)=idepth
              nzero=nzero+1
           else
              nsurvive=nsurvive+1
           end if
        endif
     end do

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(nzero,nzero_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     nzero=nzero_all
     call MPI_ALLREDUCE(nsurvive,nsurvive_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     nsurvive=nsurvive_all
#endif
     if(r%verbose.and.g%myid==1)write(*,*)'level=',idepth,'nmove=',nzero,'survived=',nsurvive
     idepth=idepth+1

  end do
  ! End loop over peak levels

  mergelevel_max=idepth-2 ! last level has no more clumps, also idepth=idepth+1 still happens on last level.
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(mergelevel_max,mergelevel_max_global,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  mergelevel_max=mergelevel_max_global
#endif

  ! Count surviving peaks
  nsurvive=0
  do ipeak=1,c%npeak
     if(alive(ipeak)>0)then
        nsurvive=nsurvive+1
     endif
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nsurvive,nsurvive_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nsurvive=nsurvive_all
#endif
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
     do ipeak=1,c%npeak
        if(alive(ipeak)>0)then
           if (c%saddle_dens(ipeak)>0)then
              relevance_peak=c%max_dens(ipeak)/c%saddle_dens(ipeak)
           else
              relevance_peak=c%max_dens(ipeak)/c%density_threshold
           end if
           c%relevance(ipeak)=relevance_peak
        else
           c%relevance(ipeak)=0
        endif
     end do

     ! Merge all peaks to deepest level
     do ilev=idepth-2,0,-1
        call open_cache_clump(s,pack_size=storage_size(dummy_merge_clump)/32,&
             pack=pack_fetch_merge,unpack=unpack_fetch_merge)
        do ipeak=1,c%npeak
           if(c%lev_peak(ipeak)==ilev)then
              merge_to=c%new_peak(ipeak)
              call get_peak(s,merge_to,jpeak,flush_cache=.false.,fetch_cache=.true.)
              c%new_peak(ipeak)=c%new_peak(jpeak)
           endif
        end do
        call close_cache(s,m%grid_dict)
     end do

  endif

  if(action.EQ.'saddleden')then

     ! Compute peak index for the halo
     do ipeak=1,c%npeak
        c%ind_halo(ipeak)=c%new_peak(ipeak)
     end do

     do ilev=idepth-2,0,-1
        call open_cache_clump(s,pack_size=storage_size(dummy_merge_clump)/32,&
             pack=pack_fetch_merge,unpack=unpack_fetch_merge)
        do ipeak=1,c%npeak
           if(c%lev_peak(ipeak)==ilev)then
              merge_to=c%ind_halo(ipeak)
              call get_peak(s,merge_to,jpeak,flush_cache=.false.,fetch_cache=.true.)
              c%ind_halo(ipeak)=c%ind_halo(jpeak)
           endif
        end do
        call close_cache(s,m%grid_dict)
     end do

     ! Compute halo masses
     c%halo_mass=0
     c%n_cells_halo=0
     call open_cache_clump(s,pack_size=storage_size(dummy_halo_clump)/32,&
          init=init_flush_halo,flush=pack_flush_halo,combine=unpack_flush_halo)
     do ipeak=1,c%npeak
        merge_to=c%ind_halo(ipeak)
        call get_peak(s,merge_to,jpeak,flush_cache=.true.,fetch_cache=.false.)
        c%halo_mass(jpeak)=c%halo_mass(jpeak)+c%clump_mass(ipeak)
        c%n_cells_halo(jpeak)=c%n_cells_halo(jpeak)+c%n_cells(ipeak)
     end do
     call close_cache(s,m%grid_dict)

     ! Assign back halo mass to peak halo mass
     call open_cache_clump(s,pack_size=storage_size(dummy_halo_clump)/32,&
          pack=pack_fetch_halo,unpack=unpack_fetch_halo)
     do ipeak=1,c%npeak
        merge_to=c%ind_halo(ipeak)
        call get_peak(s,merge_to,jpeak,flush_cache=.false.,fetch_cache=.true.)
        c%halo_mass(ipeak)=c%halo_mass(jpeak)
     end do
     call close_cache(s,m%grid_dict)

  endif

  end associate

end subroutine merge_clumps
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_merge(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_merge_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_merge_clump)::msg

  msg%npeak=c%new_peak(local_peak_id)
  msg%nhalo=c%ind_halo(local_peak_id)
  msg%mdens=c%max_dens(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_merge
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_merge(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_merge_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_merge_clump)::msg

  msg=transfer(msg_array,msg)

  c%max_dens(local_peak_id)=msg%mdens
  c%new_peak(local_peak_id)=msg%npeak
  c%ind_halo(local_peak_id)=msg%nhalo

end subroutine unpack_fetch_merge
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_halo(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_halo_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_halo_clump)::msg

  msg%ncell=c%n_cells_halo(local_peak_id)
  msg%mhalo=c%halo_mass(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_halo
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_halo(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_halo_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_halo_clump)::msg

  msg=transfer(msg_array,msg)

  c%n_cells_halo(local_peak_id)=msg%ncell
  c%halo_mass(local_peak_id)=msg%mhalo

end subroutine unpack_fetch_halo
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_halo(c,local_peak_id)
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  c%n_cells_halo(local_peak_id)=0
  c%halo_mass(local_peak_id)=0d0

end subroutine init_flush_halo
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_halo(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_halo_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_halo_clump)::msg

  msg%ncell=c%n_cells_halo(local_peak_id)
  msg%mhalo=c%halo_mass(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_halo
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_halo(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_halo_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_halo_clump)::msg

  msg=transfer(msg_array,msg)

  c%halo_mass(local_peak_id)=c%halo_mass(local_peak_id)+msg%mhalo
  c%n_cells_halo(local_peak_id)=c%n_cells_halo(local_peak_id)+msg%ncell

end subroutine unpack_flush_halo
!################################################################
!################################################################
!################################################################
!################################################################
subroutine compute_clump_properties(s,rtype)
  use amr_commons, only: dp,ndim
  use clfind_commons
  use ramses_commons, only: ramses_t
  use cache_commons, only: msg_prop_clump
  use cache
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  type(ramses_t)::s
  integer::rtype
  !----------------------------------------------------------------------------
  ! This subroutine performs a loop over all cells above the threshold and
  ! collects the  relevant information. After some MPI communications,
  ! all necessary peak-patch properties are computed
  !----------------------------------------------------------------------------
  type(msg_prop_clump)::dummy_prop_clump
  integer(kind=8)::global_peak_id
  integer::ipart,grid,peak_nr,ilevel,ipeak,plevel,igrid,itest,icelln,idim,ind
  real(dp),dimension(1:ndim)::xcell,accel
  real(dp)::dx_loc,tot_mass
  real(dp)::zero=0
  ! variables needed temporarily store cell properties
  real(dp)::d=0, vol=0
  ! variables related to the size of a cell on a given level
  integer::nx_loc
  logical::periodic
#ifndef WITHOUTMPI
  integer::i
  real(dp)::tot_mass_tot
#endif

  associate(g=>s%g,r=>s%r,m=>s%m,c=>s%c)

  if(g%myid==1.and.r%verbose)write(*,*)'Entering compute clump properties'

  !-----------------------------------------------------------------------
  ! Loop over local peaks and compute peak cell coordinates, velocities...
  !-----------------------------------------------------------------------
  c%peak_pos=0d0; c%peak_vel=0d0; c%peak_acc=0d0

  do ipeak=1,c%npeak
     ilevel=c%peak_level(ipeak)
     igrid=c%peak_grid(ipeak)
     ind=c%peak_cell(ipeak)
     dx_loc=r%boxlen/2**ilevel
     ! Peak cell coordinates and acceleration
     xcell(1)=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx_loc-m%skip(1)
#if NDIM>1
     xcell(2)=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx_loc-m%skip(2)
#endif
#if NDIM>2
     xcell(3)=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx_loc-m%skip(3)
#endif
     c%peak_pos(ipeak,1:ndim)=xcell(1:ndim)
#ifdef HYDRO
     if (r%hydro.AND.rtype.eq.4)then
        c%peak_vel(ipeak,1:ndim)=m%grid(igrid)%uold(ind,2:ndim+1)/m%grid(igrid)%uold(ind,1)
     endif
#endif
#ifdef GRAV
     c%peak_acc(ipeak,1:ndim)=m%grid(igrid)%f(ind,1:ndim)
#endif
  end do

  !--------------------------------------------------------------------------
  ! Loop over all cells above the threshold and compute mass, volume...
  !--------------------------------------------------------------------------
  c%min_dens=huge(zero); c%n_cells=0
  c%halo_mass=0d0; c%clump_mass=0d0; c%clump_vol=0d0

  call open_cache_clump(s,storage_size(dummy_prop_clump)/32,&
       init=init_flush_prop,flush=pack_flush_prop,combine=unpack_flush_prop)
  do itest=1,c%ntest
     ilevel=c%level(itest)
     igrid=c%grid(itest)
     ind=c%cell(itest)
     global_peak_id=m%grid(igrid)%flag1(ind)

     ! Save peak patch id into flag2 because flag1 will become halo patch id
     m%grid(igrid)%flag2(ind)=m%grid(igrid)%flag1(ind)

     if (global_peak_id /=0 ) then
        call get_peak(s,global_peak_id,peak_nr,flush_cache=.true.,fetch_cache=.false.)
        
        ! Cell density
#ifdef GRAV
        d=m%grid(igrid)%rho(ind)
#endif
        ! Cell volume
        dx_loc=r%boxlen/2**ilevel
        vol=dx_loc**ndim
        
        ! Number of leaf cells per clump
        c%n_cells(peak_nr)=c%n_cells(peak_nr)+1
        
        ! Clump min density
        c%min_dens(peak_nr)=min(c%min_dens(peak_nr),d)
        
        ! Clump mass
        c%clump_mass(peak_nr)=c%clump_mass(peak_nr)+vol*d
        
        ! Clump volume
        c%clump_vol(peak_nr)=c%clump_vol(peak_nr)+vol
        
     end if
  end do
  call close_cache(s,m%grid_dict)

  ! Initialize halo mass to clump mass
  c%halo_mass(1:c%npeak)=c%clump_mass(1:c%npeak)

  ! Calculate total mass above threshold
  tot_mass=sum(c%clump_mass(1:c%npeak))

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(tot_mass,tot_mass_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  tot_mass=tot_mass_tot
#endif

  end associate

end subroutine compute_clump_properties
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_prop(c,local_peak_id)
  use amr_commons, only: dp,ndim
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id
  real(dp)::zero

  c%n_cells(local_peak_id)=0
  c%min_dens(local_peak_id)=huge(zero)
  c%clump_mass(local_peak_id)=0d0
  c%clump_vol(local_peak_id)=0d0

end subroutine init_flush_prop
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_prop(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg%ncell=c%n_cells(local_peak_id)
  msg%dens=c%min_dens(local_peak_id)
  msg%mass=c%clump_mass(local_peak_id)
  msg%vol=c%clump_vol(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_prop
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_prop(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg=transfer(msg_array,msg)

  c%n_cells(local_peak_id)=c%n_cells(local_peak_id)+msg%ncell
  c%min_dens(local_peak_id)=min(c%min_dens(local_peak_id),msg%dens)
  c%clump_mass(local_peak_id)=c%clump_mass(local_peak_id)+msg%mass
  c%clump_vol(local_peak_id)=c%clump_vol(local_peak_id)+msg%vol

end subroutine unpack_flush_prop
!################################################################
!################################################################
!################################################################
!################################################################
subroutine trim_clumps(s)
  use amr_commons, only: dp,ndim
  use clfind_commons
  use ramses_commons, only: ramses_t
  use cache_commons, only: msg_prop_clump
  use cache
  implicit none
  type(ramses_t)::s
  !----------------------------------------------------------------------------
  ! This subroutine remove all clumps and halosthat are considered irrelevant.
  ! They are removed because their relevance (or peakiness) is below the
  ! relevance threshold or because their mass is too small.
  ! The flag1 and flag2 arrays are modified accordingly.
  !----------------------------------------------------------------------------
  type(msg_prop_clump)::dummy_prop_clump
  integer(kind=8)::global_peak_id,global_halo_id
  integer::ipeak,jpeak,igrid,ilevel,ind,itest

  associate(g=>s%g,r=>s%r,m=>s%m,c=>s%c)

  if(g%myid==1.and.r%verbose)write(*,*)'Entering trim_clumps'

  call open_cache_clump(s,storage_size(dummy_prop_clump)/32,&
       pack=pack_fetch_prop,unpack=unpack_fetch_prop)
  do itest=1,c%ntest
     ilevel=c%level(itest)
     igrid=c%grid(itest)
     ind=c%cell(itest)
     global_peak_id=m%grid(igrid)%flag2(ind)
     if (global_peak_id /=0 ) then
        call get_peak(s,global_peak_id,ipeak,flush_cache=.false.,fetch_cache=.true.)
        if(c%relevance(ipeak).LE.c%relevance_threshold)then
           m%grid(igrid)%flag2(ind)=0
        endif
        if(c%clump_mass(ipeak).LE.c%mass_threshold)then
           m%grid(igrid)%flag2(ind)=0
        endif
     endif
     if(c%saddle_threshold>0)then
        global_halo_id=m%grid(igrid)%flag1(ind)
        if (global_halo_id /=0 ) then
           call get_peak(s,global_halo_id,ipeak,flush_cache=.false.,fetch_cache=.true.)
           if(c%relevance(ipeak).LE.c%relevance_threshold)then
              m%grid(igrid)%flag1(ind)=0
           endif
           if(c%halo_mass(ipeak).LE.c%mass_threshold)then
              m%grid(igrid)%flag1(ind)=0
           endif
        endif
     endif
  end do
  call close_cache(s,m%grid_dict)

  end associate

end subroutine trim_clumps
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_prop(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg%dens=c%relevance(local_peak_id)
  msg%mass=c%clump_mass(local_peak_id)
  msg%vol=c%halo_mass(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_prop
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_prop(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg=transfer(msg_array,msg)

  c%relevance(local_peak_id)=msg%dens
  c%clump_mass(local_peak_id)=msg%mass
  c%halo_mass(local_peak_id)=msg%vol

end subroutine unpack_fetch_prop
!################################################################
!################################################################
!################################################################
!################################################################
subroutine central_in_halos(s)
  use amr_commons, only: dp,ndim
  use clfind_commons
  use ramses_commons, only: ramses_t
  use cache_commons, only: msg_maxmass_clump, msg_prop_clump
  use cache
  implicit none
  type(ramses_t)::s
  !----------------------------------------------------------------------------
  ! This subroutine find the 3 most massive clumps in each halo if any.
  !----------------------------------------------------------------------------
  type(msg_maxmass_clump)::dummy_maxmass_clump
  type(msg_prop_clump)::dummy_prop_clump
  integer(kind=8)::global_peak_id,merge_to
  integer::ipeak,jpeak,jpeak1,jpeak2,jpeak3
  real(dp)::mass1,mass2,mass3

  associate(g=>s%g,r=>s%r,m=>s%m,c=>s%c)

  if(g%myid==1.and.r%verbose)write(*,*)'Entering central_halos'

  !--------------------------------------
  ! Find most massive peaks in each halo
  !--------------------------------------

  ! Identify most massive peak within each halo
  c%ind_halo_1=0
  c%ind_max_mass=0
  c%max_peak_mass=0
  call open_cache_clump(s,pack_size=storage_size(dummy_maxmass_clump)/32,&
       init=init_flush_maxmass,flush=pack_flush_maxmass,combine=unpack_flush_maxmass)
  do ipeak=1,c%npeak
     global_peak_id=ipeak+c%npeak_cum(g%myid-1)
     merge_to=c%ind_halo(ipeak)
     call get_peak(s,merge_to,jpeak,flush_cache=.true.,fetch_cache=.false.)
     if(c%clump_mass(ipeak).gt.c%max_peak_mass(jpeak))then
        c%max_peak_mass(jpeak)=c%clump_mass(ipeak)
        c%ind_max_mass(jpeak)=global_peak_id
     endif
  end do
  call close_cache(s,m%grid_dict)
  c%ind_halo_1=c%ind_max_mass

  ! Identify second massive peak within each halo
  c%ind_halo_2=0
  c%ind_max_mass=0
  c%max_peak_mass=0
  call open_cache_clump(s,pack_size=storage_size(dummy_maxmass_clump)/32,&
       pack=pack_fetch_maxmass,unpack=unpack_fetch_maxmass,&
       init=init_flush_maxmass,flush=pack_flush_maxmass,combine=unpack_flush_maxmass)
  do ipeak=1,c%npeak
     global_peak_id=ipeak+c%npeak_cum(g%myid-1)
     merge_to=c%ind_halo(ipeak)
     call get_peak(s,merge_to,jpeak,flush_cache=.true.,fetch_cache=.true.)
     if(c%ind_halo_1(jpeak).NE.global_peak_id)then
        if(c%clump_mass(ipeak).gt.c%max_peak_mass(jpeak))then
           c%max_peak_mass(jpeak)=c%clump_mass(ipeak)
           c%ind_max_mass(jpeak)=global_peak_id
        endif
     endif
  end do
  call close_cache(s,m%grid_dict)
  c%ind_halo_2=c%ind_max_mass

  ! Identify third massive peak within each halo
  c%ind_halo_3=0
  c%ind_max_mass=0
  c%max_peak_mass=0
  call open_cache_clump(s,pack_size=storage_size(dummy_maxmass_clump)/32,&
       pack=pack_fetch_maxmass,unpack=unpack_fetch_maxmass,&
       init=init_flush_maxmass,flush=pack_flush_maxmass,combine=unpack_flush_maxmass)
  do ipeak=1,c%npeak
     global_peak_id=ipeak+c%npeak_cum(g%myid-1)
     merge_to=c%ind_halo(ipeak)
     call get_peak(s,merge_to,jpeak,flush_cache=.true.,fetch_cache=.true.)
     if(c%ind_halo_1(jpeak).NE.global_peak_id.and.c%ind_halo_2(jpeak).NE.global_peak_id)then
        if(c%clump_mass(ipeak).gt.c%max_peak_mass(jpeak))then
           c%max_peak_mass(jpeak)=c%clump_mass(ipeak)
           c%ind_max_mass(jpeak)=global_peak_id
        endif
     endif
  end do
  call close_cache(s,m%grid_dict)
  c%ind_halo_3=c%ind_max_mass

  ! Write masses of relevant centrals if any
  call open_cache_clump(s,pack_size=storage_size(dummy_prop_clump)/32,&
       pack=pack_fetch_central,unpack=unpack_fetch_central)
  do ipeak=1,c%npeak
     global_peak_id=ipeak+c%npeak_cum(g%myid-1)
     if(c%ind_halo(ipeak)==global_peak_id.AND.&
          & c%halo_mass(ipeak) > c%mass_threshold.AND. &
          & c%relevance(ipeak) > c%relevance_threshold)then

        ! Get 3 most massive peak patches and lock them in cache memory
        mass1=0
        mass2=0
        mass3=0
        if(c%ind_halo_1(ipeak).NE.0)then
           global_peak_id=c%ind_halo_1(ipeak)
           call get_peak(s,global_peak_id,jpeak1,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
           if(c%clump_mass(jpeak1) > c%mass_threshold.AND. &
                & c%relevance(jpeak1) > c%relevance_threshold)then
              mass1 = c%clump_mass(jpeak1)
           endif
           if(c%ind_halo_2(ipeak).NE.0)then
              global_peak_id=c%ind_halo_2(ipeak)
              call get_peak(s,global_peak_id,jpeak2,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
              if(c%clump_mass(jpeak2) > c%mass_threshold.AND. &
                   & c%relevance(jpeak2) > c%relevance_threshold)then
                 mass2 = c%clump_mass(jpeak2)
              endif
              if(c%ind_halo_3(ipeak).NE.0)then
                 global_peak_id=c%ind_halo_3(ipeak)
                 call get_peak(s,global_peak_id,jpeak3,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
                 if(c%clump_mass(jpeak3) > c%mass_threshold.AND. &
                      & c%relevance(jpeak3) > c%relevance_threshold)then
                    mass3 = c%clump_mass(jpeak3)
                 endif
              endif
           endif
        endif

        ! Check that we have a most massive clump
        if(mass1.EQ.0)then
           write(*,*)'Problem: we need at least one central galaxy.'
        endif
        
        ! Remove centrals if not massive enough according to several criteria
        if(mass3.LT.0.1*mass1)then
           mass3=0
        endif
        if(mass2.LT.0.1*mass1)then
           mass2=0
        endif

        ! Unlock 3 most mossive peak patches
        if(c%ind_halo_1(ipeak).NE.0)then
           call unlock_cache_clump(s,jpeak1)
           if(c%ind_halo_2(ipeak).NE.0)then
              call unlock_cache_clump(s,jpeak2)
              if(c%ind_halo_3(ipeak).NE.0)then
                 call unlock_cache_clump(s,jpeak3)
              endif
           endif
        endif

        ! Set index of removed central to 0
        if(mass3==0)then
           c%ind_halo_3(ipeak)=0
        endif
        if(mass2==0)then
           c%ind_halo_2(ipeak)=0
        endif
        
!!$        if(mass2>0.or.mass3>0)then
!!$           write(*,*)c%ind_halo(ipeak),c%ind_halo_1(ipeak),c%ind_halo_2(ipeak),c%ind_halo_3(ipeak),&
!!$                & c%halo_mass(ipeak),mass1,mass2,mass3
!!$        endif
!!$
     endif
  end do
  call close_cache(s,m%grid_dict)

  end associate

end subroutine central_in_halos
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_maxmass(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_maxmass_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_maxmass_clump)::msg

  msg%ind1=c%ind_halo_1(local_peak_id)
  msg%ind2=c%ind_halo_2(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_maxmass
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_maxmass(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_maxmass_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_maxmass_clump)::msg

  msg=transfer(msg_array,msg)

  c%ind_halo_1(local_peak_id)=msg%ind1
  c%ind_halo_2(local_peak_id)=msg%ind2

end subroutine unpack_fetch_maxmass
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_maxmass(c,local_peak_id)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  c%max_peak_mass(local_peak_id)=0d0
  c%ind_max_mass(local_peak_id)=0

end subroutine init_flush_maxmass
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_maxmass(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_maxmass_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_maxmass_clump)::msg

  msg%mass=c%max_peak_mass(local_peak_id)
  msg%ind=c%ind_max_mass(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_maxmass
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_maxmass(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_maxmass_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_maxmass_clump)::msg

  msg=transfer(msg_array,msg)

  if(msg%mass.GT.c%max_peak_mass(local_peak_id))then
     c%max_peak_mass(local_peak_id)=msg%mass
     c%ind_max_mass(local_peak_id)=msg%ind
  endif

end subroutine unpack_flush_maxmass
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_central(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg%mass=c%clump_mass(local_peak_id)
  msg%dens=c%relevance(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_central
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_central(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg=transfer(msg_array,msg)

  c%clump_mass(local_peak_id)=msg%mass
  c%relevance(local_peak_id)=msg%dens

end subroutine unpack_fetch_central
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine particle_clump_properties(s,p)
  use amr_parameters, only: ndim,nbin,twotondim,dp
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  !==================================================================
  ! This routine computes various clump properties.
  ! In particular, it computes for each particle its parent peak id.
  ! This is used to compute mass profiles for each halo.
  ! This is also stored in the peak_part and peak_star files.
  ! Written by Romain Teyssier (mini-ramses version in June 2024).
  !==================================================================
  ! Local variables
  type(msg_prop_clump)::dummy_prop_clump
  integer::i,ipeak,ipart,icell,ind,idim,ibin,ilevel
  integer(kind=8)::global_peak_id,global_halo_id
  integer::halo_nr,peak_nr,no_halo
  real(dp)::dist,xx,rad,dx_loc,r2
  real(dp),dimension(1:ndim)::xpart

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c)

  !---------------------------------
  ! Reads peak id of input particles
  !---------------------------------
  call particle_peak_id(s,p,no_halo)

  !--------------------------------------------
  ! Sort particles according to global clump id
  !--------------------------------------------
  call quick_sort_int_int(p%workp(1),p%sortp(1),p%npart)

  !-----------------------------------------
  ! Compute peak velocity based on particles
  !-----------------------------------------
  c%particle_mass=0
  c%peak_vel=0
  call open_cache_clump(s,storage_size(dummy_prop_clump)/32,&
       pack=pack_fetch_part,unpack=unpack_fetch_part,&
       init=init_flush_part,flush=pack_flush_part,combine=unpack_flush_part)
  do i=1+no_halo,p%npart
     ! Get peak id
     ipart=p%sortp(i)
     global_peak_id=p%workp(i)
     if (global_peak_id /=0 ) then
        call get_peak(s,global_peak_id,peak_nr,flush_cache=.true.,fetch_cache=.true.)
        xpart(1:ndim)=p%xp(ipart,1:ndim)-c%peak_pos(peak_nr,1:ndim)
        ! In case of periodic boundaries
        r2=0d0
        do idim=1,ndim
           if(xpart(idim)> r%boxlen*0.5)xpart(idim)=xpart(idim)-r%boxlen
           if(xpart(idim)<-r%boxlen*0.5)xpart(idim)=xpart(idim)+r%boxlen
           r2=r2+xpart(idim)**2
        end do
        ilevel=c%peak_level(peak_nr)
        dx_loc=r%boxlen/2**ilevel
        ! Keep only particles within a 2-cell radius
        if(r2<=4d0*dx_loc**2)then
           c%peak_vel(peak_nr,1:ndim)=c%peak_vel(peak_nr,1:ndim)+p%mp(ipart)*p%vp(ipart,1:ndim)
           c%particle_mass(peak_nr)=c%particle_mass(peak_nr)+p%mp(ipart)
        endif
     endif
  end do
  call close_cache(s,m%grid_dict)

  ! Compute specific quantities
  do ipeak=1,c%npeak
     if (c%particle_mass(ipeak)>0)then
        c%peak_vel(ipeak,1:ndim)=c%peak_vel(ipeak,1:ndim)/c%particle_mass(ipeak)
     end if
  end do

  end associate

end subroutine particle_clump_properties
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_part(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg%pos(1:ndim)=c%peak_pos(local_peak_id,1:ndim)
  msg%vol=c%peak_level(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_part(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg=transfer(msg_array,msg)

  c%peak_pos(local_peak_id,1:ndim)=msg%pos(1:ndim)
  c%peak_level(local_peak_id)=msg%vol

end subroutine unpack_fetch_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_part(c,local_peak_id)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  c%particle_mass(local_peak_id)=0
  c%peak_vel(local_peak_id,1:ndim)=0d0

end subroutine init_flush_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_part(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg%mass=c%particle_mass(local_peak_id)
  msg%vel(1:ndim)=c%peak_vel(local_peak_id,1:ndim)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_part(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg=transfer(msg_array,msg)

  c%particle_mass(local_peak_id)=c%particle_mass(local_peak_id)+msg%mass
  c%peak_vel(local_peak_id,1:ndim)=c%peak_vel(local_peak_id,1:ndim)+msg%vel(1:ndim)

end subroutine unpack_flush_part
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine particle_split_centrals(s,p)
  use amr_parameters, only: ndim,nbin,twotondim,dp
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  !==================================================================
  ! This routine computes various clump properties.
  ! In particular, it computes for each particle its parent peak id.
  ! This is used to compute mass profiles for each halo.
  ! This is also stored in the peak_part and peak_star files.
  ! Written by Romain Teyssier (mini-ramses version in June 2024).
  !==================================================================
  ! Local variables
  type(msg_prop_clump)::dummy_prop_clump
  type(msg_mbin_clump)::dummy_mbin_clump
  integer::i,ipart,icell,ind,idim,ibin,ilevel
  integer(kind=8)::global_peak_id,global_halo_id
  integer::ipeak,jpeak,no_halo
  real(dp)::pi,grav,radius,velocity,distmin
  real(dp)::dist,dist1,dist2,dist3
  real(dp)::xdist1,xdist2,xdist3
  real(dp)::vdist1,vdist2,vdist3
  real(dp),dimension(1:ndim)::xpart,vpart

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c)

  ! Constants
  pi=ACOS(-1.0D0)
  grav=1d0
  if(s%r%cosmo)grav=3d0/8d0/pi*s%g%omega_m*s%g%aexp

  !---------------------------------
  ! Reads halo id of input particles
  !---------------------------------
  call particle_halo_id(s,p,no_halo)

  !--------------------------------------------
  ! Sort particles according to global clump id
  !--------------------------------------------
  call quick_sort_int_int(p%workp(1),p%sortp(1),p%npart)

  !-----------------------------------------------------------
  ! Assign particle to central using clustering in phase space
  !-----------------------------------------------------------
  call open_cache_clump(s,storage_size(dummy_prop_clump)/32,&
       pack=pack_fetch_split,unpack=unpack_fetch_split)
  do i=1+no_halo,p%npart
     ! Get halo id
     ipart=p%sortp(i)
     global_halo_id=p%workp(i)
     if (global_halo_id /=0 ) then
        call get_peak(s,global_halo_id,ipeak,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
        ! Compute halo radius
        radius=(c%halo_mass(ipeak)/4d0/pi*3d0/200d0)**(1d0/3d0)
        ! Get first central peak id
        global_peak_id=c%ind_halo_1(ipeak)
        dist1=1e10
        if(global_peak_id>0)then
           call get_peak(s,global_peak_id,jpeak,flush_cache=.false.,fetch_cache=.true.)
           ! Compute Euclidian distance in configuration space
           velocity=sqrt(grav*c%clump_mass(jpeak)/radius)
           dist1=cmp_distance(p%xp(ipart,1:ndim),c%peak_pos(jpeak,1:ndim),radius, &
                &             p%vp(ipart,1:ndim),c%peak_vel(jpeak,1:ndim),velocity,r%boxlen)
        endif
        ! Get second central peak id
        global_peak_id=c%ind_halo_2(ipeak)
        dist2=1e10
        if(global_peak_id>0)then
           call get_peak(s,global_peak_id,jpeak,flush_cache=.false.,fetch_cache=.true.)
           ! Compute Euclidian distance in configuration space
           velocity=sqrt(grav*c%clump_mass(jpeak)/radius)
           dist2=cmp_distance(p%xp(ipart,1:ndim),c%peak_pos(jpeak,1:ndim),radius, &
                &             p%vp(ipart,1:ndim),c%peak_vel(jpeak,1:ndim),velocity,r%boxlen)
        endif
        ! Get third central peak id
        global_peak_id=c%ind_halo_3(ipeak)
        dist3=1e10
        if(global_peak_id>0)then
           call get_peak(s,global_peak_id,jpeak,flush_cache=.false.,fetch_cache=.true.)
           ! Compute Euclidian distance in configuration space
           velocity=sqrt(grav*c%clump_mass(jpeak)/radius)
           dist3=cmp_distance(p%xp(ipart,1:ndim),c%peak_pos(jpeak,1:ndim),radius, &
                &             p%vp(ipart,1:ndim),c%peak_vel(jpeak,1:ndim),velocity,r%boxlen)
        endif
        ! Assign particle to closest central in phase space
        distmin=min(dist1,min(dist2,dist3))
        if(dist3.EQ.distmin)p%workp(i)=c%ind_halo_3(ipeak)
        if(dist2.EQ.distmin)p%workp(i)=c%ind_halo_2(ipeak)
        if(dist1.EQ.distmin)p%workp(i)=c%ind_halo_1(ipeak)
        ! Unlock halo
        call unlock_cache_clump(s,ipeak)
     endif
  end do
  call close_cache(s,m%grid_dict)

  !-------------------------------------------------------
  ! Assign peak to central using clustering in phase space
  !-------------------------------------------------------
  call open_cache_clump(s,storage_size(dummy_prop_clump)/32, &
       pack=pack_fetch_split,unpack=unpack_fetch_split)
  do i=1,c%npeak
     if(c%clump_mass(i) > c%mass_threshold.AND. &
          & c%relevance(i) > c%relevance_threshold)then
     global_halo_id=c%ind_halo(i)
     if(global_halo_id>0)then
        call get_peak(s,global_halo_id,ipeak,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
        ! Compute halo radius
        radius=(c%halo_mass(ipeak)/4d0/pi*3d0/200d0)**(1d0/3d0)
        ! Get first central peak id
        global_peak_id=c%ind_halo_1(ipeak)
        dist1=1e10
        if(global_peak_id>0)then
           call get_peak(s,global_peak_id,jpeak,flush_cache=.false.,fetch_cache=.true.)
           ! Compute Euclidian distance in configuration space
           velocity=sqrt(grav*c%clump_mass(jpeak)/radius)
           dist1=cmp_distance(c%peak_pos(i,1:ndim),c%peak_pos(jpeak,1:ndim),radius, &
                &             c%peak_vel(i,1:ndim),c%peak_vel(jpeak,1:ndim),velocity,r%boxlen)
        endif
        ! Get second central peak id
        global_peak_id=c%ind_halo_2(ipeak)
        dist2=1e10
        if(global_peak_id>0)then
           call get_peak(s,global_peak_id,jpeak,flush_cache=.false.,fetch_cache=.true.)
           ! Compute Euclidian distance in configuration space
           velocity=sqrt(grav*c%clump_mass(jpeak)/radius)
           dist2=cmp_distance(c%peak_pos(i,1:ndim),c%peak_pos(jpeak,1:ndim),radius, &
                &             c%peak_vel(i,1:ndim),c%peak_vel(jpeak,1:ndim),velocity,r%boxlen)
        endif
        ! Get third central peak id
        global_peak_id=c%ind_halo_3(ipeak)
        dist3=1e10
        if(global_peak_id>0)then
           call get_peak(s,global_peak_id,jpeak,flush_cache=.false.,fetch_cache=.true.)
           ! Compute Euclidian distance in configuration space
           velocity=sqrt(grav*c%clump_mass(jpeak)/radius)
           dist3=cmp_distance(c%peak_pos(i,1:ndim),c%peak_pos(jpeak,1:ndim),radius, &
                &             c%peak_vel(i,1:ndim),c%peak_vel(jpeak,1:ndim),velocity,r%boxlen)
        endif
        ! Assign particle to closest central in phase space
        distmin=min(dist1,min(dist2,dist3))
        if(dist3.EQ.distmin)c%ind_central(i)=c%ind_halo_3(ipeak)
        if(dist2.EQ.distmin)c%ind_central(i)=c%ind_halo_2(ipeak)
        if(dist1.EQ.distmin)c%ind_central(i)=c%ind_halo_1(ipeak)
        ! Unlock halo
        call unlock_cache_clump(s,ipeak)
     endif
     endif
  end do
  call close_cache(s,m%grid_dict)

  !--------------------------------------------------------
  ! Compute particle mass profile around their central peak
  !--------------------------------------------------------
  call open_cache_clump(s,storage_size(dummy_mbin_clump)/32,&
       pack=pack_fetch_mbin,unpack=unpack_fetch_mbin,&
       init=init_flush_mbin,flush=pack_flush_mbin,combine=unpack_flush_mbin)
  c%mass_bin=0d0
  do i=1+no_halo,p%npart
     ! Get central global peak id
     ipart=p%sortp(i)
     global_peak_id=p%workp(i)
     if (global_peak_id /=0 ) then
        call get_peak(s,global_peak_id,ipeak,flush_cache=.true.,fetch_cache=.true.)
        ! Compute halo maximum radius
        radius=2d0*(c%halo_mass(ipeak)/4d0/pi*3d0/200d0)**(1d0/3d0)
        ! Compute particle radius
        dist=0d0
        xpart(1:ndim)=p%xp(ipart,1:ndim)-c%peak_pos(ipeak,1:ndim)
        ! In case of periodic boundaries
        do idim=1,ndim
           if(xpart(idim)> r%boxlen*0.5)xpart(idim)=xpart(idim)-r%boxlen
           if(xpart(idim)<-r%boxlen*0.5)xpart(idim)=xpart(idim)+r%boxlen
           dist=dist+xpart(idim)**2
        end do
        dist=sqrt(dist)
        do ibin=1,nbin
           ! We use a simple linear binning as the mass is usually propto r
           if(dist<=dble(ibin)/dble(nbin)*radius)then
              c%mass_bin(ipeak,ibin)=c%mass_bin(ipeak,ibin)+p%mp(ipart)
              exit
           endif
        end do
     endif
  end do
  call close_cache(s,m%grid_dict)

  !------------------------
  ! Compute cumulative mass
  !------------------------
  do ipeak=1,c%npeak
     if(c%ind_central(ipeak).EQ.ipeak+c%npeak_cum(g%myid-1).AND. &
          & c%halo_mass(ipeak) > c%mass_threshold.AND. &
          & c%relevance(ipeak) > c%relevance_threshold)then
        do ibin=1,nbin-1
           c%mass_bin(ipeak,ibin+1)=c%mass_bin(ipeak,ibin+1)+c%mass_bin(ipeak,ibin)
        end do
     endif
  end do

  end associate

end subroutine particle_split_centrals
!################################################################
!################################################################
!################################################################
!################################################################
function cmp_distance(x1,x2,radius,v1,v2,velocity,boxlen)
  use amr_parameters, only: ndim, dp
  real(dp),dimension(1:ndim)::x1,x2,v1,v2
  real(dp)::radius,velocity,boxlen
  real(dp)::cmp_distance

  integer::idim
  real(dp)::xdist,vdist
  real(dp),dimension(1:ndim)::xpart,vpart

  xdist=0d0
  xpart(1:ndim)=x1(1:ndim)-x2(1:ndim)
  ! In case of periodic boundaries
  do idim=1,ndim
     if(xpart(idim)> boxlen*0.5)xpart(idim)=xpart(idim)-boxlen
     if(xpart(idim)<-boxlen*0.5)xpart(idim)=xpart(idim)+boxlen
     xdist=xdist+xpart(idim)**2
  end do
  ! Rescale distance in configuration space
  xdist=xdist/radius**2
  ! Compute Euclidian distance in velocity space
  vdist=0d0
  vpart(1:ndim)=v1(1:ndim)-v2(1:ndim)
  do idim=1,ndim
     vdist=vdist+vpart(idim)**2
  end do
  ! Rescale distance in velocity space
  vdist=vdist/velocity**2
  cmp_distance=xdist+vdist
end function cmp_distance
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_split(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg%mass=c%clump_mass(local_peak_id)
  msg%dens=c%halo_mass(local_peak_id)
  msg%vel=c%peak_vel(local_peak_id,1:ndim)
  msg%pos=c%peak_pos(local_peak_id,1:ndim)
  msg%ind(1)=c%ind_halo_1(local_peak_id)
  msg%ind(2)=c%ind_halo_2(local_peak_id)
  msg%ind(3)=c%ind_halo_3(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_split
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_split(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg=transfer(msg_array,msg)

  c%clump_mass(local_peak_id)=msg%mass
  c%halo_mass(local_peak_id)=msg%dens
  c%peak_vel(local_peak_id,1:ndim)=msg%vel
  c%peak_pos(local_peak_id,1:ndim)=msg%pos
  c%ind_halo_1(local_peak_id)=msg%ind(1)
  c%ind_halo_2(local_peak_id)=msg%ind(2)
  c%ind_halo_3(local_peak_id)=msg%ind(3)

end subroutine unpack_fetch_split
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_mbin(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_mbin_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_mbin_clump)::msg

  msg%pos(1:ndim)=c%peak_pos(local_peak_id,1:ndim)
  msg%mass=c%halo_mass(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_mbin
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_mbin(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_mbin_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_mbin_clump)::msg

  msg=transfer(msg_array,msg)

  c%peak_pos(local_peak_id,1:ndim)=msg%pos(1:ndim)
  c%halo_mass(local_peak_id)=msg%mass

end subroutine unpack_fetch_mbin
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_mbin(c,local_peak_id)
  use amr_commons, only: nbin
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  c%mass_bin(local_peak_id,1:nbin)=0d0

end subroutine init_flush_mbin
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_mbin(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: nbin
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_mbin_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_mbin_clump)::msg

  msg%mbin(1:nbin)=c%mass_bin(local_peak_id,1:nbin)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_mbin
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_mbin(c,local_peak_id,msg_size,msg_array)
  use amr_commons, only: nbin
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_mbin_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_mbin_clump)::msg

  msg=transfer(msg_array,msg)

  c%mass_bin(local_peak_id,1:nbin)=c%mass_bin(local_peak_id,1:nbin)+msg%mbin(1:nbin)

end subroutine unpack_flush_mbin
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine particle_peak_id(s,p,no_peak)
  use amr_parameters, only: ndim,nbin,twotondim,dp
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use nbors_utils
  use cache_commons
  use cache
  use boundaries, only: init_bound_flag2
  use marshal, only: pack_fetch_flag2, unpack_fetch_flag2
  use hilbert
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::no_peak
  !==================================================================
  ! This routine reads from the grid peak map (flag2) the peak id
  ! of the input particle object. It could be dark matter or stars.
  ! Written by Romain Teyssier (mini-ramses version in June 2024).
  !==================================================================
  ! Local variables
  integer,dimension(1:ndim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_cell
  integer::i,ipart,icell,ind,idim,ibin,ilevel
  integer(kind=8)::global_peak_id
  integer::local_peak_id,ipeak,jpeak,merge_to
  integer::halo_nr,peak_nr
  real(dp)::dx_loc,rmin,rmax,dist,xx
  type(oct),pointer::gridp
  type(msg_int4)::dummy_int4
  logical::ok_level,ok_leaf

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c)

  ! Open cache for array flag1 (fetch)
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
       hilbert=m%domain,pack_size=storage_size(dummy_int4)/32,&
       pack=pack_fetch_flag2,unpack=unpack_fetch_flag2,bound=init_bound_flag2)

  ! Loop over particles
  no_peak=0
  do ilevel=r%levelmin,r%nlevelmax

     ! Mesh spacing in that level
     dx_loc=r%boxlen/2**ilevel

     do ipart=p%headp(ilevel),p%tailp(ilevel)

        ok_level=.true.

        ! Find parent cell at level ilevel
        do idim=1,ndim
           ckey(idim)=int(p%xp(ipart,idim)/dx_loc)
        end do

        ! Get parent cell at level ilevel using cache
        hash_cell(0)=ilevel+1
        hash_cell(1:ndim)=ckey(1:ndim)
        call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)

        ! If cell does not exist at current level, then find cell at coarser level
        if(.not.associated(gridp))then

           ! NGP at level ilevel-1
           do idim=1,ndim
              ckey(idim)=int(p%xp(ipart,idim)/dx_loc/2)
           end do

           ! Get parent cell at level ilevel-1 using cache
           hash_cell(0)=ilevel
           hash_cell(1:ndim)=ckey(1:ndim)
           call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
           if(.not.associated(gridp))ok_level=.false.

        end if

        if(.not. ok_level)then
           write(*,*)"Something went wrong in particle_peak_id"
           write(*,*)"Current level grid and coarser grid both dont exist..."
           stop
        endif

        ! Read flag2 value
        global_peak_id=gridp%flag2(icell)
        if (global_peak_id==0)no_peak=no_peak+1

        ! Store global peak id in workp array
        p%sortp(ipart)=ipart
        p%workp(ipart)=global_peak_id

     end do
     ! End loop over particles
  end do
  ! End loop over levels

  call close_cache(s,m%grid_dict)

  end associate

end subroutine particle_peak_id
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine particle_halo_id(s,p,no_peak)
  use amr_parameters, only: ndim,nbin,twotondim,dp
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use nbors_utils
  use cache_commons
  use cache
  use boundaries, only: init_bound_flag
  use marshal, only: pack_fetch_flag, unpack_fetch_flag
  use hilbert
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::no_peak
  !==================================================================
  ! This routine reads from the grid peak map (flag1) the peak id
  ! of the input particle object. It could be dark matter or stars.
  ! Written by Romain Teyssier (mini-ramses version in June 2024).
  !==================================================================
  ! Local variables
  integer,dimension(1:ndim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_cell
  integer::i,ipart,icell,ind,idim,ibin,ilevel
  integer(kind=8)::global_peak_id
  integer::local_peak_id,ipeak,jpeak,merge_to
  integer::halo_nr,peak_nr
  real(dp)::dx_loc,rmin,rmax,dist,xx
  type(oct),pointer::gridp
  type(msg_int4)::dummy_int4
  logical::ok_level,ok_leaf

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c)

  ! Open cache for array flag1 (fetch)
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
       hilbert=m%domain,pack_size=storage_size(dummy_int4)/32,&
       pack=pack_fetch_flag,unpack=unpack_fetch_flag,bound=init_bound_flag)

  ! Loop over particles
  no_peak=0
  do ilevel=r%levelmin,r%nlevelmax

     ! Mesh spacing in that level
     dx_loc=r%boxlen/2**ilevel

     do ipart=p%headp(ilevel),p%tailp(ilevel)

        ok_level=.true.

        ! Find parent cell at level ilevel
        do idim=1,ndim
           ckey(idim)=int(p%xp(ipart,idim)/dx_loc)
        end do

        ! Get parent cell at level ilevel using cache
        hash_cell(0)=ilevel+1
        hash_cell(1:ndim)=ckey(1:ndim)
        call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)

        ! If cell does not exist at current level, then find cell at coarser level
        if(.not.associated(gridp))then

           ! NGP at level ilevel-1
           do idim=1,ndim
              ckey(idim)=int(p%xp(ipart,idim)/dx_loc/2)
           end do

           ! Get parent cell at level ilevel-1 using cache
           hash_cell(0)=ilevel
           hash_cell(1:ndim)=ckey(1:ndim)
           call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
           if(.not.associated(gridp))ok_level=.false.

        end if

        if(.not. ok_level)then
           write(*,*)"Something went wrong in particle_halo_id"
           write(*,*)"Current level grid and coarser grid both dont exist..."
           stop
        endif

        ! Read flag1 value
        global_peak_id=gridp%flag1(icell)
        if (global_peak_id==0)no_peak=no_peak+1

        ! Store global halo id in workp array
        p%sortp(ipart)=ipart
        p%workp(ipart)=global_peak_id

     end do
     ! End loop over particles
  end do
  ! End loop over levels

  call close_cache(s,m%grid_dict)

  end associate

end subroutine particle_halo_id
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module clump_merger_module
