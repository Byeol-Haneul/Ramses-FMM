module clump_merger_module
contains
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
  use sparse_matrix
  implicit none
  type(ramses_t)::s

  associate(g=>s%g,m=>s%m,c=>s%c)

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

  ! Deallocate peak patch arrays
  deallocate(c%peak_cell)
  deallocate(c%peak_grid)
  deallocate(c%peak_level)
  deallocate(c%max_dens)

  deallocate(c%n_cells)
  deallocate(c%n_cells_halo)
  deallocate(c%lev_peak)
  deallocate(c%new_peak)
  deallocate(c%ind_halo)
  deallocate(c%halo_mass)
  deallocate(c%clump_mass)
  deallocate(c%relevance)

  deallocate(c%clump_size)
  deallocate(c%peak_pos)
  deallocate(c%center_of_mass)
  deallocate(c%min_dens)
  deallocate(c%av_dens)
  deallocate(c%clump_vol)

  ! Deallocate sparse density matrix
  call sparse_kill(c%sparse_saddle_dens)

  ! Deallocate hash table
  deallocate(c%hkey,c%gkey,c%nkey)

  end associate

end subroutine deallocate_peak_patch_arrays
!################################################################
!################################################################
!################################################################
!################################################################
subroutine allocate_peak_patch_arrays(s)
  use amr_parameters, ONLY: ndim, dp
  use clfind_commons
  use ramses_commons, ONLY: ramses_t
  use sparse_matrix
  implicit none
  type(ramses_t)::s
  
  integer::bit_length,ncode,igrid,ind
  integer::itest,global_peak_id,local_peak_id

  associate(g=>s%g,m=>s%m,c=>s%c)

  !-------------------------------
  ! Allocate peak-patch_properties
  !-------------------------------
  allocate(c%n_cells(1:c%npeak_max))
  allocate(c%n_cells_halo(1:c%npeak_max))
  allocate(c%lev_peak(1:c%npeak_max))
  allocate(c%new_peak(c%npeak_max))
  allocate(c%ind_halo(1:c%npeak_max))
  allocate(c%halo_mass(1:c%npeak_max))
  allocate(c%clump_mass(1:c%npeak_max))
  allocate(c%relevance(1:c%npeak_max))

  allocate(c%clump_size(1:c%npeak_max,1:ndim))
  allocate(c%peak_pos(1:c%npeak_max,1:ndim))
  allocate(c%center_of_mass(1:c%npeak_max,1:ndim))
  allocate(c%min_dens(1:c%npeak_max))
  allocate(c%av_dens(1:c%npeak_max))
  allocate(c%clump_vol(1:c%npeak_max))

  !-------------------------------------------
  ! Initialize sparse matrix for saddle points
  !-------------------------------------------
  call sparse_initialize(c%npeak_max,c%sparse_saddle_dens)

  !--------------------
  ! Allocate hash table
  !--------------------
  ncode=c%npeak_max-c%npeak
  do bit_length=1,32
     ncode=ncode/2
     if(ncode<=1) exit
  end do
  c%nhash=c%prime(bit_length+1)
  allocate(c%hkey(1:c%nhash))
  c%hfree=c%npeak+1
  c%hcollision=0
  allocate(c%gkey(c%npeak+1:c%npeak_max))
  allocate(c%nkey(c%npeak+1:c%npeak_max))
  c%hkey=0; c%gkey=0; c%nkey=0

  !------------------------------------------------
  ! Initialize the hash table with interior patches
  !------------------------------------------------
  do itest=1,c%ntest
     igrid=c%grid(itest)
     ind=c%cell(itest)
     global_peak_id=m%grid(igrid)%flag1(ind) ! global peak id
     if (global_peak_id>0)then
        call get_local_peak_id(s,global_peak_id,local_peak_id)
     end if
  end do

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
subroutine get_local_peak_id(s,global_peak_id,local_peak_id)
  use amr_commons
  use ramses_commons, only: ramses_t
  use clfind_commons
  implicit none
  integer::global_peak_id,local_peak_id
  type(ramses_t)::s
  
  integer::ihash,ikey,jkey

  associate(g=>s%g,c=>s%c)

  if(    global_peak_id > c%npeak_cum(g%myid-1) .and. &
       & global_peak_id <= c%npeak_cum(g%myid))then
     local_peak_id=global_peak_id-c%npeak_cum(g%myid-1)
  else
     ihash=MOD(global_peak_id,c%nhash)+1 ! compute the simple prime hash key
     if(c%hkey(ihash)==0)then ! hash table is empty
        c%hkey(ihash)=c%hfree
        local_peak_id=c%hfree
        c%gkey(c%hfree)=global_peak_id
        c%hfree=c%hfree+1
        if(c%hfree.eq.c%npeak_max)then
           write(*,*)'Too many peaks'
           write(*,*)'Increase npeak_max'
           stop
        endif
     else
        ikey=c%hkey(ihash) ! collision in the hash table
        do while(ikey>0)
           jkey=ikey
           if(c%gkey(ikey)==global_peak_id)exit
           ikey=c%nkey(ikey)
        end do
        if(ikey==0)then ! peak doesn't already exist
           c%nkey(jkey)=c%hfree
           local_peak_id=c%hfree
           c%gkey(c%hfree)=global_peak_id
           c%hfree=c%hfree+1
           c%hcollision=c%hcollision+1
           if(c%hfree.eq.c%npeak_max)then
              write(*,*)'Too many peaks'
              write(*,*)'Increase npeak_max'
              stop
           endif
        else ! peak already exists
           local_peak_id=ikey
        end if
     end if
  end if

  end associate

end subroutine get_local_peak_id
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_local_peak_cpu(s,local_peak_id,peak_cpu)
  use amr_commons
  use clfind_commons
  use ramses_commons, only: ramses_t
  implicit none
  integer::local_peak_id,peak_cpu
  type(ramses_t)::s
  ! get the mpi-domain a peak belongs to from its LOCAL id  
  integer::icpu,global_peak_id

  associate(g=>s%g,c=>s%c)

  if(local_peak_id <= c%npeak) then
     peak_cpu = g%myid
  else
     global_peak_id = c%gkey(local_peak_id)
     peak_cpu = g%ncpu
     do icpu = 1,g%ncpu
        if(    global_peak_id > c%npeak_cum(icpu-1) .and. &
             & global_peak_id <= c%npeak_cum(icpu))then
           peak_cpu = icpu
        endif
     end do
  end if

  end associate

end subroutine get_local_peak_cpu
!################################################################
!################################################################
!################################################################
!################################################################
subroutine build_peak_communicator(s)
  use amr_commons
  use ramses_commons, only: ramses_t
  use clfind_commons
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s

#ifndef WITHOUTMPI
  integer::info,ipeak,icpu
  integer,dimension(1:s%g%ncpu,1:s%g%ncpu)::npeak_alltoall
  integer,dimension(1:s%g%ncpu,1:s%g%ncpu)::npeak_alltoall_tot
  integer,dimension(1:s%g%ncpu)::ipeak_alltoall

  associate(g=>s%g,c=>s%c)

  npeak_alltoall=0
  do ipeak=c%npeak+1,c%hfree-1
     call get_local_peak_cpu(s,ipeak,icpu)
     npeak_alltoall(g%myid,icpu)=npeak_alltoall(g%myid,icpu)+1
  end do
  call MPI_ALLREDUCE(npeak_alltoall,npeak_alltoall_tot,g%ncpu*g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  npeak_alltoall=npeak_alltoall_tot
  if(.not. allocated(c%peak_send_cnt))then
     allocate(c%peak_send_cnt(1:g%ncpu),c%peak_send_oft(1:g%ncpu))
     allocate(c%peak_recv_cnt(1:g%ncpu),c%peak_recv_oft(1:g%ncpu))
  endif
  c%peak_send_cnt=0; c%peak_send_oft=0; c%peak_send_tot=0
  c%peak_recv_cnt=0; c%peak_recv_oft=0; c%peak_recv_tot=0
  do icpu=1,g%ncpu
     c%peak_send_cnt(icpu)=npeak_alltoall(g%myid,icpu)
     c%peak_recv_cnt(icpu)=npeak_alltoall(icpu,g%myid)
     c%peak_send_tot=c%peak_send_tot+c%peak_send_cnt(icpu)
     c%peak_recv_tot=c%peak_recv_tot+c%peak_recv_cnt(icpu)
     if(icpu<g%ncpu)then
        c%peak_send_oft(icpu+1)=c%peak_send_oft(icpu)+npeak_alltoall(g%myid,icpu)
        c%peak_recv_oft(icpu+1)=c%peak_recv_oft(icpu)+npeak_alltoall(icpu,g%myid)
     endif
  end do
  if(allocated(c%peak_send_buf))then
     deallocate(c%peak_send_buf,c%peak_recv_buf)
  endif
  allocate(c%peak_send_buf(1:c%peak_send_tot))
  allocate(c%peak_recv_buf(1:c%peak_recv_tot))
  ipeak_alltoall=0
  do ipeak=c%npeak+1,c%hfree-1
     call get_local_peak_cpu(s,ipeak,icpu)
     ipeak_alltoall(icpu)=ipeak_alltoall(icpu)+1
     c%peak_send_buf(c%peak_send_oft(icpu)+ipeak_alltoall(icpu))=c%gkey(ipeak)
  end do
  call MPI_ALLTOALLV(c%peak_send_buf,c%peak_send_cnt,c%peak_send_oft,MPI_INTEGER, &
       &             c%peak_recv_buf,c%peak_recv_cnt,c%peak_recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  
  end associate

#endif
end subroutine build_peak_communicator
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine merge_clumps(s,action)
  use amr_commons, only: dp, ndim
  use ramses_commons, only: ramses_t
  use sparse_matrix
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
  integer::j,i,merge_to,ipart,igrid,ind,itest
  integer::current,nmove,ipeak,jpeak,iter
  integer::nsurvive,nzero,idepth
  integer::ilev,ilevel,global_peak_id,mergelevel_max
  real(dp)::value_iij,zero=0,relevance_peak
  real(dp)::d,dx_loc,vol,mass_min
  integer,dimension(1:s%c%npeak_max)::alive,ind_sort
  real(dp),dimension(1:s%c%npeak_max)::peakd
  logical::do_merge=.false.

#ifndef WITHOUTMPI
  integer::mergelevel_max_global
  integer::nmove_all,nsurvive_all,nzero_all
#endif

  associate(g=>s%g,r=>s%r,m=>s%m,c=>s%c)

  mass_min=g%mp_min*r%mass_threshold

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
        if(c%relevance(i)>r%relevance_threshold)then
           alive(i)=1
        else
           alive(i)=0
        endif
     endif
  end do

  ! Compute clump mass
  c%clump_mass=0
  do itest=1,c%ntest
     ilevel=c%level(itest)
     igrid=c%grid(itest)
     ind=c%cell(itest)
     global_peak_id=m%grid(igrid)%flag1(ind)
     if (global_peak_id>0) then
        call get_local_peak_id(s,global_peak_id,ipeak)
#ifdef GRAV
        d=m%grid(igrid)%rho(ind) ! Cell density
#endif
        dx_loc=r%boxlen/2**ilevel ! Cell size
        vol=dx_loc**ndim  ! Cell volume
        c%clump_mass(ipeak)=c%clump_mass(ipeak)+vol*d
     end if
  end do
#ifndef WITHOUTMPI
  ! Update communicator with possibly new remote peaks imported locally
  call build_peak_communicator(s)
  ! Collect results from all MPI domains
  call virtual_peak_dp(s,c%clump_mass,'sum')
  call boundary_peak_dp(s,c%clump_mass)
#endif

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

     ! Compute maximum saddle density for each clump
     call get_max(s,c%hfree-1,c%sparse_saddle_dens)

#ifndef WITHOUTMPI
     ! Create new local duplicated peaks and update communicator
     call virtual_saddle_max(s)
     call build_peak_communicator(s)

     ! Set up bounday values
     call boundary_peak_dp(s,c%sparse_saddle_dens%maxval)
     call boundary_peak_int(s,c%sparse_saddle_dens%maxloc)
     call boundary_peak_dp(s,c%max_dens)
     call boundary_peak_int(s,c%new_peak)
     call boundary_peak_int(s,alive)
#endif

     ! Merge peaks
     nmove=c%npeak_tot
     iter=0
     do while(nmove>0)
        nmove=0
        do i=c%npeak,1,-1
           ipeak=ind_sort(i)
           merge_to=c%new_peak(ipeak)
           if(alive(ipeak)>0)then
              if(action.EQ.'relevance')then
                 if(c%sparse_saddle_dens%maxval(ipeak)>0)then
                    relevance_peak=c%max_dens(ipeak)/c%sparse_saddle_dens%maxval(ipeak)
                 else
                    relevance_peak=c%max_dens(ipeak)/r%density_threshold
                 end if
                 do_merge=(relevance_peak<r%relevance_threshold.OR.c%clump_mass(ipeak)<mass_min)
              endif
              if(action.EQ.'saddleden')then
                 do_merge=(c%sparse_saddle_dens%maxval(ipeak)>r%saddle_threshold)
              endif
              if(do_merge)then
                 if(c%sparse_saddle_dens%maxloc(ipeak)>0)then
                    call get_local_peak_id(s,c%sparse_saddle_dens%maxloc(ipeak),jpeak)
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
        ! Update boundary conditions for new_peak array
        call boundary_peak_int(s,c%new_peak)
        iter=iter+1
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(nmove,nmove_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        nmove=nmove_all
#endif
        if(g%myid==1.and.r%verbose)write(*,*)'niter=',iter,'nmove=',nmove
     end do

     ! Transfer matrix elements of merged peaks to surviving peaks
     ! Create new local duplicated (remote) peaks and update communicator
     do ipeak=1,c%hfree-1
        if(alive(ipeak)>0)then
           merge_to=c%new_peak(ipeak)
           if(ipeak.LE.c%npeak)then
              global_peak_id=c%npeak_cum(g%myid-1)+ipeak
           else
              global_peak_id=c%gkey(ipeak)
           endif
           if(merge_to.NE.global_peak_id)then
              call get_local_peak_id(s,merge_to,jpeak)
              current=c%sparse_saddle_dens%first(ipeak) ! first element of line ipeak
              do while(current>0) ! walk the line
                 j=c%sparse_saddle_dens%col(current)
                 value_iij=c%sparse_saddle_dens%val(current) ! value of the matrix
                 ! Copy the value of density only if larger
                 if(value_iij>get_value(jpeak,j,c%sparse_saddle_dens))then
                    call set_value(jpeak,j,value_iij,c%sparse_saddle_dens)
                    call set_value(j,jpeak,value_iij,c%sparse_saddle_dens)
                 end if
                 current=c%sparse_saddle_dens%next(current)
              end do
              call set_value(jpeak,jpeak,zero,c%sparse_saddle_dens)
           end if
        endif
     end do
     call build_peak_communicator(s)

     if(action.EQ.'relevance')then

        ! Update merged peak mass using halo mass
        c%halo_mass=0
        do ipeak=1,c%npeak
           if(alive(ipeak)>0)then
              merge_to=c%new_peak(ipeak)
              call get_local_peak_id(s,merge_to,jpeak)
              c%halo_mass(jpeak)=c%halo_mass(jpeak)+c%clump_mass(ipeak)
           endif
        end do
        call build_peak_communicator(s)
        call virtual_peak_dp(s,c%halo_mass,'sum')
        call boundary_peak_dp(s,c%halo_mass)

     endif

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
     call boundary_peak_int(s,alive)

     if(action.EQ.'relevance')then

        ! Assign halo mass to alive peak mass
        do ipeak=1,c%npeak
           if(alive(ipeak)>0)then
              c%clump_mass(ipeak)=c%halo_mass(ipeak)
           endif
        end do

     endif

     ! Remove all matrix elements corresponding to merged peaks
     do ipeak=1,c%hfree-1
        current=c%sparse_saddle_dens%first(ipeak) ! first element of line ipeak
        do while(current>0) ! walk the line
           j=c%sparse_saddle_dens%col(current)
           if(alive(ipeak)==0 .or. alive(j)==0)then
              call set_value(ipeak,j,zero,c%sparse_saddle_dens)
           endif
           current=c%sparse_saddle_dens%next(current)
        end do
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

  ! Compute maximum saddle density for each surviving clump
  ! Create new local duplicated peaks and update communicator
  call get_max(s,c%hfree-1,c%sparse_saddle_dens)

#ifndef WITHOUTMPI
  call virtual_saddle_max(s)
  call build_peak_communicator(s)
#endif

  ! Set up bounday values
  call boundary_peak_dp(s,c%sparse_saddle_dens%maxval)
  call boundary_peak_int(s,c%sparse_saddle_dens%maxloc)
  call boundary_peak_dp(s,c%max_dens)
  call boundary_peak_int(s,c%new_peak)
  call boundary_peak_int(s,alive)

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
           if (c%sparse_saddle_dens%maxval(ipeak)>0)then
              relevance_peak=c%max_dens(ipeak)/c%sparse_saddle_dens%maxval(ipeak)
           else
              relevance_peak=c%max_dens(ipeak)/r%density_threshold
           end if
           c%relevance(ipeak)=relevance_peak
        else
           c%relevance(ipeak)=0
        endif
     end do

     ! Merge all peaks to deepest level
     do ilev=idepth-2,0,-1
        do ipeak=1,c%npeak
           if(c%lev_peak(ipeak)==ilev)then
              merge_to=c%new_peak(ipeak)
              call get_local_peak_id(s,merge_to,jpeak)
              c%new_peak(ipeak)=c%new_peak(jpeak)
           endif
        end do
        call build_peak_communicator(s)
        call boundary_peak_int(s,c%new_peak)
     end do

     ! Update flag1 field
     do itest=1,c%ntest
        igrid=c%grid(itest)
        ind=c%cell(itest)
        global_peak_id=m%grid(igrid)%flag1(ind)
        if (global_peak_id>0)then
           call get_local_peak_id(s,global_peak_id,ipeak)
           merge_to=c%new_peak(ipeak)
           call get_local_peak_id(s,merge_to,jpeak)
           m%grid(igrid)%flag1(ind)=merge_to
        end if
     end do
     call build_peak_communicator(s)

  endif

  if(action.EQ.'saddleden')then

     ! Compute peak index for the halo
     do ipeak=1,c%npeak
        c%ind_halo(ipeak)=c%new_peak(ipeak)
     end do
     call boundary_peak_int(s,c%ind_halo)
     do ilev=idepth-2,0,-1
        do ipeak=1,c%npeak
           if(c%lev_peak(ipeak)==ilev)then
              merge_to=c%ind_halo(ipeak)
              call get_local_peak_id(s,merge_to,jpeak)
              c%ind_halo(ipeak)=c%ind_halo(jpeak)
           endif
        end do
        call build_peak_communicator(s)
        call boundary_peak_int(s,c%ind_halo)
     end do

     ! Compute halo masses
     c%halo_mass=0
     c%n_cells_halo=0
     do ipeak=1,c%npeak
        merge_to=c%ind_halo(ipeak)
        call get_local_peak_id(s,merge_to,jpeak)
        c%halo_mass(jpeak)=c%halo_mass(jpeak)+c%clump_mass(ipeak)
        c%n_cells_halo(jpeak)=c%n_cells_halo(jpeak)+c%n_cells(ipeak)
     end do
     call build_peak_communicator(s)
     call virtual_peak_dp(s,c%halo_mass,'sum')
     call boundary_peak_dp(s,c%halo_mass)
     call virtual_peak_int(s,c%n_cells_halo,'sum')
     call boundary_peak_int(s,c%n_cells_halo)
     ! Assign back halo mass to peak halo mass
     do ipeak=1,c%npeak
        merge_to=c%ind_halo(ipeak)
        call get_local_peak_id(s,merge_to,jpeak)
        c%halo_mass(ipeak)=c%halo_mass(jpeak)
     end do

  endif
  end associate
end subroutine merge_clumps
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_max(s,imax,mat)
  use ramses_commons, only: ramses_t
  use sparse_matrix
  type(ramses_t)::s
  type(sparse_mat)::mat
  integer::imax
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! get maximum in i-th line by walking the linked list
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::current,i,icol

  do i=1,imax
     
     mat%maxval(i)=0
     mat%maxloc(i)=0
     
     ! walk the line...
     current=mat%first(i)
     do  while( current /= 0 )
        if(mat%maxval(i)<mat%val(current))then
           mat%maxval(i)=mat%val(current)
           icol=mat%col(current)
           if(icol<=s%c%npeak)then
              mat%maxloc(i)=s%c%npeak_cum(s%g%myid-1)+icol
           else
              mat%maxloc(i)=s%c%gkey(icol)
           endif
        end if
        current=mat%next(current)
     end do

  end do

end subroutine get_max
!################################################################
!################################################################
!################################################################
!################################################################
subroutine virtual_peak_int(s,xx,action)
  use amr_commons
  use ramses_commons, only: ramses_t
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  integer,dimension(1:s%c%npeak_max)::xx
  character(len=3)::action
  
#ifndef WITHOUTMPI
  integer,allocatable,dimension(:)::int_peak_send_buf,int_peak_recv_buf
  integer::ipeak,icpu,info,j
  integer,dimension(1:s%g%ncpu)::ipeak_alltoall

  associate(g=>s%g,c=>s%c)

  allocate(int_peak_send_buf(1:c%peak_send_tot))
  allocate(int_peak_recv_buf(1:c%peak_recv_tot))

  ipeak_alltoall=0
  do ipeak=c%npeak+1,c%hfree-1
     call get_local_peak_cpu(s,ipeak,icpu)
     ipeak_alltoall(icpu)=ipeak_alltoall(icpu)+1
     int_peak_send_buf(c%peak_send_oft(icpu)+ipeak_alltoall(icpu))=xx(ipeak)
  end do
  call MPI_ALLTOALLV(int_peak_send_buf,c%peak_send_cnt,c%peak_send_oft,MPI_INTEGER, &
       &             int_peak_recv_buf,c%peak_recv_cnt,c%peak_recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)

  select case (action)
  case('sum')
     do j=1,c%peak_recv_tot
        ipeak=c%peak_recv_buf(j)-c%npeak_cum(g%myid-1)
        xx(ipeak)=xx(ipeak)+int_peak_recv_buf(j)
     end do
  case('min')
     do j=1,c%peak_recv_tot
        ipeak=c%peak_recv_buf(j)-c%npeak_cum(g%myid-1)
        xx(ipeak)=MIN(xx(ipeak),int_peak_recv_buf(j))
     end do
  case('max')
     do j=1,c%peak_recv_tot
        ipeak=c%peak_recv_buf(j)-c%npeak_cum(g%myid-1)
        xx(ipeak)=MAX(xx(ipeak),int_peak_recv_buf(j))
     end do
  end select
  
  deallocate(int_peak_send_buf,int_peak_recv_buf)

  end associate

#endif
end subroutine virtual_peak_int
!################################################################
!################################################################
!################################################################
!################################################################
subroutine virtual_peak_dp(s,xx,action)
  use amr_commons, only: dp
  use ramses_commons, only: ramses_t
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  real(dp),dimension(1:s%c%npeak_max)::xx
  character(len=3)::action
  
#ifndef WITHOUTMPI
  real(kind=8),allocatable,dimension(:)::dp_peak_send_buf,dp_peak_recv_buf
  integer::ipeak,icpu,info,j
  integer,dimension(1:s%g%ncpu)::ipeak_alltoall

  associate(g=>s%g,c=>s%c)

  allocate(dp_peak_send_buf(1:c%peak_send_tot))
  allocate(dp_peak_recv_buf(1:c%peak_recv_tot))

  ipeak_alltoall=0
  do ipeak=c%npeak+1,c%hfree-1
     call get_local_peak_cpu(s,ipeak,icpu)
     ipeak_alltoall(icpu)=ipeak_alltoall(icpu)+1
     dp_peak_send_buf(c%peak_send_oft(icpu)+ipeak_alltoall(icpu))=xx(ipeak)
  end do
  call MPI_ALLTOALLV(dp_peak_send_buf,c%peak_send_cnt,c%peak_send_oft,MPI_DOUBLE_PRECISION, &
       &             dp_peak_recv_buf,c%peak_recv_cnt,c%peak_recv_oft,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,info)

  select case (action)
  case('sum')
     do j=1,c%peak_recv_tot
        ipeak=c%peak_recv_buf(j)-c%npeak_cum(g%myid-1)
        xx(ipeak)=xx(ipeak)+dp_peak_recv_buf(j)
     end do
  case('min')
     do j=1,c%peak_recv_tot
        ipeak=c%peak_recv_buf(j)-c%npeak_cum(g%myid-1)
        xx(ipeak)=MIN(xx(ipeak),dp_peak_recv_buf(j))
     end do
  case('max')
     do j=1,c%peak_recv_tot
        ipeak=c%peak_recv_buf(j)-c%npeak_cum(g%myid-1)
        xx(ipeak)=MAX(xx(ipeak),dp_peak_recv_buf(j))
     end do
  end select

  deallocate(dp_peak_send_buf,dp_peak_recv_buf)

  end associate

#endif
end subroutine virtual_peak_dp
!################################################################
!################################################################
!################################################################
!################################################################
subroutine virtual_saddle_max(s)
  use amr_commons
  use ramses_commons, only: ramses_t
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  
#ifndef WITHOUTMPI
  integer::info,icpu
  real(kind=8),allocatable,dimension(:)::dp_peak_send_buf,dp_peak_recv_buf
  integer,allocatable,dimension(:)::int_peak_send_buf,int_peak_recv_buf
  integer::ipeak,jpeak,j
  integer,dimension(1:s%g%ncpu)::ipeak_alltoall

  associate(g=>s%g,c=>s%c)

  allocate(int_peak_send_buf(1:c%peak_send_tot))
  allocate(int_peak_recv_buf(1:c%peak_recv_tot))
  allocate(dp_peak_send_buf(1:c%peak_send_tot))
  allocate(dp_peak_recv_buf(1:c%peak_recv_tot))

  ipeak_alltoall=0
  do ipeak=c%npeak+1,c%hfree-1
     call get_local_peak_cpu(s,ipeak,icpu)
     ipeak_alltoall(icpu)=ipeak_alltoall(icpu)+1
     dp_peak_send_buf(c%peak_send_oft(icpu)+ipeak_alltoall(icpu))=c%sparse_saddle_dens%maxval(ipeak)
     int_peak_send_buf(c%peak_send_oft(icpu)+ipeak_alltoall(icpu))=c%sparse_saddle_dens%maxloc(ipeak)
  end do
  call MPI_ALLTOALLV(dp_peak_send_buf,c%peak_send_cnt,c%peak_send_oft,MPI_DOUBLE_PRECISION, &
       &             dp_peak_recv_buf,c%peak_recv_cnt,c%peak_recv_oft,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,info)
  call MPI_ALLTOALLV(int_peak_send_buf,c%peak_send_cnt,c%peak_send_oft,MPI_INTEGER, &
       &             int_peak_recv_buf,c%peak_recv_cnt,c%peak_recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  do j=1,c%peak_recv_tot
     ipeak=c%peak_recv_buf(j)-c%npeak_cum(g%myid-1)
     if(c%sparse_saddle_dens%maxval(ipeak)<dp_peak_recv_buf(j))then
        c%sparse_saddle_dens%maxval(ipeak)=dp_peak_recv_buf(j)
        c%sparse_saddle_dens%maxloc(ipeak)=int_peak_recv_buf(j)
        call get_local_peak_id(s,int_peak_recv_buf(j),jpeak)
     endif
  end do

  deallocate(dp_peak_send_buf,dp_peak_recv_buf)
  deallocate(int_peak_send_buf,int_peak_recv_buf)

  end associate
#endif

end subroutine virtual_saddle_max
!################################################################
!################################################################
!################################################################
!################################################################
subroutine boundary_peak_int(s,xx)
  use amr_commons
  use ramses_commons, only: ramses_t
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  integer,dimension(1:s%c%npeak_max)::xx

#ifndef WITHOUTMPI
  integer,allocatable,dimension(:)::int_peak_send_buf,int_peak_recv_buf
  integer::ipeak,icpu,info,j
  integer,dimension(1:s%g%ncpu)::ipeak_alltoall

  associate(g=>s%g,c=>s%c)

  allocate(int_peak_send_buf(1:c%peak_send_tot))
  allocate(int_peak_recv_buf(1:c%peak_recv_tot))
  do j=1,c%peak_recv_tot
     ipeak=c%peak_recv_buf(j)-c%npeak_cum(g%myid-1)
     int_peak_recv_buf(j)=xx(ipeak)
  end do
  call MPI_ALLTOALLV(int_peak_recv_buf,c%peak_recv_cnt,c%peak_recv_oft,MPI_INTEGER, &
       &             int_peak_send_buf,c%peak_send_cnt,c%peak_send_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  ipeak_alltoall=0
  do ipeak=c%npeak+1,c%hfree-1
     call get_local_peak_cpu(s,ipeak,icpu)
     ipeak_alltoall(icpu)=ipeak_alltoall(icpu)+1
     xx(ipeak)=int_peak_send_buf(c%peak_send_oft(icpu)+ipeak_alltoall(icpu))
  end do

  deallocate(int_peak_send_buf,int_peak_recv_buf)

  end associate

#endif
end subroutine boundary_peak_int
!################################################################
!################################################################
!################################################################
!################################################################
subroutine boundary_peak_dp(s,xx)
  use amr_commons
  use ramses_commons, only: ramses_t
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  real(dp),dimension(1:s%c%npeak_max)::xx

#ifndef WITHOUTMPI
  real(kind=8),allocatable,dimension(:)::dp_peak_send_buf,dp_peak_recv_buf
  integer::ipeak,icpu,info,j
  integer,dimension(1:s%g%ncpu)::ipeak_alltoall

  associate(g=>s%g,c=>s%c)

  allocate(dp_peak_send_buf(1:c%peak_send_tot))
  allocate(dp_peak_recv_buf(1:c%peak_recv_tot))

  do j=1,c%peak_recv_tot
     ipeak=c%peak_recv_buf(j)-c%npeak_cum(g%myid-1)
     dp_peak_recv_buf(j)=xx(ipeak)
  end do
  call MPI_ALLTOALLV(dp_peak_recv_buf,c%peak_recv_cnt,c%peak_recv_oft,MPI_DOUBLE_PRECISION, &
       &             dp_peak_send_buf,c%peak_send_cnt,c%peak_send_oft,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,info)
  ipeak_alltoall=0
  do ipeak=c%npeak+1,c%hfree-1
     call get_local_peak_cpu(s,ipeak,icpu)
     ipeak_alltoall(icpu)=ipeak_alltoall(icpu)+1
     xx(ipeak)=dp_peak_send_buf(c%peak_send_oft(icpu)+ipeak_alltoall(icpu))
  end do

  deallocate(dp_peak_send_buf,dp_peak_recv_buf)

  end associate
  
#endif
end subroutine boundary_peak_dp

!################################################################
!################################################################
!################################################################
!################################################################
subroutine analyze_peak_memory(s)
  use ramses_commons, only: ramses_t
  use clfind_commons
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::i,j
  integer,dimension(1:s%g%ncpu)::npeaks_all,npeaks_tot
  integer,dimension(1:s%g%ncpu)::hfree_all,hfree_tot
  integer,dimension(1:s%g%ncpu)::sparse_all,sparse_tot
  integer,dimension(1:s%g%ncpu)::coll_all,coll_tot

  npeaks_all=0
  npeaks_all(s%g%myid)=s%c%npeak
  coll_all=0
  coll_all(s%g%myid)=s%c%hcollision
  hfree_all=0
  hfree_all(s%g%myid)=s%c%hfree-s%c%npeak
  sparse_all=0
  sparse_all(s%g%myid)=s%c%sparse_saddle_dens%used
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npeaks_all,npeaks_tot,s%g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(coll_all,coll_tot,s%g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(hfree_all,hfree_tot,s%g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(sparse_all,sparse_tot,s%g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  npeaks_tot=npeaks_all
  coll_tot=coll_all
  hfree_tot=hfree_all
  sparse_tot=sparse_all
#endif
  if(s%g%myid==1)then
     write(*,*)'peaks per cpu'
     do i=0,s%g%ncpu-1,10
        write(*,'(256(I8,1X))')(npeaks_tot(j),j=i+1,min(i+10,s%g%ncpu))
     end do
     write(*,*)'ghost peaks per cpu'
     do i=0,s%g%ncpu-1,10
        write(*,'(256(I8,1X))')(hfree_tot(j),j=i+1,min(i+10,s%g%ncpu))
     end do
     write(*,*)'hash table collisions'
     do i=0,s%g%ncpu-1,10
        write(*,'(256(I8,1X))')(coll_tot(j),j=i+1,min(i+10,s%g%ncpu))
     end do
     write(*,*)'sparse matrix used'
     do i=0,s%g%ncpu-1,10
        write(*,'(256(I8,1X))')(sparse_tot(j),j=i+1,min(i+10,s%g%ncpu))
     end do
  end if
end subroutine analyze_peak_memory
!################################################################
!################################################################
!################################################################
!################################################################
subroutine compute_clump_properties(s)
  use amr_commons, only: dp,ndim
  use clfind_commons
  use ramses_commons, only: ramses_t
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  type(ramses_t)::s
  !----------------------------------------------------------------------------
  ! This subroutine performs a loop over all cells above the threshold and
  ! collects the  relevant information. After some MPI communications,
  ! all necessary peak-patch properties are computed
  !----------------------------------------------------------------------------
  integer::ipart,grid,peak_nr,ilevel,global_peak_id,ipeak,plevel,igrid,itest,icelln,idim,ind
  real(dp),dimension(1:ndim)::xcell
  real(dp)::dx_loc,tot_mass
  real(dp),dimension(1:ndim)::xcen
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

  c%min_dens=huge(zero)
  c%n_cells=0; c%n_cells_halo=0
  c%halo_mass=0d0; c%clump_mass=0d0; c%clump_vol=0d0
  c%center_of_mass=0d0
  c%peak_pos=0d0

  if(g%myid==1.and.r%verbose)write(*,*)'Entering compute clump properties'
  
  !--------------------------------------------------------
  ! Loop over local peaks and compute peak cell coordinates
  !--------------------------------------------------------
  do ipeak=1,c%npeak
    ilevel=c%peak_level(ipeak)
    igrid=c%peak_grid(ipeak)
    ind=c%peak_cell(ipeak)
    dx_loc=r%boxlen/2**ilevel
     ! Peak cell coordinates
    xcell(1)=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx_loc-m%skip(1)
#if NDIM>1
    xcell(2)=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx_loc-m%skip(2)
#endif
#if NDIM>2
    xcell(3)=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx_loc-m%skip(3)
#endif
    c%peak_pos(ipeak,1:ndim)=xcell(1:ndim)
  end do
#ifndef WITHOUTMPI
  ! Scatter results to all MPI domains
  do i=1,ndim
     call boundary_peak_dp(s,c%peak_pos(1,i))
  end do
#endif

  !--------------------------------------------------------------------------
  ! loop over all cells above the threshold
  !--------------------------------------------------------------------------
  do itest=1,c%ntest
    ilevel=c%level(itest)
    igrid=c%grid(itest)
    ind=c%cell(itest)
    global_peak_id=m%grid(igrid)%flag1(ind)

    if (global_peak_id /=0 ) then
      call get_local_peak_id(s,global_peak_id,peak_nr)

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

      ! Cell coordinates
      xcell(1)=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx_loc-m%skip(1)
#if NDIM>1
      xcell(2)=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx_loc-m%skip(2)
#endif
#if NDIM>2
      xcell(3)=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx_loc-m%skip(3)
#endif
      ! In case of periodic boundaries
      do idim=1,ndim
         if ((xcell(idim)-c%peak_pos(peak_nr,idim))>r%boxlen*0.5)xcell(idim)=xcell(idim)-r%boxlen
         if ((xcell(idim)-c%peak_pos(peak_nr,idim))<-r%boxlen*0.5)xcell(idim)=xcell(idim)+r%boxlen
      end do

      ! Clump center of mass location
      c%center_of_mass(peak_nr,1:ndim)=c%center_of_mass(peak_nr,1:ndim)+vol*d*xcell(1:ndim)

    end if
  end do
  call build_peak_communicator(s)

#ifndef WITHOUTMPI
  ! Collect results from all MPI domains
  call virtual_peak_int(s,c%n_cells,'sum')
  call virtual_peak_dp(s,c%min_dens,'min')
  call virtual_peak_dp(s,c%clump_mass,'sum')
  call virtual_peak_dp(s,c%clump_vol,'sum')
  do i=1,ndim
     call virtual_peak_dp(s,c%center_of_mass(1,i),'sum')
  end do
#endif

  ! Compute specific quantities
  do ipeak=1,c%npeak
     if (c%relevance(ipeak)>0..and.c%n_cells(ipeak)>0)then
        c%center_of_mass(ipeak,1:ndim)=c%center_of_mass(ipeak,1:ndim)/c%clump_mass(ipeak)
     end if
  end do

#ifndef WITHOUTMPI
  ! Scatter results to all MPI domains
  do i=1,ndim
     call boundary_peak_dp(s,c%peak_pos(1,i))
     call boundary_peak_dp(s,c%center_of_mass(1,i))
  end do
#endif

  ! Initialize halo mass to clump mass
  c%halo_mass=c%clump_mass

  ! Calculate total mass above threshold
  tot_mass=sum(c%clump_mass(1:c%npeak))

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(tot_mass,tot_mass_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  tot_mass=tot_mass_tot
#endif

  ! Compute further properties of the clumps
  do ipeak=1,c%npeak
     if (c%relevance(ipeak)>0..and.c%n_cells(ipeak)>0)then
        c%av_dens(ipeak)=c%clump_mass(ipeak)/c%clump_vol(ipeak)
     end if
  end do

#ifndef WITHOUTMPI
  ! Scatter results to all MPI domains
  call boundary_peak_dp(s,c%av_dens)
#endif

  end associate

end subroutine compute_clump_properties
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module clump_merger_module
