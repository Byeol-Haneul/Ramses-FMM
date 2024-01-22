module clump_finder_module
  use clump_merger_module
contains
subroutine m_clump_finder(pst,create_output,keep_alive)
  use output_clump_module
  use amr_parameters, only: flen
  use ramses_commons, only: pst_t
#ifdef GRAV
  use rho_fine_module, only: m_rho_fine
#endif
  implicit none
  type(pst_t)::pst
  logical::create_output,keep_alive
  !-----------------------------------------------------------------------
  ! This is the master routine for the RAMSES clump finder.
  !-----------------------------------------------------------------------  
  character(LEN=5)::nchar
  character(LEN=flen)::filename,filedir
  integer,dimension(1:flen/4)::input_array
  integer::dummy(1)

#if NDIM==3 && defined(GRAV)

  associate(r=>pst%s%r,g=>pst%s%g)

  if(r%verbose)write(*,*)'Entering clump_finder'
  
  !-----------------------------------------------------------------------
  ! Compute rho from gas density and/or dark matter and/or star particles
  !-----------------------------------------------------------------------
  call m_rho_fine(pst,r%levelmin)
  
  !------------------------------------------
  ! Find relevant peak patches and halos
  !------------------------------------------
  call r_clump_finder(pst,r%levelmin,1)
 
  !------------------------------------------
  ! Output clumps properties to file
  !------------------------------------------
  if(create_output.and.g%output_done)then
     call title(g%ifout-1,nchar)
     filename='output_'//TRIM(nchar)//'/'
     input_array=transfer(filename,input_array)
     if(r%verbose)write(*,*)'Writing clump and halo files'
     filename=TRIM(filename)//'clump_field_descriptor.txt'
     call file_descriptor_clump(r,filename)
     call r_output_clump(pst,input_array,flen/4,dummy,0)
  endif

  !------------------------------------------
  ! Deallocate all peak arrays if needed
  !------------------------------------------
  if(.not. keep_alive)then
     call r_deallocate_clump(pst,r%levelmin,1)
  endif

  end associate

#endif
end subroutine m_clump_finder
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_clump_finder(pst,ilevel,input_size)
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
     call r_clump_finder(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call clump_finder(pst%s)
  endif
  
end subroutine r_clump_finder
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine clump_finder(s)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s

#if NDIM==3 && defined(GRAV)

  !----------------------------------------------------------------------
  ! Count and collect all cells above the prescribed density threshold.
  ! We call these cell test particles for the watershed algorithm.
  !----------------------------------------------------------------------
  call collect_test(s)

  !----------------------------------------------------------------------
  ! Count and collect all density peaks.
  ! We also compute for each test particle the coordinates of its
  ! densest neighbor.
  !----------------------------------------------------------------------
  call collect_peak(s)

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
  ! Merge all neighboring peaks above the prescribed density
  ! threshold into halos, only if their saddle point density is larger
  ! that the prescribed saddle density threshold.
  !----------------------------------------------------------------------
  if(s%r%saddle_threshold>0)then
     call merge_clumps(s,'saddleden')
  endif

#endif
end subroutine clump_finder
!################################################################
!################################################################
!################################################################
!################################################################
#if NDIM==3 && defined(GRAV)
subroutine collect_test(s)
  use amr_parameters, only: twotondim,ndim,dp
  use ramses_commons, only: ramses_t
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  !==================================================================
  ! This is the clump finder routine for collecting test particles
  ! also known as all cells above the prescribed density threshold.
  ! Count number of test particles and share info across processors
  ! Written by Ziyong Wu (mini-ramses version December 2023).
  !==================================================================
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:s%g%ncpu)::ntest_cpu_all
#endif

  integer::ind,igrid,idim,icpu,ngrid,nleaf,now_level,next_level
  integer::istep,ipart,jpart,ip,itest
  integer::ilevel,i,levelmin_part
  integer,dimension(1:s%g%ncpu)::ntest_cpu
  integer(kind=8),dimension(0:s%g%ncpu)::ntest_cum
  logical::verbose_all=.false.
  integer::action,ivar_clump
  logical::ok
  real(kind=8)::dx,vol
  real(kind=8)::d,dx_loc
  real(dp),allocatable,dimension(:)::dens
  integer,allocatable,dimension(:)::isort
  integer,allocatable,dimension(:)::iswap

#ifdef GRAV

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c)

  if(g%myid==1.and.r%verbose)write(*,*)'Entering collect_test'

  !---------------------------------------------------------
  ! Count cells above threshold. We name them test particles
  !---------------------------------------------------------
  c%ntest = 0
  do ilevel=r%levelmin,r%nlevelmax ! Loop over levels
     do igrid=m%head(ilevel),m%tail(ilevel) ! Loop over grids
        do ind=1,twotondim ! Loop over cells
           ok = .not. m%grid(igrid)%refined(ind) ! Select leaf cells
           d = m%grid(igrid)%rho(ind)
           ok = ok .and. d > r%density_threshold
           m%grid(igrid)%flag1(ind) = 0
           if(ok)then
              c%ntest=c%ntest+1
           endif
        end do
     end do
  end do
  
  !-------------------------------------------------
  ! Compute number of test particles across all CPUs
  !-------------------------------------------------
  ntest_cpu=0
  ntest_cpu(g%myid)=c%ntest

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntest_cpu,ntest_cpu_all,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntest_cpu=ntest_cpu_all
#endif
  ntest_cum=0
  do icpu=1,g%ncpu
     ntest_cum(icpu)=ntest_cum(icpu-1)+int(ntest_cpu(icpu),kind=8)
  end do
  c%ntest_tot=ntest_cum(g%ncpu)
  if(g%myid==1)then
     if(c%ntest_tot.gt.0.and.(r%verbose.or.r%clump_info))then
        write(*,'(" Total number of cells above threshold=",I12)')c%ntest_tot
     endif
  end if

  if (c%ntest>0) then

     !------------------------------------------
     ! Allocate arrays for cells above threshold
     !------------------------------------------
     allocate(c%grid(c%ntest),c%cell(c%ntest),c%level(c%ntest),c%hash(1:c%ntest,0:ndim))
     c%grid=0;  c%cell=0; c%level=0; c%hash=0

     !---------------------------------------------------
     ! Compute and store arrays for cells above threshold
     !---------------------------------------------------
     allocate(dens(c%ntest))
     itest=0
     do ilevel=r%levelmin,r%nlevelmax ! Loop over levels
        do igrid=m%head(ilevel),m%tail(ilevel) ! Loop over grids
           do ind=1,twotondim ! Loop over cells
              ok=.not.m%grid(igrid)%refined(ind) ! Select leaf cells
              d=m%grid(igrid)%rho(ind)
              ok=ok.and.d>r%density_threshold
              if(ok)then
                 itest=itest+1
                 dens(itest)=d
                 c%grid(itest)=igrid
                 c%cell(itest)=ind
                 c%level(itest)=ilevel
              endif
           end do
        end do
     end do
     
     !-----------------------------------------------------------------------
     ! Sort cells above threshold according to their density
     !-----------------------------------------------------------------------
     allocate(isort(c%ntest))
     do i=1,c%ntest
        dens(i)=-dens(i)
        isort(i)=i
     end do
     call quick_sort_dp(dens,isort,c%ntest)
     deallocate(dens)
     ! Swap arrays to sort them permanently
     allocate(iswap(1:c%ntest))
     do i=1,c%ntest
        itest=isort(i)
        iswap(i)=c%grid(itest)
     end do
     c%grid=iswap
     do i=1,c%ntest
        itest=isort(i)
        iswap(i)=c%cell(itest)
     end do
     c%cell=iswap
     do i=1,c%ntest
        itest=isort(i)
        iswap(i)=c%level(itest)
     end do
     c%level=iswap
     deallocate(iswap,isort)
     
  endif

  end associate
#endif  
  
end subroutine collect_test
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine collect_peak(s)
  use amr_parameters, only: twotondim,ndim
  use amr_commons, only:oct,nbor
  use ramses_commons, only: ramses_t
  use cache_commons
  use cache
  use nbors_utils
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  !===================================================================
  ! This is the clump finder routine for collecting densest
  ! neighbors. Only scan cells above the density threshold.
  ! Density peaks are cells without densest neighbors.
  ! Count number of density peaks and share info across processors.
  ! Store the hash key of the densest neighbor for later usage.
  ! Written by Ziyong Wu (mini-ramses version December 2023).
  !==================================================================
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:s%g%ncpu)::npeak_cpu_all
#endif
  type(msg_int4_small_realdp)::dummy_int4_small_realdp
  type(oct),pointer::gridn
  integer::ilevel
  integer::icpu,next_level,now_level,icelln,idim,igrid,ind,itest,j,k
  integer::ipart,jpart,ip,i,icellp,icellpm,ipeak
  integer,dimension(1:s%g%ncpu)::npeak_cpu
  integer,dimension(1:ndim)::ckey,ckey_nbor
  integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
  real(dp),dimension(1:ndim)::xcen,xnei
  integer, parameter::nSnei=56
  type(nbor),dimension(1:nSnei) :: grid_nbor
  integer(kind=8),dimension(1:nSnei)::icell_nbor,level_nbor
  real(dp),dimension(1:ndim,1:nSnei)::xSnei
  real(dp)::dens_nbor,density_max,x,y,z
  logical::ok,ok_peak

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c)    

  if(g%myid==1.and.r%verbose)write(*,*)'Entering collect_peak'

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
           !             if((i==1.or.i==4).and.(j==1.or.j==4).and.(k==1.or.k==4)) ok=.false. ! edge
           if((i==2.or.i==3).and.(j==2.or.j==3).and.(k==2.or.k==3)) ok=.false. ! centre
           if(ok)then
              ind = ind+1
              x = (i-1)+0.5d0 - 2
              y = (j-1)+0.5d0 - 2
              z = (k-1)+0.5d0 - 2
              xSnei(1,ind) = x/2d0
              xSnei(2,ind) = y/2d0
              xSnei(3,ind) = z/2d0
           endif
        enddo
     enddo
  enddo
  
  !----------------------------------------
  ! Compute hash key of densest neighbor
  !----------------------------------------
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
       hilbert=m%domain,pack_size=storage_size(dummy_int4_small_realdp)/32,&
       pack=pack_fetch_saddle, unpack=unpack_fetch_saddle)
  
  c%npeak=0
  do itest=1,c%ntest
     ilevel=c%level(itest)
     igrid=c%grid(itest)
     ind=c%cell(itest)
     
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
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
        
        ! If missing, get neighboring cell at ilevel-1
        if(.not.associated(gridn))then
           call unlock_cache(s,gridn)
           ckey_nbor(1:ndim)=int(xnei(1:ndim)/2.0)
           hash_nbor(0)=ilevel
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
           
           ! If refined, get neighboring cell at ilevel+1
        else if (gridn%refined(icelln))then
           call unlock_cache(s,gridn)
           ckey_nbor(1:ndim)=int(xnei(1:ndim)*2.0)
           hash_nbor(0)=ilevel+2
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
        endif
        
        grid_nbor(j)%p => gridn
        icell_nbor(j) = icelln
        level_nbor(j) = hash_nbor(0)
        
     end do
     
     density_max=1.0001*m%grid(igrid)%rho(ind)
     ok_peak=.true.
     
     do j=1,nSnei
        gridn => grid_nbor(j)%p ! Gather neighboring grid
        icelln = icell_nbor(j)
        dens_nbor = gridn%rho(icelln)
        if(dens_nbor > density_max)then
           ok_peak=.false.
           density_max=dens_nbor
           ! Store hash key of densest neighbor
           c%hash(itest,0)=level_nbor(j)
           c%hash(itest,1)=2*gridn%ckey(1)+MOD((icelln-1)  ,2)
#if NDIM>1
           c%hash(itest,2)=2*gridn%ckey(2)+MOD((icelln-1)/2,2)
#endif
#if NDIM>2
           c%hash(itest,3)=2*gridn%ckey(3)+MOD((icelln-1)/4,2)
#endif
        endif
     end do
     
     if(ok_peak)then
        c%npeak=c%npeak+1
        c%hash(itest,0:ndim)=0
     endif
     
     ! Unlock neighboring grids
     do j=1,nSnei
        gridn => grid_nbor(j)%p
        call unlock_cache(s,gridn)
     end do
     
  end do

  call close_cache(s,m%grid_dict)
  
  !------------------------------------------------
  ! Compute total number of peaks across all CPUs
  ! Determine offset for global peak IDs
  !------------------------------------------------
  allocate(c%npeak_cum(0:g%ncpu))
  npeak_cpu=0
  npeak_cpu(g%myid)=c%npeak
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npeak_cpu,npeak_cpu_all,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  npeak_cpu=npeak_cpu_all
#endif
  c%npeak_cum=0
  do icpu=1,g%ncpu
     c%npeak_cum(icpu)=c%npeak_cum(icpu-1)+int(npeak_cpu(icpu),kind=8)
  end do
#ifdef WITHOUTMPI
  c%npeak_tot=int(c%npeak,kind=8)
#else
  c%npeak_tot=c%npeak_cum(g%ncpu)
#endif
  if (g%myid==1.and.c%npeak_tot>0) &
       & write(*,'(" Total number of density peaks found=",I10)')c%npeak_tot
  
  ! Compute the size of the peak-based arrays
  c%npeak_max=MAX(4*MAXVAL(npeak_cpu),1000)
  allocate(c%peak_cell(c%npeak_max))
  allocate(c%peak_grid(c%npeak_max))
  allocate(c%max_dens(c%npeak_max))
  c%max_dens=0d0; c%peak_cell=0; c%peak_grid=0

end associate

end subroutine collect_peak
!################################################################
!################################################################
!################################################################
!################################################################
subroutine collect_patch(s)
  use amr_parameters, only: twotondim,ndim
  use amr_commons,only: oct
  use ramses_commons, only: ramses_t
  use cache_commons
  use cache
  use nbors_utils
  use boundaries, only: init_bound_flag
  use marshal, only: pack_fetch_flag, unpack_fetch_flag
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  !==================================================================
  ! This is the clump finder routine for segmenting the density
  ! field into peak patches around each density peak.
  ! - loop over cells in descending density order
  ! - propagate peak id from densest neighbor
  ! - nmove is the number of peak id's passed along
  ! - done when nmove_tot=0 (for single core, only one sweep is necessary)
  ! Written by Ziyong Wu (mini-ramses version December 2023).
  !==================================================================
#ifndef WITHOUTMPI
  integer::info
  integer(kind=8)::nmove_all,nzero_all
#endif
  integer(kind=8),dimension(0:ndim)::hash_nbor
  type(msg_twin_realdp)::dummy_int4
  type(oct),pointer::gridn
  integer::icelln,igrid,ind,ipeak,istep,itest,nmove,nzero
  integer(kind=8)::nmove_tot,nzero_tot

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c)

  if(g%myid==1.and.r%verbose)write(*,*)'Entering collect_patch'

  !----------------------------------------------------------------------
  ! Flag peaks with global peak id using flag1 array
  !----------------------------------------------------------------------
  do igrid=1,r%ngridmax
     m%grid(igrid)%flag1(1:twotondim)=0
  end do
  ipeak = 0
  do itest=1,c%ntest
     if(c%hash(itest,0)==0)then
        ipeak=ipeak+1
        igrid=c%grid(itest)
        ind=c%cell(itest)
        m%grid(igrid)%flag1(ind)=ipeak+c%npeak_cum(g%myid-1)
        c%peak_grid(ipeak)=igrid
        c%peak_cell(ipeak)=ind
        c%max_dens(ipeak)=m%grid(igrid)%rho(ind)
     endif
  end do
  
  !-----------------------------------------------------------------------
  ! Determine peak-patches around each peak
  !-----------------------------------------------------------------------
  if (g%myid==1.and.c%ntest_tot>0)write(*,*)'Finding peak patches'
  nmove_tot=1
  istep=0
  do while (nmove_tot.gt.0)
     
     call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
          hilbert=m%domain,pack_size=storage_size(dummy_int4)/32,&
          pack=pack_fetch_flag,unpack=unpack_fetch_flag,bound=init_bound_flag)
     
     nmove=0
     nzero=0
     do itest=1,c%ntest
        hash_nbor=c%hash(itest,0:ndim)
        if(hash_nbor(0)>0)then
           igrid=c%grid(itest)
           ind=c%cell(itest)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.false.,fetch_cache=.true.)
           if(m%grid(igrid)%flag1(ind).ne.gridn%flag1(icelln))nmove=nmove+1
           m%grid(igrid)%flag1(ind)=gridn%flag1(icelln)
           if(m%grid(igrid)%flag1(ind).eq.0)nzero=nzero+1
        endif
     end do
     
     call close_cache(s,m%grid_dict)
     
     istep=istep+1
     nmove_tot=int(nmove,kind=8)
     nzero_tot=int(nzero,kind=8)
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(nmove_tot,nmove_all,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(nzero_tot,nzero_all,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
     nmove_tot=nmove_all
     nzero_tot=nzero_all
#endif
     if(c%ntest_tot>0.and.r%verbose.and.g%myid==1)write(*,*)"istep=",istep,"nmove=",nmove_tot
  end do
  
  end associate
    
end subroutine collect_patch
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine collect_saddle(s)
  use amr_parameters, only: twotondim,ndim
  use amr_commons, only:oct,nbor
  use ramses_commons, only: ramses_t
  use cache_commons
  use cache
  use nbors_utils
  use sparse_matrix
  implicit none
  type(ramses_t)::s
  !===================================================================
  ! This is the clump finder routine for collecting saddle points
  ! between each neighboring patches. Store the saddle points
  ! in the saddle point matrix in sparse matrix format.
  ! Written by Ziyong Wu (mini-ramses version December 2023).
  !==================================================================
  type(msg_int4_small_realdp)::dummy_int4_small_realdp
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
           !             if((i==1.or.i==4).and.(j==1.or.j==4).and.(k==1.or.k==4)) ok=.false. ! edge
           if((i==2.or.i==3).and.(j==2.or.j==3).and.(k==2.or.k==3)) ok=.false. ! centre
           if(ok)then
              ind = ind+1
              x = (i-1)+0.5d0 - 2
              y = (j-1)+0.5d0 - 2
              z = (k-1)+0.5d0 - 2
              xSnei(1,ind) = x/2d0
              xSnei(2,ind) = y/2d0
              xSnei(3,ind) = z/2d0
           endif
        enddo
     enddo
  enddo
  
  !----------------------------------------
  ! Compute hash key of densest neighbor
  !----------------------------------------
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
       hilbert=m%domain,pack_size=storage_size(dummy_int4_small_realdp)/32,&
       pack=pack_fetch_saddle, unpack=unpack_fetch_saddle)
  
  do itest=1,c%ntest
     ilevel=c%level(itest)
     igrid=c%grid(itest)
     ind=c%cell(itest)

     peak_cen = m%grid(igrid)%flag1(ind)
     dens_cen = m%grid(igrid)%rho(ind)
     
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
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
        
        ! If missing, get neighboring cell at ilevel-1
        if(.not.associated(gridn))then
           call unlock_cache(s,gridn)
           ckey_nbor(1:ndim)=int(xnei(1:ndim)/2.0)
           hash_nbor(0)=ilevel
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
           
           ! If refined, get neighboring cell at ilevel+1
        else if (gridn%refined(icelln))then
           call unlock_cache(s,gridn)
           ckey_nbor(1:ndim)=int(xnei(1:ndim)*2.0)
           hash_nbor(0)=ilevel+2
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
        endif
        
        grid_nbor(j)%p => gridn
        icell_nbor(j) = icelln
        
     end do
     
     do j=1,nSnei
        gridn => grid_nbor(j)%p ! Gather neighboring grid
        icelln = icell_nbor(j)

        peak_nbor = gridn%flag1(icelln)
        dens_nbor = gridn%rho(icelln)
        
        ok = peak_cen/=0
        ok = ok .and. peak_nbor/=0
        ok = ok .and. peak_cen/=peak_nbor
        dens_ave = 0.5*(dens_cen+dens_nbor)

        if(ok)then ! if all criteria met, replace saddle density array value
           call get_local_peak_id(s,peak_cen,ipeak)
           call get_local_peak_id(s,peak_nbor,jpeak)
           if (get_value(ipeak,jpeak,c%sparse_saddle_dens) < dens_ave)then
              call set_value(ipeak,jpeak,dens_ave,c%sparse_saddle_dens)
           end if
           if (get_value(jpeak,ipeak,c%sparse_saddle_dens) < dens_ave)then
              call set_value(jpeak,ipeak,dens_ave,c%sparse_saddle_dens)
           end if
        end if
     
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
#endif
end module clump_finder_module

