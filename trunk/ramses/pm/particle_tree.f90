!################################################################
!################################################################
!################################################################
!################################################################
! subroutine memory_sort_level(ind_com,np,ilevel,icpu)
!   use pm_commons
!   use amr_commons

!   ! sort all the particles in memory according to their level

!   implicit none
!   integer::np,icpu,ilevel
!   integer,dimension(1:nvector)::ind_com
  
!   integer::i,idim,igrid
!   integer,dimension(1:nvector),save::ind_list,ind_part
!   logical,dimension(1:nvector),save::ok=.true.
!   integer::current_property

!   ! Compute parent grid index
!   do i=1,np
!      igrid=emission(icpu,ilevel)%fp(ind_com(i),1)
!      ind_list(i)=emission(icpu,ilevel)%igrid(igrid)
!   end do

!   ! Add particle to parent linked list
!   call remove_free(ind_part,np)
!   call add_list(ind_part,ind_list,ok,np)

!   ! Scatter particle level and identity
!   do i=1,np
!      levelp(ind_part(i))=emission(icpu,ilevel)%fp(ind_com(i),2)
!      idp   (ind_part(i))=emission(icpu,ilevel)%fp(ind_com(i),3)
!   end do

!   ! Scatter particle position and velocity
!   do idim=1,ndim
!   do i=1,np
!      xp(ind_part(i),idim)=emission(icpu,ilevel)%up(ind_com(i),idim     )
!      vp(ind_part(i),idim)=emission(icpu,ilevel)%up(ind_com(i),idim+ndim)
!   end do
!   end do

!   current_property = twondim+1

!   ! Scatter particle mass
!   do i=1,np
!      mp(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
!   end do
!   current_property = current_property+1

! #ifdef OUTPUT_PARTICLE_POTENTIAL
!   ! Scatter particle phi
!   do i=1,np
!      ptcl_phi(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
!   end do
!   current_property = current_property+1
! #endif

! end subroutine memory_sort_level
! !################################################################
! !################################################################
! !################################################################
! !################################################################

! subroutine get_particle_levels

! end subroutine get_particle_levels
! !################################################################
! !################################################################
! !################################################################
! !################################################################

subroutine sort_particles(ilevel, use_histograms)
  use pm_commons,  only: npart, part_level_offset, &
                         nbins, bin_keys, part_hkey, &
                         part_ind_permutation
  use amr_commons, only: dp, myid, levelmin, nlevelmax, ncpu
  use sort,        only: lsd_radix_sort_particles, apply_particle_permutation
  use hilbert,     only: hilbert_for_particle 
  use particle_communication, only: build_communicator
  implicit none

  integer, intent(in) :: ilevel
  logical, intent(in) :: use_histograms
  
  integer :: ndata_remote, ndata_local, local_oft
  integer :: ilev, offset, np, ip
  integer, dimension(1:ncpu, 1:4) :: communicator
  integer, allocatable, dimension(:) :: refined
  
  ! This routine moves all particles that sit in refined cells at level ilevel
  ! to level ilevel + 1. It then sorts the remaining (true) ilevel particles 
  ! by hilbert key.

  offset = part_level_offset(ilevel)
  np = npart - part_level_offset(ilevel)

  ! Compute hilbert keys (probably move outside of this routine)
  call hilbert_for_particle(offset, npart - offset, 0, ilevel)
  
  ! Compute a permutation that sorts ALL particles starting from offset
  call lsd_radix_sort_particles(offset, npart - offset, ilevel, ilevel, .true.)

  if (use_histograms)then
     call compute_particle_histogram(offset, npart - offset)          
     call build_communicator(communicator, ndata_remote, &
                             nbins, ndata_local, local_oft, &
                             bin_keys(1:nbins, 2), bin_keys(1:nbins, 1), bin_keys(1:nbins, 0), &
                             ilevel)
     allocate(refined(1:nbins))     
     call communicate_refinements(communicator, ndata_remote, &
                                  nbins, ndata_local, local_oft, refined, &
                                  bin_keys(1:nbins, 2), bin_keys(1:nbins, 1), bin_keys(1:nbins, 0), &
                                  ilevel)
     call reshuffle_particles(ilevel, np, nbins, refined, use_histograms)
  else

     ! Need to apply particle permutation here already to have parts sorted
     ! in memory. Using the index insidet build_communicator and
     ! communicate_refinements is possible but will let the code deviate more
     ! from the histogrammed case.
     call apply_particle_permutation(offset, npart - offset, ilevel) 
     call build_communicator(communicator, ndata_remote, &
                             np, ndata_local, local_oft, &
                             part_hkey(offset + 1: offset + np, 2), &
                             part_hkey(offset + 1: offset + np, 1), &
                             part_hkey(offset + 1: offset + np, 0), ilevel)
     allocate(refined(1:np))
     call communicate_refinements(communicator, ndata_remote, &
                                  np, ndata_local, local_oft, refined, &
                                  part_hkey(offset + 1: offset + np, 2), &
                                  part_hkey(offset + 1: offset + np, 1), &
                                  part_hkey(offset + 1: offset + np, 0), ilevel)

     call reshuffle_particles(ilevel, np, np, refined, use_histograms)
     print*, myid, 'refined total on level ', ilevel, sum(refined)
  end if

  ! Compute NEW number of particles in ilevel
  np = part_level_offset(ilevel + 1) - part_level_offset(ilevel)  

  ! Re-sort remaining (ilevel particles)
  call lsd_radix_sort_particles(offset, np, ilevel, ilevel, .true.)
  call apply_particle_permutation(offset, np, ilevel)
  deallocate(refined)

!  if (ilevel == nlevelmax)then
!     do ilev=levelmin,nlevelmax
!        print*, 'nparts on (going out) ',myid, ilev, part_level_offset(ilev + 1) - part_level_offset(ilev)
!     end do
!  end if
  
end subroutine sort_particles
!################################################################
!################################################################
!################################################################
!################################################################
subroutine communicate_refinements(communicator, ndata_remote, ndata, ndata_local, ndata_local_oft, &
                                   refined, keys2, keys1, keys0, ilevel)

  use amr_commons,   only: ncpu, myid, bound_key_level, son, nvector, nlevelmax
  use particle_communication, only: part_data_to_domain_i8, part_data_to_domain_i4, &
                                    domain_data_to_part_i4
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in) ::  ilevel, ndata
  integer, intent(in) :: ndata_remote, ndata_local, ndata_local_oft
  integer(kind=8), dimension(1:ndata), intent(in) :: keys2, keys1, keys0
  integer, dimension(1:ncpu, 1:4), intent(in) :: communicator
  integer, dimension(1:ndata), intent(inout) :: refined
  
  ! This routine sorts particles between ilevel and ilevel + 1
  ! (formerly known as kill_tree_fine).

  ! The routine can only run after:
  ! call compute_particle_histogram(ilevel)
  ! call build_communicator(...)

  ! The "bins" here are the histogram bins and correspond to 
  ! actual grid cells. "local_... variables" denote properties in the
  ! MPI process which hosts the particles, while "remote_... variables" 
  ! are used for properties in the MPI process which hosts the corresponding
  ! leaf-cell.


  integer  :: idata, ioft, nd
  integer, dimension(1:nvector) :: dummy_int
  
  integer,         allocatable, dimension(:  ) :: remote_refined
  integer(kind=8), allocatable, dimension(:,:) :: remote_keys

  allocate(remote_refined(1:ndata_remote))
  allocate(remote_keys(1:ndata_remote, 0:2))

  call part_data_to_domain_i8(communicator, keys0, remote_keys(:,0))
  call part_data_to_domain_i8(communicator, keys1, remote_keys(:,1))
  call part_data_to_domain_i8(communicator, keys2, remote_keys(:,2))

  ! Probe local cells for refinement (abuse refined to store cell index)
  do ioft = ndata_local_oft, ndata_local_oft + ndata_local -1, nvector
     nd = min(ndata_local_oft + ndata_local - ioft, nvector) 
     call get_cell_index_from_hilbertkey(refined(ioft + 1 : ioft + nd), &
          dummy_int(1 : nd), keys2(ioft + 1: ioft + nd), &
          keys1(ioft + 1: ioft + nd), keys0(ioft + 1: ioft + nd), nd, ilevel)
  end do

  ! Mark data corresponding to refined cells
  do idata = ndata_local_oft + 1, ndata_local_oft + ndata_local
     if (son(refined(idata)) > 0)then
        refined(idata) = 1
     else
        refined(idata) = 0
     end if
  end do
  
  ! Probe remote cells for refinement (abuse remote_refined to store cell index) 
  do ioft = 0, ndata_remote - 1 , nvector
     nd = min(ndata_remote - ioft, nvector) 
     call get_cell_index_from_hilbertkey(remote_refined(ioft + 1 : ioft + nd), &
          dummy_int(1 : nd), remote_keys(ioft + 1: ioft + nd,2), &
          remote_keys(ioft + 1: ioft + nd,1), remote_keys(ioft + 1: ioft + nd,0), nd, ilevel)
  end do

  ! Mark bins corresponding to refined cells
  do idata = 1, ndata_remote
     if (son(remote_refined(idata))>0)then
        remote_refined(idata) = 1
     else
        remote_refined(idata) = 0
     end if
  end do

  ! Send refinement information back
  call domain_data_to_part_i4(communicator, remote_refined, refined)

  deallocate(remote_refined)
  deallocate(remote_keys)
  
end subroutine communicate_refinements
!################################################################
!################################################################
!################################################################
!################################################################

!################################################################
!################################################################
!################################################################
!################################################################
subroutine reshuffle_particles(ilevel, np, ndata, refined, use_histograms)
  use pm_commons,     only: part_level_offset, part_ind_permutation, part_hkey, &
                            bin_keys, part_ind_permutation2
  use sort,           only: gt_3keys, apply_particle_permutation
  use amr_parameters, only: nlevelmax
  implicit none

  
  integer, intent(in) :: ilevel, np, ndata
  integer, dimension(1:ndata), intent(in) :: refined
  logical, intent(in) :: use_histograms
  
  ! Reshuffle particles in memory according to their level
  ! by applying a couting sort on the particles.

  integer  :: offset, ibin, ip, ipart
  integer  :: unrefined_pos, refined_pos
  logical  :: unrefined

  if (ilevel == nlevelmax) return
  if (ndata == 0) return
  
  offset = part_level_offset(ilevel)

  ! Find starting indices for refined particles
  unrefined_pos = offset; refined_pos = offset          

  if (use_histograms)then
     ibin = 1; unrefined = (refined(ibin) == 0)     
     do ip = offset + 1, offset + np  
        ipart = part_ind_permutation(ip)
        if (gt_3keys(part_hkey(ipart,0:2), bin_keys(ibin,0:2)))then
           ibin=ibin+1
        end if
        if (refined(ibin) == 0) then 
           refined_pos = refined_pos + 1
        end if
     end do
  else
     refined_pos = refined_pos + np - sum(refined)
  end if

  ! set "level boundary" in particle array and rearrange particles
  part_level_offset(ilevel+1) = refined_pos

  if (use_histograms)then
     ibin = 1; unrefined = (refined(ibin) == 0)
     do ip = offset + 1, offset + np
        ipart = part_ind_permutation(ip)
        if (gt_3keys(part_hkey(ipart,0:2), bin_keys(ibin,0:2))) then
           ibin = ibin + 1
        end if
        if (refined(ibin) == 1)then
           refined_pos = refined_pos + 1
           part_ind_permutation2(refined_pos) = ipart
        else
           unrefined_pos = unrefined_pos + 1
           part_ind_permutation2(unrefined_pos) = ipart
        end if
     end do
  else
     do ip = offset + 1, offset + np
        ipart = part_ind_permutation(ip)
        if (refined(ipart - offset) == 1)then
           refined_pos = refined_pos + 1
           part_ind_permutation2(refined_pos) = ipart
        else
           unrefined_pos = unrefined_pos + 1
           part_ind_permutation2(unrefined_pos) = ipart
        end if
     end do
  end if
  
  part_ind_permutation(offset + 1:offset + np) = &
       part_ind_permutation2(offset +1:offset + np)
  call apply_particle_permutation(offset, np, ilevel)

end subroutine reshuffle_particles
!################################################################
!################################################################
!################################################################
!################################################################
! subroutine part_to_cell_i8(part_array,sortind)
!   use pm_parameters, only: npartmax
!   implicit none
! #ifndef WITHOUTMPI
!   include 'mpif.h'
! #endif
  
!   integer(kind=8),dimension(1:npartmax)::part_array, sortind
  
!   ! Communication routine that sends a particle-based quantity
!   ! to the MPI domain that owns the respective cell
  
  
  
  

! #ifndef WITHOUTMPI
!   do j=1,part_recv_tot
!      ipart=part_recv_buf(j)-ipart_start(myid)
!      int_part_recv_buf(j)=xx(ipart)
!   end do
!   call MPI_ALLTOALLV(int8_part_recv_buf,part_recv_cnt,part_recv_oft,MPI_INTEGER, &
!        &             int8_part_send_buf,part_send_cnt,part_send_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  
  
!   deallocate(int8_part_send_buf,int8_part_recv_buf)
! #endif
  
! end subroutine part_to_cell_i8


! !################################################################
! !################################################################
! !################################################################
! !################################################################
!subroutine part_to_cell_dp
  

!  implicit none
  

!end subroutine part_to_cell_dp
! !################################################################
! !################################################################
! !################################################################
! !################################################################

! !################################################################
! !################################################################
! !################################################################
! !################################################################
!subroutine get_cell_index_from_hilbert(cell_index,hkey0,hkey1,hkey2,np,ilevel)
!   use amr_commons
!   ! assuming the last 3 digits are in cartesian order
!   do i=1,np
!      grid_index=hash_get(hkey0(i))
!      ind=mod(hkey0(i),8)
!      iskip=ncoarse+(ind-1)*ngridmax
!      cell_index(i)=iskip+ind_grid
!   end do
! end subroutine get_cell_index

subroutine get_cell_index_from_hilbertkey(cell_index,cell_levl,hilbert_key2,hilbert_key1,hilbert_key0,np,ilevel)
  use amr_commons, only: nlevelmax, nvector, myid, ncoarse, ngridmax, xg
  use hilbert,     only: hilbert3d_reverse
  implicit none
  integer, intent(in)::np,ilevel
  integer(kind=8),dimension(1:nvector)::x,y,z
  integer(kind=8),dimension(1:nvector)::hilbert_key2,hilbert_key1,hilbert_key0
  integer,dimension(1:nvector)::cell_levl, cell_index
  integer,dimension(1:nvector)::cell_levl2, cell_index2
  integer :: i
  call hilbert3d_reverse(x,y,z,hilbert_key2,hilbert_key1,hilbert_key0,ilevel,np)
  call get_cell_index_from_cartesian_hash(cell_index,cell_levl,x,y,z,ilevel,np,ilevel)
  ! call get_cell_index_from_cartesian(cell_index2,cell_levl2,x,y,z,ilevel,np,ilevel)

  ! do i = 1, np
  !    if (cell_index(i) .ne. cell_index2(i) .or. cell_levl(i) .ne. cell_levl2(i) ) then
  !       print*, 'problem in cind:', i
  !       print*, x(i),y(i),z(i)
  !       print*, cell_index(i), cell_index2(i)
  !       print*,cell_levl(i), cell_levl2(i), ilevel
  !       print*,'xg',xg(mod(cell_index2(i) - ncoarse,ngridmax),1:3) 
  !       stop
  !    end if
  ! end do
     
end subroutine get_cell_index_from_hilbertkey

subroutine get_cell_index_from_cartesian(cell_index,cell_levl,xx,yy,zz,ilevel,n,bit_length)
  use amr_commons
  implicit none

  integer, intent(in)::n,ilevel,bit_length
  integer,dimension(1:nvector)::cell_index,cell_levl
  integer(kind=8),dimension(1:nvector)::xx,yy,zz
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  integer::i,j,ind,iskip,igrid,ind_cell,igrid0
  integer(kind=8)::ii,jj,kk

  if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  else
     write(*,*)"nx=ny=nz != 1,3 is not supported."
     stop
  end if

  if (bit_length>21)then
     print*, 'bit length too big for now'
  end if
  
  ind_cell=0
  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
  do i=1,n
     igrid=igrid0
     do j=1,ilevel 
        ii=ISHFT(xx(i),-bit_length+j)
        jj=ISHFT(yy(i),-bit_length+j)
        kk=ISHFT(zz(i),-bit_length+j)
        ii=mod(ii,2)
        jj=mod(jj,2)
        kk=mod(kk,2)
        ind=1+ii+2*jj+4*kk
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        igrid=son(ind_cell)
        if(igrid==0.or.j==ilevel)exit
     end do
     cell_index(i)=ind_cell
     cell_levl(i)=j
  end do
end subroutine get_cell_index_from_cartesian


subroutine get_cell_index_from_cartesian_hash(cell_index,cell_levl,xx,yy,zz,ilevel,n,bit_length)
  use amr_commons
  use hash, only: hash_get
  implicit none

  integer, intent(in)::n,ilevel,bit_length
  integer,intent(inout), dimension(1:nvector)::cell_index,cell_levl
  integer(kind=8),intent(in),dimension(1:nvector)::xx,yy,zz
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  integer :: i
  integer(kind=8), dimension(0:ndim) :: key
  
  if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  else
     write(*,*)"nx=ny=nz != 1,3 is not supported."
     stop
  end if

  if (bit_length>21)then
     print*, 'bit length too big for now'
  end if

  cell_levl(1:n) = ilevel

  ! Probe for cells starting from ilevel, if cell not present, try coarser
  do i = 1, n
     key(0) = ilevel
     key(1) = xx(i)
     key(2) = yy(i)
     key(3) = zz(i)
     cell_index(i) = hash_get(cell_dict, key)

     do while (cell_index(i) == 0 .and. cell_levl(i) > 1)
        cell_levl(i) = cell_levl(i) - 1
        key(0) = cell_levl(i)
        key(1) = ISHFT(key(1), -1)
        key(2) = ISHFT(key(2), -1)
        key(3) = ISHFT(key(3), -1)
        cell_index(i) = hash_get(cell_dict, key)
     end do
  end do

  ! Do check if all went well
  do i = 1, n
     if (cell_index(i) == 0) then
        write(*,*)"Problem in get_cell_index_from_cartesian"
        stop
     end if
  end do
end subroutine get_cell_index_from_cartesian_hash

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_cell_index(cell_index,cell_levl,xpart,ilevel,n)
  use amr_commons
  implicit none

  integer::n,ilevel
  integer,dimension(1:nvector)::cell_index,cell_levl
  real(dp),dimension(1:nvector,1:3)::xpart

  !----------------------------------------------------------------------------
  ! This routine returns the index and level of the cell, (at maximum level
  ! ilevel), in which the input the position specified by xpart lies
  !----------------------------------------------------------------------------

  real(dp)::xx,yy,zz
  integer::i,j,ii,jj,kk,ind,iskip,igrid,ind_cell,igrid0

  if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  else
     write(*,*)"nx=ny=nz != 1,3 is not supported."
     call clean_stop
  end if

  ind_cell=0
  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
  do i=1,n
     xx = xpart(i,1)/boxlen + (nx-1)/2.0
     yy = xpart(i,2)/boxlen + (ny-1)/2.0
     zz = xpart(i,3)/boxlen + (nz-1)/2.0

     if(xx<0.)xx=xx+dble(nx)
     if(xx>dble(nx))xx=xx-dble(nx)
     if(yy<0.)yy=yy+dble(ny)
     if(yy>dble(ny))yy=yy-dble(ny)
     if(zz<0.)zz=zz+dble(nz)
     if(zz>dble(nz))zz=zz-dble(nz)

     igrid=igrid0
     do j=1,ilevel 
        ii=1; jj=1; kk=1
        if(xx<xg(igrid,1))ii=0
        if(yy<xg(igrid,2))jj=0
        if(zz<xg(igrid,3))kk=0
        ind=1+ii+2*jj+4*kk
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        igrid=son(ind_cell)
        if(igrid==0.or.j==ilevel)exit
     end do
     cell_index(i)=ind_cell
     cell_levl(i)=j
  end do
end subroutine get_cell_index

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine compute_particle_histogram(offset, np)
  use pm_commons
  use amr_commons
  use sort,        only: gt_3keys
  implicit none
  integer, intent(in) :: offset, np

  !----------------------------------------------------------------------------
  ! This routine computes a particle histogram for np particles in memory, 
  ! starting from offset+1 to offset+np. There must be a precomputed array 
  ! part_ind_permutation which sorts the particles by hilbert key.
  !----------------------------------------------------------------------------

  integer,                         save :: ibin, ipart, ip
  integer(kind=8), dimension(0:2), save :: current_bin_key


  ! if there is nothing to do...
  nbins = 0
  if (.not. np > 0)return
  
  ! Count the number of bins
  nbins = 1
  current_bin_key(0:2) = part_hkey(part_ind_permutation(offset + 1),0:2)
  do ip = offset + 2, offset + np
     ipart = part_ind_permutation(ip)
     if (gt_3keys(part_hkey(ipart,0:2), current_bin_key(0:2)))then
        nbins=nbins+1
        current_bin_key(0:2) = part_hkey(ipart,0:2)
     end if
  end do
  
  ! Allocate histograms
  if(allocated(bin_count))then
     deallocate(bin_keys)
     deallocate(bin_count)
     deallocate(bin_start_offset)
     deallocate(bin_mass)
  end if
  if (.not. allocated(bin_keys))then
     allocate(bin_keys(nbins,0:2))
     allocate(bin_count(nbins))
     allocate(bin_start_offset(nbins+1))
     allocate(bin_mass(nbins))
  end if
  bin_mass = 0; bin_count = 0.d0

  ! Label every bin by a key and sum up the particles per bin, store 
  ! the offset of the first particle in each bin in the particle array

  ! First particle in first bin
  ibin=1
  bin_keys(ibin, 0:2) = part_hkey(part_ind_permutation(offset + 1), 0:2)  
  bin_count(ibin) = 1.d0
  bin_start_offset(ibin) = offset

  ! All other particles/bins
  current_bin_key(0:2) = part_hkey(part_ind_permutation(offset + 1), 0:2)
  do ip = offset + 2, offset + np
     ipart = part_ind_permutation(ip)
     if (gt_3keys(part_hkey(ipart,0:2), current_bin_key(0:2)))then
        ibin = ibin + 1
        bin_start_offset(ibin) = ip - 1 
        bin_keys(ibin,0:2) = part_hkey(ipart, 0:2)
        current_bin_key(0:2) = part_hkey(ipart, 0:2)
     end if
     bin_count(ibin) = bin_count(ibin) + 1.d0
  end do
  bin_start_offset(nbins+1) = offset + np

end subroutine compute_particle_histogram





subroutine count_parts
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  ! ugly routine to count all particles in the simulation level by level
  ! used for debugging
  integer::ilevel
  integer::igrid,jgrid,i,ngrid,ncache
  integer::ig,ip,npart1,npart2,npart2_tot,icpu,info
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nlevelmax)::npts

  npts=0  
  do ilevel=1,nlevelmax
     npart2=0
     ! Loop over cpus
     do icpu=1,ncpu
        igrid=headl(icpu,ilevel)
        ig=0
        ip=0
        ! Loop over grids
        do jgrid=1,numbl(icpu,ilevel)
           npart1=numbp(igrid)  ! Number of particles in the grid
           npart2=npart2+npart1
           igrid=next(igrid)   ! Go to next grid                                                              
        
        end do
     end do
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(npart2,npart2_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
     npart2_tot=npart2
#endif          
     npts(ilevel)=npart2_tot
  end do
  if (myid==1)print*,'total'
  if (myid==1)print*,npts(1:nlevelmax)
  
  

   npts=0  
   do ilevel=1,nlevelmax
      npart2=0        
      
      ncache=active(ilevel)%ngrid
      do igrid=1,ncache,nvector
         ngrid=MIN(nvector,ncache-igrid+1)
         do i=1,ngrid
            ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
         end do
         do i=1,ngrid
            npart2=npart2+numbp(ind_grid(i))
         end do
        
      end do

#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(npart2,npart2_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
      npart2_tot=npart2
#endif          
      npts(ilevel)=npart2_tot
   end do
   if (myid==1)print*,'active'
   if (myid==1)print*,npts(1:nlevelmax)
   
   
end subroutine count_parts
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine check_sorted(offset, np)
  use pm_commons,   only : part_hkey, part_ind_permutation
  use sort,         only : ge_3keys
  use amr_commons,  only : myid
  implicit none
  integer, intent(in) :: offset, np

  !----------------------------------------------------------------------------
  ! Simple helper routine to check whether the np particles starting form offset
  ! are sorted by hilbert key
  !----------------------------------------------------------------------------
  logical,                         save :: ok
  integer,                         save :: ipart, ip
  integer(kind=8), dimension(0:2), save :: current_key

  ok =.true.
  
  !  current_key(0:2) = part_hkey(part_ind_permutation(offset + 1),0:2)
    current_key(0:2) = part_hkey(offset + 1,0:2)
  do ip = offset + 2, offset + np
     !ipart = part_ind_permutation(ip)
     ipart = ip
     if (.not. ge_3keys(part_hkey(ipart,0:2), current_key(0:2)))then
        ok=.false.
        print*, "Detected unsorted particles on process", myid
        print*, part_hkey(ipart,0),current_key(0)
        print*, part_hkey(ipart,1),current_key(1)
        print*, part_hkey(ipart,2),current_key(2)
     end if
     current_key(0:2) = part_hkey(ipart,0:2)
  end do
  if (.not. ok)stop
end subroutine check_sorted




subroutine write_ascii_parts
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  ! ugly routine to count all particles in the simulation level by level
  ! used for debugging
  integer::ilevel, ilun
  integer::igrid,jgrid,i,ngrid,ncache, ipart, jpart, next_part
  integer::ig,ip,npart1,npart2,npart2_tot,icpu,info
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nlevelmax)::npts
  character(len=80)::filename
  character(LEN=5)::nchar
  
  ilun=myid+100
  
  call title(ilun,nchar)
  filename='OLDParts'//nchar
 
  open(unit=ilun,file=filename,form='formatted')
  
  do ilevel=levelmin,levelmin
     ! Loop over cpus
     do icpu=1,ncpu
        igrid=headl(icpu,ilevel)
        ! Loop over grids
        do jgrid=1,numbl(icpu,ilevel)
           npart1=numbp(igrid)  ! Number of particles in the grid
           if(npart1>0)then
              ipart=headp(igrid)
              ! Loop over particles
              do jpart=1,npart1
                 ! Save next particle   <--- Very important !!!
                 next_part=nextp(ipart)
                 write(ilun,"(6(F20.15),2(I10))")xp(ipart,1:3),vp(ipart,1:3),idp(ipart), levelp(ipart) 
                 ipart=next_part  ! Go to next particle
              end do
           endif
           igrid=next(igrid)   ! Go to next grid                                                              
        end do
     end do
  end do
  close(ilun)

  filename='NEWParts'//nchar
  open(unit=ilun,file=filename,form='formatted')
  do ipart=1,npart
     write(ilun,"(6(F20.15),2(I10))")xp(ipart,1:3),vp(ipart,1:3), &
          idp(ipart), levelp(ipart)
  end do

  close(ilun)
end subroutine write_ascii_parts
! ################################################################################
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
