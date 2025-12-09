module sink_merger_module
#ifndef WITHOUTMPI
  use mpi
#endif
  use constants
  use flag_utils
  use amr_commons
  use mdl_parameters
  use cache_commons
  use amr_parameters, only: ndim, twotondim
  use marshal, only: pack_fetch_flag, unpack_fetch_flag
  use sink_evolution_module, only: sink_B_spline_weights_TSC
  use nbors_utils, only: get_parent_cell

  type :: collision_pair_t
     integer :: id1
     integer :: id2
  end type collision_pair_t

  type :: sink_data_t
     integer :: id
     real(kind=8) :: mass
     real(kind=8), dimension(1:3) :: position
     real(kind=8), dimension(1:3) :: velocity
     logical :: exists
     integer :: cpu_owner
  end type sink_data_t

  type :: out_merger_t
     integer :: dummy_field = 0  ! Add this output type
  end type out_merger_t
contains

  !==============================================================================
  ! MAIN RECURSIVE FUNCTION
  !==============================================================================
  recursive subroutine r_sink_merger(pst,ilevel,input_size)
    use mdl_module
    use ramses_commons, only: pst_t
    use mdl_parameters
    implicit none
    type(pst_t)::pst
    integer,VALUE::input_size
    integer::ilevel

    integer::rID

    if(pst%nLower>0)then
       rID = mdl_send_request(pst%s%mdl,MDL_SINK_MERGER,pst%iUpper+1,input_size,0,ilevel)
       call r_sink_merger(pst%pLower,ilevel,input_size)
       call mdl_get_reply(pst%s%mdl,rID,0)
    else
       call sink_merger_optimized(pst,ilevel)
    endif

  end subroutine r_sink_merger

  !==============================================================================
  ! OPTIMIZED: Single efficient traversal
  !==============================================================================
  subroutine sink_merger_optimized(pst, ilevel)
    use ramses_commons, only: pst_t, ramses_t
    use pm_commons, only: part_t
    use mdl_parameters
    use amr_parameters, only: ndim
    implicit none
    type(pst_t)::pst
    integer::ilevel

    type(collision_pair_t),dimension(:),allocatable::local_collision_pairs
    integer,dimension(:),allocatable::all_id1, all_id2
    integer,dimension(:),allocatable::unique_ids
    type(sink_data_t),dimension(:),allocatable::global_sink_data
    integer::i, j, ierr, nBHnei
    integer::n_count_local, n_count_total, n_valid_mergers, n_unique
    real(kind=8)::dx_loc, factG

    if(pst%s%r%verbose) write(*,*) 'Entering sink merger...' 

#ifndef WITHOUTMPI
    ! Set parameters
    dx_loc = pst%s%r%boxlen/2**ilevel
    factG = 1.0d0
    nBHnei = 3**ndim
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! Step 1: Deposit sink IDs to grid
    call sink_id_deposition(pst%s, pst%s%sink, ilevel, dx_loc, nBHnei)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

    ! Step 2: Collect local collisions
    allocate(local_collision_pairs(1000))
    call collect_collisions(pst%s, pst%s%sink, ilevel, dx_loc, nBHnei, local_collision_pairs, n_count_local)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if(pst%s%r%verbose) write(*,*) 'Proc:', pst%s%g%myid
    if(pst%s%r%verbose) write(*,*) 'Local collision count:', n_count_local


    ! Step 3: Gather all collisions across MPI
    call gather_all_collisions(local_collision_pairs, n_count_local, all_id1, all_id2, n_count_total)
    if(n_count_total == 0) then
       if(pst%s%r%verbose .and. pst%s%g%myid == 1) write(*,*) 'No collisions detected'
       deallocate(local_collision_pairs)
       return
    endif

    if(pst%s%r%verbose .and. pst%s%g%myid == 1) write(*,*) 'Total collision count:', n_count_total

    ! Step 4: Build unique sink ID list
    call build_unique_sink_list(all_id1, all_id2, n_count_total, unique_ids, n_unique)

    if(pst%s%r%verbose .and. pst%s%g%myid == 1) write(*,*) 'Unique sinks involved:', n_unique

    ! Step 5: Gather sink data with MPI
    allocate(global_sink_data(n_unique))
    call mpi_gather_all_sink_data(pst%s%sink, unique_ids, n_unique, global_sink_data)

    ! Step 6: Apply filtering with pre-gathered data
    call filter_collision_pairs_fast(all_id1, all_id2, n_count_total, n_valid_mergers, &
         dx_loc, factG, global_sink_data, n_unique)

    if(pst%s%r%verbose .and. pst%s%g%myid == 1) write(*,*) 'Valid mergers after filtering:', n_valid_mergers

    ! Step 7: Execute mergers
    if(n_valid_mergers > 0) then
       call execute_mergers_batch(pst%s, pst%s%sink, all_id1, all_id2, n_count_total, &
            global_sink_data, n_unique)
    endif

    deallocate(local_collision_pairs, all_id1, all_id2, unique_ids, global_sink_data)

    if(pst%s%r%verbose .and. pst%s%g%myid == 1) write(*,*) 'Sink merger process complete, executed:', n_valid_mergers
#endif

  end subroutine sink_merger_optimized

  !==============================================================================
  ! Sink ID deposition
  !==============================================================================
  subroutine sink_id_deposition(s,p,ilevel,dx_loc,nBHnei) 
    use amr_parameters, only: ndim, twotondim
    use amr_commons, only: oct
    use ramses_commons, only: ramses_t
    use pm_parameters
    use pm_commons, only: part_t
    use nbors_utils
    use cache_commons
    use cache
    use hilbert
    implicit none
    type(ramses_t)::s
    type(part_t)::p
    integer::ilevel,nBHnei
    real(kind=8)::dx_loc

    real(kind=8),dimension(1:ndim)::xcen
    integer(kind=8),dimension(0:ndim)::hash_nbor
    integer::ipart,icelln,j
    type(oct),pointer::gridn
    
    real(kind=8),dimension(1:ndim,1:nBHnei)::xBHnei
    integer,dimension(1:ndim,1:nBHnei)::ckeynei
    real(kind=8),dimension(1:nBHnei)::vol
    real(kind=8)::weight
    integer::ind,igrid,new_id
    type(msg_int4)::dummy_int4

#ifdef HYDRO
#if NDIM==3
    associate(r=>s%r,g=>s%g,m=>s%m)

    ! Set flag1 to max possible index
    do igrid=m%head(ilevel),m%tail(ilevel)
       do ind=1,twotondim
          m%grid(igrid)%flag1(ind)=huge(1)
       end do
    end do

    call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
         hilbert=m%domain,pack_size=storage_size(dummy_int4)/32,&
         pack=pack_fetch_flag, unpack=unpack_fetch_flag,&
         init=init_flush_idsinkmin, flush=pack_flush_idsinkmin, combine=unpack_flush_idsinkmin)

    hash_nbor(0) = ilevel+1
    do ipart = p%headp(r%nlevelmax), p%tailp(r%nlevelmax)

       ! Skip zero-mass sinks
       if(p%mp(ipart) <= 0.0d0) cycle

       ! Compute TSC weights
       do j=1,ndim
          xcen(j) = p%xp(ipart,j)/dx_loc
       end do
       call sink_B_spline_weights_TSC(s,xcen(1:ndim),xBHnei,ckeynei,vol,ilevel)

       do j = 1, nBHnei
          if(vol(j) <= 0.0d0) cycle

          hash_nbor(1:ndim) = ckeynei(1:ndim,j)
          call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.false.)

          if(.not.associated(gridn)) cycle

          weight = vol(j)
          if(weight > 0.0d0) then
             new_id = int(p%idp(ipart), kind=4)
             if(new_id < gridn%flag1(icelln)) then
                gridn%flag1(icelln) = new_id
             endif
          endif
       end do
    end do
    call close_cache(s,m%grid_dict)

    end associate
#endif
#endif
  end subroutine sink_id_deposition

  !==============================================================================
  ! Cache initialization routines
  !==============================================================================
  subroutine init_flush_idsinkmin(grid,hash_key)
    use amr_parameters, only: ndim,twotondim
    use amr_commons, only: oct
    use cache_commons, only: msg_int4
    implicit none
    type(oct)::grid
    integer(kind=8),dimension(0:ndim)::hash_key
    integer::ind
  
    grid%lev=hash_key(0)
    grid%ckey(1:ndim)=hash_key(1:ndim)
    do ind=1,twotondim
       grid%flag1(ind)=huge(1)  ! Initialize all cells to large value
    end do
  end subroutine init_flush_idsinkmin

  subroutine pack_flush_idsinkmin(grid,msg_size,msg_array)
    use amr_parameters, only: twotondim
    use amr_commons, only: oct
    use cache_commons, only: msg_int4
    implicit none
    type(oct)::grid
    integer::msg_size
    integer,dimension(1:msg_size),optional::msg_array
    integer::ind
    type(msg_int4)::msg

    do ind=1,twotondim
       msg%int4(ind)=grid%flag1(ind)
    end do
  
    msg_array=transfer(msg,msg_array)
  end subroutine pack_flush_idsinkmin

  subroutine unpack_flush_idsinkmin(grid,msg_size,msg_array,hash_key)
    use amr_parameters, only: ndim,twotondim
    use amr_commons, only: oct
    use cache_commons, only: msg_int4
    implicit none
    type(oct)::grid
    integer::msg_size
    integer,dimension(1:msg_size),optional::msg_array
    integer(kind=8),dimension(0:ndim)::hash_key
    integer::ind
    type(msg_int4)::msg

    msg=transfer(msg_array,msg)
  
    do ind=1,twotondim
       if(msg%int4(ind) < grid%flag1(ind)) then
          grid%flag1(ind)=msg%int4(ind)  ! Take minimum sink ID
       endif
    end do
  end subroutine unpack_flush_idsinkmin
  
  !==============================================================================
  ! Gather all collisions across MPI - FIXED VERSION
  !==============================================================================
  subroutine gather_all_collisions(local_pairs, n_local, all_id1, all_id2, n_total)
    type(collision_pair_t),dimension(:)::local_pairs
    integer::n_local
    integer,dimension(:),allocatable::all_id1, all_id2
    integer::n_total

    integer::i, ierr, nprocs, myrank
    integer,dimension(:),allocatable::recvcounts, displs
    integer,dimension(:),allocatable::temp_id1, temp_id2  ! Temporary arrays

#ifndef WITHOUTMPI
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    call MPI_ALLREDUCE(n_local, n_total, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    if(n_total == 0) return
  
    allocate(recvcounts(nprocs))
    allocate(displs(nprocs))

    call MPI_ALLGATHER(n_local, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
  
    displs(1) = 0
    do i = 2, nprocs
       displs(i) = displs(i-1) + recvcounts(i-1)
    end do

    ! Allocate temporary arrays for local data
    allocate(temp_id1(n_local))
    allocate(temp_id2(n_local))
  
    ! Copy from derived type to temporary arrays
    do i = 1, n_local
       temp_id1(i) = local_pairs(i)%id1
       temp_id2(i) = local_pairs(i)%id2
    end do

    ! Now use the temporary arrays for MPI communication
    allocate(all_id1(n_total))
    allocate(all_id2(n_total))
  
    call MPI_ALLGATHERV(temp_id1, n_local, MPI_INTEGER, &
         all_id1, recvcounts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHERV(temp_id2, n_local, MPI_INTEGER, &
         all_id2, recvcounts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    deallocate(temp_id1, temp_id2, recvcounts, displs)
#endif
  end subroutine gather_all_collisions

  !==============================================================================
  ! Build unique sink ID list from collision pairs
  !==============================================================================
  subroutine build_unique_sink_list(all_id1, all_id2, n_pairs, unique_ids, n_unique)
    integer,dimension(:)::all_id1, all_id2
    integer::n_pairs
    integer,dimension(:),allocatable::unique_ids
    integer::n_unique
    
    integer::i, j
    logical::found
    integer,dimension(:),allocatable::temp_ids
    
    allocate(temp_ids(2*n_pairs))
    n_unique = 0
    
    do i = 1, n_pairs
       if(all_id1(i) == 0) cycle
       
       found = .false.
       do j = 1, n_unique
          if(temp_ids(j) == all_id1(i)) then
             found = .true.
             exit
          endif
       end do
       if(.not. found) then
          n_unique = n_unique + 1
          write(*,*) 'Sink Collision ID:', all_id1(i)
          temp_ids(n_unique) = all_id1(i)
       endif
       found = .false.
       do j = 1, n_unique
          if(temp_ids(j) == all_id2(i)) then
             found = .true.
             exit
          endif
       end do
       if(.not. found) then
          n_unique = n_unique + 1
          temp_ids(n_unique) = all_id2(i)
       endif
    end do
    
    allocate(unique_ids(n_unique))
    unique_ids(1:n_unique) = temp_ids(1:n_unique)
    deallocate(temp_ids)
    
  end subroutine build_unique_sink_list

  !==============================================================================
  ! MPI-BASED SINK DATA GATHERING USING STANDARD MPI CALLS
  !==============================================================================
  subroutine mpi_gather_all_sink_data(p, unique_ids, n_unique, sink_data)
    use pm_commons, only: part_t
    implicit none
    type(part_t)::p
    integer,dimension(:)::unique_ids
    integer::n_unique
    type(sink_data_t),dimension(:)::sink_data

    integer::i, ipart, j, ierr, myrank, nprocs
    real(kind=8),dimension(:,:),allocatable::all_masses
    real(kind=8),dimension(:,:,:),allocatable::all_positions, all_velocities
    logical,dimension(:,:),allocatable::all_exists
    integer,dimension(:,:),allocatable::all_owners

#ifndef WITHOUTMPI
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

    ! Allocate arrays for gathering data from all processes
    allocate(all_masses(n_unique, nprocs))
    allocate(all_positions(n_unique, 3, nprocs))  
    allocate(all_velocities(n_unique, 3, nprocs))
    allocate(all_exists(n_unique, nprocs))
    allocate(all_owners(n_unique, nprocs))

    ! Initialize local data sections
    all_masses(:, myrank+1) = 0.0d0
    all_exists(:, myrank+1) = .false.
    all_owners(:, myrank+1) = -1
    all_positions(:, :, myrank+1) = 0.0d0
    all_velocities(:, :, myrank+1) = 0.0d0

    ! Fill local sink data
    do ipart = 1, p%npart
       do i = 1, n_unique
          if(unique_ids(i) == p%idp(ipart)) then
             all_exists(i, myrank+1) = .true.
             all_masses(i, myrank+1) = p%mp(ipart)
             all_positions(i, 1:3, myrank+1) = p%xp(ipart, 1:3)
             all_velocities(i, 1:3, myrank+1) = p%vp(ipart, 1:3)
             all_owners(i, myrank+1) = myrank
             exit
          endif
       end do
    end do

    ! Gather all data from all processes
    call MPI_ALLGATHER(MPI_IN_PLACE, n_unique, MPI_DOUBLE_PRECISION, &
         all_masses, n_unique, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER(MPI_IN_PLACE, n_unique, MPI_LOGICAL, &
         all_exists, n_unique, MPI_LOGICAL, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER(MPI_IN_PLACE, n_unique, MPI_INTEGER, &
         all_owners, n_unique, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER(MPI_IN_PLACE, 3*n_unique, MPI_DOUBLE_PRECISION, &
         all_positions, 3*n_unique, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER(MPI_IN_PLACE, 3*n_unique, MPI_DOUBLE_PRECISION, &
         all_velocities, 3*n_unique, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

    ! Combine data - each process builds the complete sink dataset
    do i = 1, n_unique
       sink_data(i)%id = unique_ids(i)
       sink_data(i)%exists = .false.
       sink_data(i)%mass = 0.0d0
       sink_data(i)%position = 0.0d0
       sink_data(i)%velocity = 0.0d0
       sink_data(i)%cpu_owner = -1

       ! Find first valid entry across all processes
       do j = 1, nprocs
          if(all_exists(i, j)) then
             sink_data(i)%exists = .true.
             sink_data(i)%mass = all_masses(i, j)
             sink_data(i)%position(1:3) = all_positions(i, 1:3, j)
             sink_data(i)%velocity(1:3) = all_velocities(i, 1:3, j)
             sink_data(i)%cpu_owner = all_owners(i, j)
             exit
          endif
       end do
    end do

    deallocate(all_masses, all_positions, all_velocities, all_exists, all_owners)
#endif

  end subroutine mpi_gather_all_sink_data

  !==============================================================================
  ! FAST FILTERING with pre-gathered sink data
  !==============================================================================
  subroutine filter_collision_pairs_fast(all_id1, all_id2, n_total, n_filtered, &
                                        dx_loc, factG, sink_data, n_sinks)
    integer,dimension(:)::all_id1, all_id2
    integer::n_total, n_filtered, n_sinks
    real(kind=8)::dx_loc, factG
    type(sink_data_t),dimension(:)::sink_data

    integer::i, j, sink1_idx, sink2_idx
    logical::should_merge

    ! Sort collision pairs
    call sort_collision_pairs(all_id1, all_id2, n_total)

    ! Remove exact duplicates
    do i = 2, n_total
      if(all_id1(i) == 0) cycle
      if(all_id1(i) == all_id1(i-1) .and. all_id2(i) == all_id2(i-1)) then
        all_id1(i) = 0
        all_id2(i) = 0
      endif
    end do

    ! Apply merger criteria using pre-gathered data
    do i = 1, n_total
      if(all_id1(i) == 0) cycle

      ! Find sink data indices
      sink1_idx = 0
      sink2_idx = 0
      do j = 1, n_sinks
         if(sink_data(j)%id == all_id1(i)) sink1_idx = j
         if(sink_data(j)%id == all_id2(i)) sink2_idx = j
      end do

      ! Check if both sinks exist and have mass
      if(sink1_idx == 0 .or. sink2_idx == 0 .or. &
         .not. sink_data(sink1_idx)%exists .or. .not. sink_data(sink2_idx)%exists .or. &
         sink_data(sink1_idx)%mass <= 0.0d0 .or. sink_data(sink2_idx)%mass <= 0.0d0) then
        all_id1(i) = 0
        all_id2(i) = 0
        cycle
      endif

      ! Apply physical criteria
      should_merge = check_merger_criteria(sink_data(sink1_idx), sink_data(sink2_idx), dx_loc, factG)
      if(.not. should_merge) then
        all_id1(i) = 0
        all_id2(i) = 0
      endif
    end do

    ! Remove triple mergers and duplicates
    do i = 1, n_total
      if(all_id2(i) == 0) cycle
      do j = 1, n_total
        if(all_id1(j) == 0) cycle
        if(i == j) cycle
        if(all_id2(i) == all_id1(j)) then
          all_id1(j) = 0
          all_id2(j) = 0
          exit
        endif
      end do
    end do

    do i = 1, n_total
      if(all_id1(i) == 0) cycle
      do j = i+1, n_total
        if(all_id1(j) == 0) cycle
        if(all_id1(i) == all_id1(j)) then
          all_id1(j) = 0
          all_id2(j) = 0
        endif
      end do
    end do

    ! Count filtered pairs
    n_filtered = 0
    do i = 1, n_total
      if(all_id1(i) /= 0) n_filtered = n_filtered + 1
    end do

  end subroutine filter_collision_pairs_fast

  !==============================================================================
  ! FAST BATCH MERGER EXECUTION - MPI AWARE
  !==============================================================================
  subroutine execute_mergers_batch(s, p, all_id1, all_id2, n_total, sink_data, n_sinks)
    use ramses_commons, only: ramses_t
    use pm_commons, only: part_t
    implicit none
    type(ramses_t)::s
    type(part_t)::p
    integer,dimension(:)::all_id1, all_id2
    integer::n_total, n_sinks
    type(sink_data_t),dimension(:)::sink_data

    integer::i, j, myrank, ierr
    integer::id_keep, id_delete
    integer::sink1_idx, sink2_idx
    real(kind=8)::mass1, mass2, total_mass
    real(kind=8),dimension(1:3)::com_position, momentum

#ifndef WITHOUTMPI
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    do i = 1, n_total
       id_keep = all_id1(i)
       id_delete = all_id2(i)

       ! Skip invalid pairs
       if(id_keep == 0) cycle

       ! Find the sink data for this pair
       sink1_idx = 0
       sink2_idx = 0
       do j = 1, n_sinks
          if(sink_data(j)%id == id_keep) sink1_idx = j
          if(sink_data(j)%id == id_delete) sink2_idx = j
       end do

       ! Skip if we couldn't find both sinks
       if(sink1_idx == 0 .or. sink2_idx == 0) cycle

       ! Calculate merged properties (ALL CPUs do this for consistency)
       mass1 = sink_data(sink1_idx)%mass
       mass2 = sink_data(sink2_idx)%mass
       total_mass = mass1 + mass2

       com_position(1:3) = (mass1 * sink_data(sink1_idx)%position(1:3) + &
            &               mass2 * sink_data(sink2_idx)%position(1:3)) / total_mass
       momentum(1:3) = mass1 * sink_data(sink1_idx)%velocity(1:3) + &
            &          mass2 * sink_data(sink2_idx)%velocity(1:3)

       ! Each CPU only updates the sinks it owns
       if(sink_data(sink1_idx)%cpu_owner == myrank) then
          ! Find and update the local sink we're keeping
          do j = 1, p%npart
             if(p%idp(j) == id_keep) then
                p%mp(j) = total_mass
                p%xp(j,1:3) = com_position(1:3)
                p%vp(j,1:3) = momentum(1:3) / total_mass
                if(s%r%verbose) write(*,*) 'CPU', myrank, ': Updated sink', id_keep, ' mass=', total_mass
                exit
             endif
          end do
       endif

       if(sink_data(sink2_idx)%cpu_owner == myrank) then
          ! Find and delete the local sink we're removing
          do j = 1, p%npart
             if(p%idp(j) == id_delete) then
                p%mp(j) = 0.0d0  ! Mark for deletion
                p%xp(j,1:3) = com_position(1:3)
                p%vp(j,1:3) = momentum(1:3) / total_mass
                if(s%r%verbose) write(*,*) 'CPU', myrank, ': Deleted sink', id_delete
                exit
             endif
          end do
       endif
    end do
#endif
  end subroutine execute_mergers_batch

  !==============================================================================
  ! Collect collision pairs
  !==============================================================================
  subroutine collect_collisions(s,p,ilevel,dx_loc,nBHnei,all_pairs,n_local)
    use amr_parameters, only: ndim, twotondim
    use ramses_commons, only: ramses_t
    use pm_parameters
    use pm_commons, only: part_t
    use cache_commons
    use cache
    implicit none
    type(ramses_t)::s
    type(part_t)::p
    integer::ilevel,nBHnei,n_local
    real(kind=8)::dx_loc
    type(collision_pair_t),dimension(:)::all_pairs

    real(kind=8),dimension(1:ndim)::xcen
    integer(kind=8),dimension(0:ndim)::hash_nbor
    integer::ipart,icelln,j,my_id,neighbor_id
    type(oct),pointer::gridn

    real(kind=8),dimension(1:ndim,1:nBHnei)::xBHnei
    integer,dimension(1:ndim,1:nBHnei)::ckeynei
    real(kind=8),dimension(1:nBHnei)::vol
    type(msg_int4)::dummy_int4

#ifdef HYDRO
#if NDIM==3
    associate(r=>s%r,g=>s%g,m=>s%m)

    if(s%r%accretion_type==0)return

    if(r%verbose .and. g%myid == 1)write(*,*)'Collecting collision pairs (using TSC)...'

    n_local = 0

    hash_nbor(0) = ilevel+1
    call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
         pack=pack_fetch_flag,unpack=unpack_fetch_flag,&
         hilbert=m%domain,pack_size=storage_size(dummy_int4)/32)

    do ipart = p%headp(r%nlevelmax), p%tailp(r%nlevelmax)
       my_id = p%idp(ipart)

       ! Skip zero-mass sinks
       if(p%mp(ipart) <= 0.0d0) cycle

       do j=1,ndim
          xcen(j) = p%xp(ipart,j)/dx_loc
       end do

       call sink_B_spline_weights_TSC(s,xcen(1:ndim),xBHnei,ckeynei,vol,ilevel)

       do j = 1, nBHnei
          if(vol(j) <= 0.0d0) cycle

          hash_nbor(1:ndim) = ckeynei(1:ndim,j)
          call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.false.,fetch_cache=.true.)

          if(associated(gridn)) then
             neighbor_id = gridn%flag1(icelln)

             if(neighbor_id < my_id) then
                n_local = n_local + 1
                if(n_local <= size(all_pairs)) then
                   write(*,*)'NEIGHBOR ID', neighbor_id
                   write(*,*)'MY ID', my_id
                   all_pairs(n_local)%id1 = neighbor_id
                   all_pairs(n_local)%id2 = my_id
                endif
             endif
          endif
       end do
    end do

    call close_cache(s,m%grid_dict)

    if(r%verbose .and. g%myid == 1)write(*,*)'Collected',n_local,'collision pairs'

    end associate
#endif
#endif
  end subroutine collect_collisions

  !==============================================================================
  ! Check merger criteria (binding energy)
  !==============================================================================
  logical function check_merger_criteria(sink1,sink2,dx_loc,factG)
    implicit none
    type(sink_data_t)::sink1,sink2
    real(kind=8)::dx_loc,factG

    real(kind=8)::dx,dy,dz,separation,dvx,dvy,dvz,vel_squared
    real(kind=8)::total_mass,binding_criteria

    check_merger_criteria = .false.

    dx = sink1%position(1) - sink2%position(1)
    dy = sink1%position(2) - sink2%position(2)
    dz = sink1%position(3) - sink2%position(3)
    separation = sqrt(dx*dx + dy*dy + dz*dz)

    dvx = sink1%velocity(1) - sink2%velocity(1)
    dvy = sink1%velocity(2) - sink2%velocity(2)
    dvz = sink1%velocity(3) - sink2%velocity(3)
    vel_squared = dvx*dvx + dvy*dvy + dvz*dvz

    total_mass = sink1%mass + sink2%mass

    if(separation < dx_loc) then
       binding_criteria = (factG * total_mass) / dx_loc * (1.0d0 - (separation/dx_loc)**2)

       if(vel_squared < binding_criteria) then
          check_merger_criteria = .true.
       endif
    endif

  end function check_merger_criteria

  !==============================================================================
  ! Sort collision pairs by id1, then id2 (bubble sort)
  !==============================================================================
  subroutine sort_collision_pairs(all_id1, all_id2, n)
    implicit none
    integer,dimension(:)::all_id1, all_id2
    integer::n
    
    integer::i, j
    integer::temp_id1, temp_id2
    logical::swapped
    
    do i = n, 1, -1
      swapped = .false.
      do j = 1, i-1
        if(all_id1(j) == 0 .or. all_id1(j+1) == 0) cycle
        
        if(all_id1(j) > all_id1(j+1) .or. &
           (all_id1(j) == all_id1(j+1) .and. all_id2(j) > all_id2(j+1))) then
          
          temp_id1 = all_id1(j)
          temp_id2 = all_id2(j)
          
          all_id1(j) = all_id1(j+1)
          all_id2(j) = all_id2(j+1)
          
          all_id1(j+1) = temp_id1
          all_id2(j+1) = temp_id2
          
          swapped = .true.
        endif
      end do
      
      if(.not.swapped) exit
    end do
    
  end subroutine sort_collision_pairs

end module sink_merger_module
