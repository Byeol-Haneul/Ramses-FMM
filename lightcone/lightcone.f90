module lightcone_module
  use amr_parameters, only: flen
  use ramses_commons, only: pst_t, ramses_t
  implicit none
contains
  subroutine m_output_lightcone(pst)  
    use mdl_module, only: mdl_mkdir

    type(pst_t), intent(in) :: pst    
    character(LEN=5) :: nchar
    integer :: dummy(1)
    integer, dimension(1:flen/4) :: input_array
    character(LEN=flen) :: filedir
    real(kind=8) :: z1, z2

    associate(r=>pst%s%r, g=>pst%s%g, mdl=>pst%s%mdl)

    z2 = 1 / g%aexp_old - 1.0d0
    z1 = 1 / g%aexp - 1.0d0

    if (g%nstep_coarse < 2 .or. z1 > r%cone_z_max .or. z2 < r%cone_z_min .or. abs(z2 - z1) < 1d-6) return
    if (r%verbose) write(*,*) 'Entering output_lightcone, nstep_coarse: ', g%nstep_coarse

    call title(g%nstep_coarse, nchar)
    filedir='lightcone/'
    call mdl_mkdir(mdl, filedir)

    ! Pass step number instead of full filename
    input_array = transfer(nchar, input_array)
    call r_output_lightcone(pst, input_array, flen/4, dummy, 0)

    end associate
  end subroutine m_output_lightcone

  recursive subroutine r_output_lightcone(pst, input_array, input_size, output_array, output_size)
    use mdl_module, only: mdl_send_request, mdl_get_reply
    use mdl_parameters, only: MDL_OUTPUT_LIGHTCONE

    type(pst_t), intent(in) :: pst
    integer :: rID
    integer, value :: input_size
    integer :: output_size
    integer, dimension(1:input_size) :: input_array
    integer, dimension(1:output_size) :: output_array
    character(LEN=flen) :: part_filename, tree_filename
    character(LEN=5) :: nchar

    if (pst%nLower > 0) then
       rID = mdl_send_request(pst%s%mdl, MDL_OUTPUT_LIGHTCONE, pst%iUpper+1, input_size, output_size, input_array)
       call r_output_lightcone(pst%pLower, input_array, input_size, output_array, output_size)
       call mdl_get_reply(pst%s%mdl, rID, output_size)
    else
       call title(pst%s%g%myid, nchar)
       nchar = TRIM(transfer(input_array, nchar))
       
       ! Build filenames for both particle types
       part_filename = 'lightcone/part_'//TRIM(nchar)
       tree_filename = 'lightcone/tree_'//TRIM(nchar)
       
       ! Output regular DM particles
       call output_lightcone(pst%s, pst%s%p, part_filename)
       
       ! Output tree particles
       if (pst%s%r%tree) then
          call output_lightcone(pst%s, pst%s%tree, tree_filename)
       end if
    endif

  end subroutine r_output_lightcone

  subroutine output_lightcone(s, p, filename)
    use pm_commons, only: part_t
    use lightcone_utils
    use lightcone_buffer_module
    use lightcone_io_module
#ifndef WITHOUTMPI
    use mpi
#endif
    ! Inputs
    type(ramses_t), intent(in) :: s
    type(part_t), intent(in) :: p
    character(LEN=flen), intent(in) :: filename

    ! Fixed parameters
    integer, parameter :: nstride = 65536, max_replicas = 10000

    ! Local variables
    real(kind=8) :: r_inner, r_outer, angle_y, angle_z
    integer :: nreplicas, &
         first_xreplica, last_xreplica, &
         first_yreplica, last_yreplica, &
         first_zreplica, last_zreplica
    real(kind=8) :: z1, z2, coverH0, omega_r
    real(kind=8) :: position(3)
    real(kind=8) :: cone_to_box_rotation(3,3), box_to_cone_rotation(3,3)
    logical, allocatable :: has_particles(:,:,:)
    type(lightcone_buffer) :: buffer
    integer :: npart, nselected, nbefore, ntotal, nthbuffer
    integer :: i, j, k, ilun
    real(kind=8) :: velocity(3)
#ifndef WITHOUTMPI
    integer :: ierr
#endif

    z2 = 1/s%g%aexp_old - 1.0d0
    z1 = 1/s%g%aexp - 1.0d0
    coverH0 = 2.9979246d+5/s%g%h0
    omega_r = 1 - s%g%omega_m - s%g%omega_l

    cone_to_box_rotation = rotation_matrix(deg2rad(s%r%cone_theta), deg2rad(s%r%cone_phi))
    box_to_cone_rotation = transpose(cone_to_box_rotation)

    angle_y = deg2rad(s%r%cone_opening_angle_y)
    angle_z = deg2rad(s%r%cone_opening_angle_z)

    ! Find the relevant replicas
    r_inner = comoving2code(s%g, comoving_distance(z1, s%g%omega_m, s%g%omega_l, omega_r, coverH0))
    r_outer = comoving2code(s%g, comoving_distance(z2, s%g%omega_m, s%g%omega_l, omega_r, coverH0))

    call compute_replica_range(cone_to_box_rotation, s%r%cone_observer, angle_y, angle_z, r_inner, r_outer, first_xreplica, last_xreplica, first_yreplica, last_yreplica, first_zreplica, last_zreplica)

    nreplicas = (last_xreplica - first_xreplica + 1) * (last_yreplica - first_yreplica + 1) * (last_zreplica - first_zreplica + 1)

    ! Check that total number of replicas does not exceed max_replicas
    if (nreplicas > max_replicas) stop 'Number of replicas exceeds max_replicas in lightcone'

    allocate(has_particles(first_xreplica:last_xreplica, first_yreplica:last_yreplica, first_zreplica:last_zreplica))
    has_particles = .false.

    ! Count selection particles and fill has_particles array
    nselected = 0
    do i = first_xreplica, last_xreplica
      do j = first_yreplica, last_yreplica
        do k = first_zreplica, last_zreplica
          do npart = 1, p%npart

            position = p%xp(npart, :)
            position(1) = position(1) + i * s%r%boxlen
            position(2) = position(2) + j * s%r%boxlen
            position(3) = position(3) + k * s%r%boxlen

            position = box_to_cone_coordinates(box_to_cone_rotation, s%r%cone_observer, position)

            if (is_in_lightcone_sector(position, r_inner, r_outer, angle_y, angle_z)) then
              nselected = nselected + 1
              has_particles(i,j,k) = .true.
            end if

          end do
        end do
      end do
    end do

    ! Selection and writing starts here
    ! Process i needs to know the number of particles in the preceding processes
    ! e.g process 5 needs to know N_1 + N_2 + N_3 + N_4 so that it can start writing at the correct position
    nbefore = 0 ! Number of particles in the preceding processes
    ntotal = nselected ! Total number of particles across all processes
#ifndef WITHOUTMPI
    ! Each process gets the total number of particles in the preceding processes
    call MPI_EXSCAN(nselected, nbefore, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    ntotal = nbefore + nselected
    ! Last process broadcasts the total number of particles (its nbefore + nselected) to all other processes
    call MPI_BCAST(ntotal, 1, MPI_INTEGER, s%g%ncpu-1, MPI_COMM_WORLD, ierr)
#endif

    if (s%r%verbose .and. s%g%myid == 1) write(*, *) 'Found ', ntotal, ' lightcone particles across ', nreplicas, ' replicas'

    call init_lightcone_buffer(buffer, nstride) ! Allocate the buffer

    ilun = 3 * s%g%ncpu + s%g%myid + 103
    call open_lightcone_file(ilun, filename)

    ! Select particles and write in chunks
    nthbuffer = 0 ! We need to keep track of the number of buffers previously written to correctly offset the write position
    do i = first_xreplica, last_xreplica
       do j = first_yreplica, last_yreplica
          do k = first_zreplica, last_zreplica

             if (.not. has_particles(i, j, k)) cycle

             do npart = 1, p%npart
                position = p%xp(npart, :)
                position(1) = position(1) + i * s%r%boxlen
                position(2) = position(2) + j * s%r%boxlen
                position(3) = position(3) + k * s%r%boxlen

                position = box_to_cone_coordinates(box_to_cone_rotation, s%r%cone_observer, position)

                if (is_in_lightcone_sector(position, r_inner, r_outer, angle_y, angle_z)) then
                   ! Transform velocity to cone coordinates
                   velocity = matmul(box_to_cone_rotation, p%vp(npart, :))
                   call add_to_buffer(buffer, p%idp(npart), real(position(:), kind=4), real(velocity(:), kind=4), real(p%mp(npart), kind=4))

                   if (buffer_is_full(buffer)) then
                      nthbuffer = nthbuffer + 1
                      call write_buffer(ilun, buffer, nbefore, ntotal, nthbuffer)
                      call empty_buffer(buffer)
                   end if

                end if
             end do
          end do
       end do
    end do

    ! Write any remaining particles
    if (.not. buffer_is_empty(buffer)) then
       nthbuffer = nthbuffer + 1
       call write_buffer(ilun, buffer, nbefore, ntotal, nthbuffer)
       call empty_buffer(buffer)
    end if

    call close_lightcone_file(ilun)
    if (s%g%myid == 1) call write_lightcone_txt_file(filename, ntotal, s%g%aexp_old, s%g%aexp)

  end subroutine output_lightcone

end module lightcone_module
