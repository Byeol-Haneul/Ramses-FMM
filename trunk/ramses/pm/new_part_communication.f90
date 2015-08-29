module particle_communication
  
contains
#ifndef WITHOUTMPI
  subroutine build_communicator(communicator, recv_tot, ndata, local_data, local_data_oft, &
       keys2, keys1, keys0, ilevel)
    use amr_parameters, only: nlevelmax
    use amr_commons,    only: ncpu, myid, bound_key_level, bound_key
    implicit none

    include 'mpif.h'
    integer, intent(in) ::  ilevel, ndata
    integer, intent(inout) ::  recv_tot, local_data, local_data_oft
    integer(kind=8), dimension(1:ndata), intent(in) :: keys2, keys1, keys0
    ! Wrap send/receive counters/offsets in one array to make the passing them
    ! around a little more convenient
    integer, dimension(1:ncpu, 1:4), intent(inout) :: communicator

    !----------------------------------------------------------------------------
    ! This routine sets up the communication structure for any kind 3-integer 
    ! quantity (particles, bins, etc.). The input hilbert keys (keys) are assumed
    ! to be sorted.
    !----------------------------------------------------------------------------

    integer :: receive_cpu, idata, info, icpu, idest, isource
    integer :: countrecv, countsend  
    integer, dimension(MPI_STATUS_SIZE,2*ncpu) :: statuses
    integer, dimension(2*ncpu)                 :: reqsend, reqrecv

    if (nlevelmax>20)then 
       print*, 'problem here with precision'
       stop
    end if

    communicator = 0
    recv_tot = 0

    ! count number of bins that need to be sent to every other process
    receive_cpu = 1; local_data = 0
    do idata = 1, ndata
       ! REPLACE this by do while( gt_3_keys(keys..., bound_key_level )):
       ! do while (gt_3keys_individual_input(keys(idata,2),keys(idata,1),keys(idata,0), &
       ! bound_key_level(idata,2),bound_key_level(idata,1),bound_key_level(idata,0)))
       do while (keys0(idata) > bound_key_level(receive_cpu, ilevel) ) 
          receive_cpu = receive_cpu + 1
       end do
       communicator(receive_cpu, 1) = communicator(receive_cpu, 1) + 1
    end do

    ! send 1 integer to every other process (including itself - a bit silly, but who cares...)
    ! results in  receive counter
    countrecv = 0; countsend = 0  
    do isource = 1, ncpu        
       countrecv = countrecv + 1
       call MPI_IRECV(communicator(isource, 3), 1, MPI_INTEGER, isource - 1, &
            1234, MPI_COMM_WORLD, reqrecv(countrecv), info)
    end do
    do idest = 1, ncpu
       countsend = countsend + 1
       call MPI_ISEND(communicator(idest, 1), 1, MPI_INTEGER, idest - 1, 1234, &
            MPI_COMM_WORLD, reqsend(countsend), info)
    end do

    call MPI_WAITALL(ncpu, reqrecv, statuses, info)
    call MPI_WAITALL(ncpu, reqsend, statuses, info)

    ! "Prefix sum" to compute send offsets
    do icpu = 1, ncpu - 1
       communicator(icpu + 1, 2) = communicator(icpu, 1) + communicator(icpu, 2)
    end do

    ! No information to is sent to itself, store the number of local data
    local_data = communicator(myid, 1)
    local_data_oft = communicator(myid, 2)
    communicator(myid, 1) = 0
    communicator(myid, 3) = 0

    do icpu = 1, ncpu - 1
       communicator(icpu + 1, 4) = communicator(icpu, 3) + communicator(icpu, 4)
    end do

    recv_tot=sum(communicator(:, 3))

  end subroutine build_communicator
  !################################################################
  !################################################################
  subroutine part_data_to_domain_i4(communicator, send_data, recv_data)
    use amr_commons,   only: ncpu
    implicit none
    include 'mpif.h'

    integer, dimension(1:ncpu, 1:4), intent(in) :: communicator
    integer, dimension(1:), intent(in) :: send_data
    integer, dimension(1:), intent(inout) :: recv_data

    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)


    call MPI_IALLTOALLV(send_data, communicator(:,1), communicator(:,2), MPI_INTEGER, &
         &              recv_data, communicator(:,3), communicator(:,4), MPI_INTEGER, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine part_data_to_domain_i4
  !################################################################
  !################################################################
  subroutine part_data_to_domain_i8(communicator, send_data, recv_data)
    use amr_commons, only: ncpu
    implicit none
    include 'mpif.h'

    integer, dimension(1:ncpu, 1:4), intent(in) :: communicator
    integer(kind = 8), dimension(:), intent(in) :: send_data
    integer(kind = 8), dimension(:), intent(inout) :: recv_data

    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)

    call MPI_IALLTOALLV(send_data, communicator(:,1), communicator(:,2), MPI_INTEGER8, &
         &              recv_data, communicator(:,3), communicator(:,4), MPI_INTEGER8, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine part_data_to_domain_i8
  !################################################################
  !################################################################
  subroutine part_data_to_domain_dp(communicator, send_data, recv_data)
    use amr_commons,   only: ncpu, dp
    implicit none
    include 'mpif.h'

    integer, dimension(1:ncpu, 1:4), intent(in) :: communicator
    real(dp), dimension(:), intent(in) :: send_data
    real(dp), dimension(:), intent(inout) :: recv_data

    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)

    call MPI_IALLTOALLV(send_data, communicator(:,1), communicator(:,2), MPI_DOUBLE_PRECISION, &
         &              recv_data, communicator(:,3), communicator(:,4), MPI_DOUBLE_PRECISION, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine part_data_to_domain_dp
  !################################################################
  !################################################################
  subroutine domain_data_to_part_i4(communicator, send_data, recv_data)
    use amr_commons,   only: ncpu
    implicit none
    include 'mpif.h'

    integer, dimension(1:ncpu, 1:4), intent(in) :: communicator
    integer, dimension(:), intent(in) :: send_data
    integer, dimension(:), intent(inout) :: recv_data

    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)

    call MPI_IALLTOALLV(send_data, communicator(:,3), communicator(:,4), MPI_INTEGER, &
         &              recv_data, communicator(:,1), communicator(:,2), MPI_INTEGER, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine domain_data_to_part_i4
  !################################################################
  !################################################################
  subroutine domain_data_to_part_dp(communicator, send_data, recv_data)
    use amr_commons,   only: ncpu, dp
    implicit none
    include 'mpif.h'

    integer, dimension(1:ncpu, 1:4), intent(in) :: communicator
    real(dp), dimension(:), intent(in) :: send_data
    real(dp), dimension(:), intent(inout) :: recv_data

    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)

    call MPI_IALLTOALLV(send_data, communicator(:,3), communicator(:,4), MPI_DOUBLE_PRECISION, &
         &              recv_data, communicator(:,1), communicator(:,2), MPI_DOUBLE_PRECISION, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine domain_data_to_part_dp
  !################################################################
  !################################################################

  !################################################################
  !################################################################

#endif
end module particle_communication
