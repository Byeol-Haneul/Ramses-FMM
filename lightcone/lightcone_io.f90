module lightcone_io_module
  use ramses_commons, only: pst_t
  use amr_parameters, only: flen, dp
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none

contains

  subroutine open_lightcone_file(ilun, filename)
    integer :: ilun
    character(LEN=flen), intent(in) :: filename
#ifndef WITHOUTMPI
    integer :: ierr
#endif

#ifndef WITHOUTMPI
    call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, ilun, ierr)
#else
    write(*,*) 'Opening lightcone file without MPI'
#endif
  end subroutine open_lightcone_file

  subroutine close_lightcone_file(ilun)
    integer :: ilun
#ifndef WITHOUTMPI
    integer :: ierr
#endif

#ifndef WITHOUTMPI
    call MPI_FILE_CLOSE(ilun, ierr)
#else
    write(*,*) 'Closing lightcone file without MPI'
#endif
  end subroutine close_lightcone_file

  subroutine write_buffer(ilun, buffer, nbefore, ntotal, nthbuffer)
    use lightcone_buffer_module

    type(lightcone_buffer), intent(in) :: buffer
    integer, intent(in) :: ilun, nbefore, ntotal, nthbuffer
    integer :: idim
#ifndef WITHOUTMPI
    integer(kind=MPI_OFFSET_KIND) :: offset
    integer :: ierr, status(MPI_STATUS_SIZE)
#endif

    if (.not. buffer_is_empty(buffer)) then
#ifndef WITHOUTMPI
      ! Write positions (properties 1, 2, 3)
      do idim = 1, 3
        offset = int(calculate_write_offset(nbefore, ntotal, idim, nthbuffer, buffer%nstride), kind=MPI_OFFSET_KIND)
        call MPI_FILE_WRITE_AT(ilun, offset, buffer%xp(1:buffer%ncurrent, idim), buffer%ncurrent, MPI_REAL, status, ierr)
      end do
      ! Write velocities (properties 4, 5, 6)
      do idim = 1, 3
        offset = int(calculate_write_offset(nbefore, ntotal, idim+3, nthbuffer, buffer%nstride), kind=MPI_OFFSET_KIND)
        call MPI_FILE_WRITE_AT(ilun, offset, buffer%vp(1:buffer%ncurrent, idim), buffer%ncurrent, MPI_REAL, status, ierr)
      end do
#else
      write(*,*) 'Writing buffer to file without MPI'
#endif
    end if
  end subroutine write_buffer

  subroutine write_lightcone_txt_file(filename, ntotal, aexp_old, aexp)
    character(LEN=flen), intent(in) :: filename
    integer, intent(in) :: ntotal
    real(dp), intent(in) :: aexp_old, aexp
    integer :: ilun

    ilun = 1171 ! random number

    open(ilun, file=TRIM(filename)//".txt", form='formatted')
    rewind(ilun)
    write(ilun, *) ntotal
    write(ilun, *) aexp_old
    write(ilun, *) aexp
    close(ilun)
  end subroutine write_lightcone_txt_file

  subroutine output_lightcone_parameters(pst)
        type(pst_t), intent(in) :: pst

        write(*, *) ""
        write(*, *) "Lightcone parameters"
        write(*, *) "  opening_angle_y: ", pst%s%r%cone_opening_angle_y
        write(*, *) "  opening_angle_z: ", pst%s%r%cone_opening_angle_z
        write(*, *) "  cone_x_theta: ", pst%s%r%cone_theta
        write(*, *) "  cone_x_phi: ", pst%s%r%cone_phi
        write(*, *) "  z_max: ", pst%s%r%cone_z_max
        write(*, *) "  z_min: ", pst%s%r%cone_z_min
        write(*, *) "  observer position: ", pst%s%r%cone_observer
        write(*, *) ""
    end subroutine output_lightcone_parameters

  integer function calculate_write_offset(nbefore, ntotal, nthproperty, nthbuffer, nstride)
    ! Calculate the offset at which to write the buffer
    ! nbefore is the number of particles in the preceding processes
    ! ntotal is the total number of particles across all processes
    ! p is the property number (x=1, y=2, z=3, vx=4, vy=5, vz=6, redshift=7)
    ! b is the buffer number (i.e if current process has already written 2 buffers, b=3)
    
    ! The final format is:
    ! All x coordinates, all y coordinates, etc for all properties
    ! So we need to:
    ! - skip all previous properties
    ! - skip all previous processes for the current property
    ! - skip all previous buffers for the current process and property
    integer, intent(in) :: nbefore, ntotal, nthproperty, nthbuffer, nstride
    integer, parameter :: bytes_per_int = 4
    calculate_write_offset = bytes_per_int * ((nthproperty - 1) * ntotal + nbefore + (nthbuffer - 1) * nstride)
  end function calculate_write_offset

end module lightcone_io_module 