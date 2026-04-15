! One MPI_File_write_at wrapper per compilation unit for gfortran + mpif.h
! (avoids implicit interface clashes between INTEGER and REAL buffers).
module lightcone_mpi_file_r4_module
  implicit none
contains
  subroutine lc_mpi_file_write_at_r4(fh, offset, buf, count, ierr)
    use amr_parameters, only: sp
    implicit none
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    integer, intent(in) :: fh
    integer(8), intent(in) :: offset
    real(kind=sp), intent(in) :: buf(:)
    integer, intent(in) :: count
    integer, intent(out) :: ierr
#ifndef WITHOUTMPI
    integer :: status(MPI_STATUS_SIZE)

    call MPI_FILE_WRITE_AT(fh, int(offset, kind=MPI_OFFSET_KIND), buf(1:count), count, MPI_REAL, status, ierr)
#else
    ierr = 0
#endif
  end subroutine lc_mpi_file_write_at_r4
end module lightcone_mpi_file_r4_module
