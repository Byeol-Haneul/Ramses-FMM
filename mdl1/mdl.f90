module mdl_module

  use mdl_parameters
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_FUNPTR, C_PTR

  type :: mdl_t
     
     integer::myid
     integer::ncpu

     integer::MDL_INPUT_MAXSIZE=1
     integer,dimension(:),allocatable::mpi_input_buffer
     
     ! Communication-related objects
     integer::mail_counter
     integer::request_id,flush_id
     integer,dimension(:),allocatable::reply_id

     ! Adopted combiner rule
     integer::combiner_rule

     ! Message sizes
     integer::size_msg_array
     integer::size_request_array
     integer::size_flush_array
     integer::size_fetch_array

     ! Message arrays
     integer(kind=4),dimension(:),allocatable::recv_request_array
     integer(kind=4),dimension(:),allocatable::send_request_array
     integer(kind=4),dimension(:),allocatable::recv_fetch_array
     integer(kind=4),dimension(:),allocatable::send_fetch_array
     integer(kind=4),dimension(:),allocatable::recv_flush_array
     integer(kind=4),dimension(:),allocatable::send_flush_array

     ! Callback functions
     type(c_funptr),dimension(0:100)::callback
     type(c_ptr),dimension(0:100)::p1opaque
     
  end type mdl_t

  contains

    subroutine mdl_abort(mdl)
      type(mdl_t)::mdl
#ifndef WITHOUTMPI
      include 'mpif.h'
      integer::info
      call MPI_ABORT(MPI_COMM_WORLD,info)
#else
      stop
#endif  
    end subroutine mdl_abort

    subroutine mdl_add_service(mdl,sid,p1,service,nin,nout)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_FUNPTR, C_LOC
      type(mdl_t)::mdl
      integer::sid
      type(*),target::p1
      type(c_funptr), intent(in), value :: service
      integer :: nin, nout ! NOTE: nin is the size IN BYTES!!
      mdl%callback(sid) = service
      mdl%p1opaque(sid) = c_loc(p1)
      mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,nin/4) ! Divide by four to get the number of Integers

    end subroutine mdl_add_service

end module mdl_module