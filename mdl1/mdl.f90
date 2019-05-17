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
     integer(kind=4),dimension(0:100)::input_size, output_size
  end type mdl_t

  interface mdl_send_request
    module procedure mdl_send_request_array, mdl_send_request_scalar
  end interface mdl_send_request

  interface mdl_get_reply
    module procedure mdl_get_reply_array, mdl_get_reply_scalar
  end interface mdl_get_reply

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

    subroutine mdl_add_service(mdl,sid,p1,service,input_size,output_size)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_FUNPTR, C_LOC
      type(mdl_t)::mdl
      integer::sid
      type(*),target::p1
      type(c_funptr), intent(in), value :: service
      integer :: input_size, output_size ! NOTE: size IN BYTES!!
      mdl%callback(sid) = service
      mdl%p1opaque(sid) = c_loc(p1)
      mdl%input_size(sid) = input_size
      mdl%output_size(sid) = output_size
      mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,input_size/4) ! Divide by four to get the number of Integers

    end subroutine mdl_add_service


    subroutine mdl_req_service(mdl,target_cpu,mdl_function_id,input,input_length)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_LOC, C_F_POINTER
      implicit none
      type(mdl_t)::mdl
#ifndef WITHOUTMPI
      include 'mpif.h'
      integer::info
#endif
      integer::mdl_function_id
      integer,value,intent(in)::target_cpu
      integer,value,intent(in),optional::input_length
      type(*),intent(in),optional,target::input

#ifndef WITHOUTMPI
      integer::launch_id,launch_tag=101
      integer,dimension(MPI_STATUS_SIZE)::launch_status
      integer,dimension(1:32)::header=0
      byte,dimension(:),pointer::dummy

      write(*,*) 'NEW:',mdl_function_id

      ! Assemble MPI message
      header(1)=mdl_function_id
      header(2)=input_length/4
!      header(3)=output_size
      header(3)=1
      mdl%mpi_input_buffer(1:32)=header
      mdl%mpi_input_buffer(33) = 99
      if (present(input_length) .and. input_length>0) then
        call c_f_pointer(c_loc(input),dummy,[input_length])
        mdl%mpi_input_buffer(33:32+input_length/4) = transfer(dummy,mdl%mpi_input_buffer(33:32+input_length/4))
      endif

      ! Send input array to the target cpu
      call MPI_ISEND(mdl%mpi_input_buffer,input_length/4+32,MPI_INTEGER,target_cpu-1,launch_tag,MPI_COMM_WORLD,launch_id,info)

      ! Wait for ISEND completion to free memory in corresponding MPI buffer
      call MPI_WAIT(launch_id,launch_status,info)      
#endif
    end subroutine mdl_req_service

    subroutine mdl_send_request_array(mdl,mdl_function_id,target_cpu,input_size,output_size,input_array)
      implicit none
      type(mdl_t)::mdl
      integer,intent(in)::mdl_function_id
      integer,intent(in)::target_cpu,input_size,output_size
      integer,intent(in),dimension(1:input_size)::input_array

#ifndef WITHOUTMPI
      include 'mpif.h'
      integer::info
      integer::launch_id,launch_tag=101
      integer,dimension(MPI_STATUS_SIZE)::launch_status
      integer,dimension(1:32)::header=0

      ! Assemble MPI message
      header(1)=mdl_function_id
      header(2)=input_size
      header(3)=output_size
      mdl%mpi_input_buffer(1:32)=header
      if(input_size>0)then
        mdl%mpi_input_buffer(33:32+input_size)=input_array
      endif

      ! Send input array to the target cpu
      call MPI_ISEND(mdl%mpi_input_buffer,input_size+32,MPI_INTEGER,target_cpu-1,launch_tag,MPI_COMM_WORLD,launch_id,info)

      ! Wait for ISEND completion to free memory in corresponding MPI buffer
      call MPI_WAIT(launch_id,launch_status,info)
#endif  
    end subroutine mdl_send_request_array

    subroutine mdl_send_request_scalar(mdl,mdl_function_id,target_cpu,input_length,output_size,input)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_LOC, C_F_POINTER
      implicit none
      type(mdl_t)::mdl
      integer,intent(in)::mdl_function_id
      integer,intent(in)::target_cpu
      integer,intent(in),optional::input_length,output_size
      type(*),intent(in),optional,target::input

#ifndef WITHOUTMPI
      include 'mpif.h'
      integer::info
      integer::launch_id,launch_tag=101
      integer,dimension(MPI_STATUS_SIZE)::launch_status
      integer,dimension(1:32)::header=0
      byte,dimension(:),pointer::dummy

      write(*,*) 'NEW:',mdl_function_id

      ! Assemble MPI message
      header(1)=mdl_function_id
      header(2)=input_length/4
      header(3)=output_size
      mdl%mpi_input_buffer(1:32)=header
      mdl%mpi_input_buffer(33) = 99
      if (present(input_length) .and. input_length>0) then
        call c_f_pointer(c_loc(input),dummy,[input_length])
        mdl%mpi_input_buffer(33:32+input_length/4) = transfer(dummy,mdl%mpi_input_buffer(33:32+input_length/4))
      endif

      ! Send input array to the target cpu
      call MPI_ISEND(mdl%mpi_input_buffer,input_length/4+32,MPI_INTEGER,target_cpu-1,launch_tag,MPI_COMM_WORLD,launch_id,info)

      ! Wait for ISEND completion to free memory in corresponding MPI buffer
      call MPI_WAIT(launch_id,launch_status,info)      
#endif
    end subroutine mdl_send_request_scalar

    !##############################################################
    !##############################################################
    !##############################################################
    !##############################################################
    subroutine mdl_get_reply_array(mdl,target_cpu,output_size,output_array)
      implicit none
      type(mdl_t)::mdl
      integer::target_cpu
      integer::output_size
      integer,dimension(1:output_size)::output_array

#ifndef WITHOUTMPI  
      include 'mpif.h'
      integer::info
      integer::output_tag=203,output_id,dummy=1
      integer,dimension(MPI_STATUS_SIZE)::output_status  
      
      ! Post a RECV for the output back from target_cpu
      if(output_size>0)then
        call MPI_IRECV(output_array,output_size,MPI_INTEGER,target_cpu-1,output_tag,MPI_COMM_WORLD,output_id,info)
      else
        call MPI_IRECV(dummy,1,MPI_INTEGER,target_cpu-1,output_tag,MPI_COMM_WORLD,output_id,info)
      endif

      ! Wait for ISEND completion to free memory in corresponding MPI buffer
      call MPI_WAIT(output_id,output_status,info)
#endif  
    end subroutine mdl_get_reply_array

    subroutine mdl_get_reply_scalar(mdl,target_cpu,output,output_length)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_LOC, C_F_POINTER
      implicit none
      type(mdl_t)::mdl
#ifndef WITHOUTMPI
      include 'mpif.h'
      integer::info
#endif
      integer,value,intent(in)::target_cpu
      integer,optional::output_length
      type(*),optional,target::output
      byte,dimension(:),pointer::buffer

#ifndef WITHOUTMPI  
      integer::output_tag=203,output_id,dummy=1
      integer,dimension(MPI_STATUS_SIZE)::output_status  
      
      ! Post a RECV for the output back from target_cpu
      if(present(output)) then
        call c_f_pointer(c_loc(output),buffer,[mdl%MDL_INPUT_MAXSIZE])
        write(*,*) 'Post MPI_IRECV',buffer(4)
        call MPI_IRECV(buffer,output_length,MPI_INTEGER,target_cpu-1,output_tag,MPI_COMM_WORLD,output_id,info)
        write(*,*) 'what now'
      else
        call MPI_IRECV(dummy,1,MPI_INTEGER,target_cpu-1,output_tag,MPI_COMM_WORLD,output_id,info)
      endif

!      if(output_size>0)then
!        call MPI_IRECV(output_array,output_size,MPI_INTEGER,target_cpu-1,output_tag,MPI_COMM_WORLD,output_id,info)
!      else
!        call MPI_IRECV(dummy,1,MPI_INTEGER,target_cpu-1,output_tag,MPI_COMM_WORLD,output_id,info)
!      endif

      ! Wait for ISEND completion to free memory in corresponding MPI buffer
      write(*,*) 'Waiting for response'
      call MPI_WAIT(output_id,output_status,info)
      if (present(output_length)) then
        call MPI_GET_COUNT(output_status, MPI_INTEGER, output_length, info)
        write(*,*) 'Length is',output_length
      endif


      write(*,*) 'Got response'

#endif  

    end subroutine mdl_get_reply_scalar




end module mdl_module