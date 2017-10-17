subroutine mdl_init(r,g,m,p)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p

  ! MPI initialization
#ifndef WITHOUTMPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,g%myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,g%ncpu,ierr)
  g%myid=g%myid+1 ! Careful with this...
  if(g%myid==1)then
     write(*,'(" Launching MPI with nproc = ",I4)')ncpu
  endif
#else
  g%ncpu=1
  g%myid=1
#endif

  ! Register call-back functions

  f(1)%id = MDL_BCAST_PARAMS
  f(1)%proc => broadcast_params
  MDL_INPUT_MAXSIZE=MAX(MDL_INPUT_MAXSIZE,storage_size(r)/32)
  MDL_OUTPUT_MAXSIZE=MAX(MDL_OUTPUT_MAXSIZE,1)
  
end subroutine mdl_init

subroutine mdl_launch(mdl_function_id,target_cpu,cpu_range,input_size,input_array)
  use mdl_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  integer::mdl_function_id
  integer::target_cpu,cpu_range,input_size
  integer,dimension(1:input_size)::input_array

#ifndef WITHOUTMPI





#endif  

end subroutine mdl_launch



