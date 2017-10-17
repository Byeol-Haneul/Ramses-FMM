module mdl_parameters

  ! MDL integer arrays max sizes
  integer::MDL_INPUT_MAXSIZE=1
  integer::MDL_OUTPUT_MAXSIZE=1

  ! Parameter for calback functions indices
  integer,parameter::MDL_BCAST_PARAMS=1

end module mdl_parameters

module mdl_commons
  
  use mdl_parameters

  interface
     subroutine ramses_function(r,g,m,p,cpu_range,input_size,output_size,input,output)
       use amr_commons, only: run_t,global_t,mesh_t
       use pm_commons, only: part_t
       use mdl_parameters
       type(run_t)::r
       type(global_t)::g
       type(mesh_t)::m
       type(part_t)::p
       integer::cpu_range,input_size,output_size
       integer,dimension(1:input_size),optional::input
       integer,dimension(1:output_size),optional::output
     end subroutine ramses_function
  end interface

  procedure(ramses_function)::broadcast_params

  type call_back
     integer::id
     procedure(ramses_function),pointer,nopass::proc
  end type call_back

  type(call_back),dimension(1:100)::f

end module mdl_commons
