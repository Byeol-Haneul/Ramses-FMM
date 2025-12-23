module call_back

  interface
     recursive subroutine ramses_function(pst,input,input_size,output,output_size)
       use ramses_commons, only: pst_t
       type(pst_t)::pst
       integer,VALUE::input_size
       integer,optional::output_size
       TYPE(*)::input
       TYPE(*),optional::output
     end subroutine ramses_function
  end interface
    
  type call_back_f
     procedure(ramses_function),pointer,nopass::proc
  end type call_back_f

  interface

     subroutine cache_function(m,igrid,msg_size,msg_array)
       use amr_commons, only: mesh_t
       type(mesh_t)::m
       integer::igrid
       integer::msg_size
       integer,dimension(1:msg_size),optional::msg_array
     end subroutine cache_function
     subroutine cache_function_init(m,igrid,hash_key)
       use amr_commons, only: mesh_t
       use amr_parameters, only: ndim
       type(mesh_t)::m
       integer::igrid
       integer(kind=8),dimension(0:ndim)::hash_key
     end subroutine cache_function_init
     subroutine cache_function_unpack(m,igrid,msg_size,msg_array,hash_key)
       use amr_commons, only: mesh_t
       use amr_parameters, only: ndim
       type(mesh_t)::m
       integer::igrid
       integer::msg_size
       integer,dimension(1:msg_size),optional::msg_array
       integer(kind=8),dimension(0:ndim)::hash_key
     end subroutine cache_function_unpack

     subroutine cache_function_bound(r,g,m,igrid,igrid_ref,ibound)
       use amr_commons, only: run_t, global_t, mesh_t
       type(run_t)::r
       type(global_t)::g
       type(mesh_t)::m
       integer::igrid, igrid_ref
       integer::ibound
     end subroutine cache_function_bound

     subroutine cache_function_clump(c,local_peak_id,msg_size,msg_array)
       use clfind_commons, only: clump_t
       type(clump_t)::c
       integer::local_peak_id
       integer::msg_size
       integer,dimension(1:msg_size),optional::msg_array
     end subroutine cache_function_clump
     subroutine cache_function_init_clump(c,local_peak_id)
       use clfind_commons, only: clump_t
       type(clump_t)::c
       integer::local_peak_id
     end subroutine cache_function_init_clump
     subroutine cache_function_unpack_clump(c,local_peak_id,msg_size,msg_array)
       use clfind_commons, only: clump_t
       type(clump_t)::c
       integer::local_peak_id
       integer::msg_size
       integer,dimension(1:msg_size),optional::msg_array
     end subroutine cache_function_unpack_clump

  end interface

  type cache_f
     procedure(cache_function),pointer,nopass::proc
  end type cache_f
  type cache_init_f
     procedure(cache_function_init),pointer,nopass::proc
  end type cache_init_f
  type cache_unpack_f
     procedure(cache_function_unpack),pointer,nopass::proc
  end type cache_unpack_f
  type cache_bound_f
     procedure(cache_function_bound),pointer,nopass::proc
  end type cache_bound_f
  type cache_f_clump
     procedure(cache_function_clump),pointer,nopass::proc
  end type cache_f_clump
  type cache_init_f_clump
     procedure(cache_function_init_clump),pointer,nopass::proc
  end type cache_init_f_clump
  type cache_unpack_f_clump
     procedure(cache_function_unpack_clump),pointer,nopass::proc
  end type cache_unpack_f_clump
     
end module call_back
