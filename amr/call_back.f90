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
     subroutine cache_function(grid,msg_size,msg_array)
       use amr_commons, only: oct
       type(oct)::grid
       integer::msg_size
       integer,dimension(1:msg_size),optional::msg_array
     end subroutine cache_function
     subroutine cache_function_init(grid,hash_key)
       use amr_commons, only: oct
       use amr_parameters, only: ndim
       type(oct)::grid
       integer(kind=8),dimension(0:ndim)::hash_key
     end subroutine cache_function_init
     subroutine cache_function_unpack(grid,msg_size,msg_array,hash_key)
       use amr_commons, only: oct
       use amr_parameters, only: ndim
       type(oct)::grid
       integer::msg_size
       integer,dimension(1:msg_size),optional::msg_array
       integer(kind=8),dimension(0:ndim)::hash_key
     end subroutine cache_function_unpack

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
     
end module call_back
