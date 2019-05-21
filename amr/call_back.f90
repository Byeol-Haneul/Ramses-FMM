module call_back

  interface
!     recursive subroutine ramses_function(pst,input,input_size,output,output_size)
!       use ramses_commons, only: pst_t
!       type(pst_t)::pst
!       integer::input_size,output_size
!       integer,dimension(1:input_size)::input
!       integer,dimension(1:output_size)::output
!     end subroutine ramses_function
     recursive subroutine ramses_function(pst,input,input_size,output,output_size)
       use ramses_commons, only: pst_t
       type(pst_t)::pst
       integer::input_size
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
  end interface

  type cache_f
     procedure(cache_function),pointer,nopass::proc
  end type cache_f
     
end module call_back
