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

  procedure(cache_function)::pack_fetch_hydro,pack_fetch_phi
  procedure(cache_function)::pack_fetch_restrict_res,pack_fetch_scan,pack_fetch_mg
  procedure(cache_function)::unpack_fetch_hydro,unpack_fetch_phi
  procedure(cache_function)::unpack_fetch_restrict_res,unpack_fetch_scan,unpack_fetch_mg
  procedure(cache_function)::init_flush_restrict_mask,init_flush_restrict_res
  procedure(cache_function)::pack_flush_build_mg,pack_flush_restrict_mask,pack_flush_restrict_res
  procedure(cache_function)::unpack_flush_build_mg,unpack_flush_restrict_mask,unpack_flush_restrict_res

  type cache_f
     procedure(cache_function),pointer,nopass::proc
  end type cache_f
     
end module call_back
