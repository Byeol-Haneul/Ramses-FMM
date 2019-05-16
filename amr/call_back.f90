module call_back

  type in_make_initial_phi_t
    integer::ilevel, icount
  end type in_make_initial_phi_t

  type in_recurrence_t
    integer::ilevel
    real(kind=8)::cg
  end type in_recurrence_t

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
       integer,optional::input_size
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

  procedure(cache_function)::pack_fetch_flag,pack_fetch_refine,pack_fetch_split,pack_fetch_hydro,pack_fetch_cg,pack_fetch_phi
  procedure(cache_function)::pack_fetch_restrict_res,pack_fetch_scan,pack_fetch_mg,pack_fetch_interpol,pack_fetch_kick
  procedure(cache_function)::unpack_fetch_flag,unpack_fetch_refine,unpack_fetch_split,unpack_fetch_hydro,unpack_fetch_cg,unpack_fetch_phi
  procedure(cache_function)::unpack_fetch_restrict_res,unpack_fetch_scan,unpack_fetch_mg,unpack_fetch_interpol,unpack_fetch_kick
  procedure(cache_function)::init_flush_initflag,init_flush_derefine,init_flush_upload,init_flush_multipole,init_flush_rho
  procedure(cache_function)::init_flush_restrict_mask,init_flush_restrict_res,init_flush_godunov
  procedure(cache_function)::pack_flush_initflag,pack_flush_derefine,pack_flush_upload,pack_flush_multipole,pack_flush_rho
  procedure(cache_function)::pack_flush_build_mg,pack_flush_restrict_mask,pack_flush_restrict_res,pack_flush_godunov,pack_flush_refine
  procedure(cache_function)::pack_flush_loadbalance
  procedure(cache_function)::unpack_flush_initflag,unpack_flush_derefine,unpack_flush_upload,unpack_flush_multipole,unpack_flush_rho
  procedure(cache_function)::unpack_flush_build_mg,unpack_flush_restrict_mask,unpack_flush_restrict_res,unpack_flush_godunov,unpack_flush_refine
  procedure(cache_function)::unpack_flush_loadbalance

  type cache_f
     procedure(cache_function),pointer,nopass::proc
  end type cache_f
     
end module call_back
