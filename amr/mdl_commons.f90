module mdl_parameters

  ! Parameter for call-back functions indices
  enum,bind(C)
     enumerator::MDL_CLEAN_STOP
     enumerator::MDL_SET_ADD
     enumerator::MDL_BCAST_PARAMS
     enumerator::MDL_BCAST_GLOBAL
     enumerator::MDL_INIT_AMR
     enumerator::MDL_INIT_TIME
     enumerator::MDL_INIT_HYDRO
     enumerator::MDL_INIT_PART
     enumerator::MDL_INPUT_PART_GRAFIC
     enumerator::MDL_INPUT_PART_ASCII
     enumerator::MDL_INPUT_PART_RESTART
     enumerator::MDL_NPART_MAX
     enumerator::MDL_INIT_FLAG
     enumerator::MDL_USER_FLAG
     enumerator::MDL_ENSURE_REF_RULES
     enumerator::MDL_COLLECT_NOCT
     enumerator::MDL_NOCT_TOT
     enumerator::MDL_NOCT_MIN
     enumerator::MDL_NOCT_MAX
     enumerator::MDL_NOCT_USED_MAX
     enumerator::MDL_GATHER_NOCT_MAX
     enumerator::MDL_INIT_REFINE_BASEGRID
     enumerator::MDL_INIT_REFINE_RESTART
     enumerator::MDL_COLLECT_BOUND_KEY
     enumerator::MDL_BROADCAST_BOUND_KEY
     enumerator::MDL_LOAD_BALANCE
     enumerator::MDL_BALANCE_PART
     enumerator::MDL_REFINE_FINE
     enumerator::MDL_SMOOTH_FINE
     enumerator::MDL_INPUT_HYDRO_CONDINIT
     enumerator::MDL_INPUT_HYDRO_GRAFIC
     enumerator::MDL_UPLOAD_FINE
     enumerator::MDL_MULTIPOLE_LEAF_CELLS
     enumerator::MDL_MULTIPOLE_SPLIT_CELLS
     enumerator::MDL_RESET_RHO
     enumerator::MDL_CIC_MULTIPOLE
     enumerator::MDL_CIC_PART
     enumerator::MDL_SPLIT_PART
     enumerator::MDL_KICK_DRIFT_PART
     enumerator::MDL_MASS_MIN_PART
     enumerator::MDL_BROADCAST_MP_MIN
     enumerator::MDL_COLLECT_MULTIPOLE
     enumerator::MDL_BROADCAST_MULTIPOLE
     enumerator::MDL_OUTPUT_AMR
     enumerator::MDL_OUTPUT_HYDRO
     enumerator::MDL_OUTPUT_POISSON
     enumerator::MDL_OUTPUT_PART
     enumerator::MDL_SYNCHRO_HYDRO_FINE
     enumerator::MDL_SAVE_PHI_OLD
     enumerator::MDL_FORCE_ANALYTIC
     enumerator::MDL_GRADIENT_PHI
     enumerator::MDL_COMPUTE_EPOT
     enumerator::MDL_COMPUTE_RHOMAX
     enumerator::MDL_BROADCAST_AEXP
     enumerator::MDL_COURANT_FINE
     enumerator::MDL_GODUNOV_FINE
     enumerator::MDL_SET_UNEW
     enumerator::MDL_SET_UOLD
     enumerator::MDL_GRAVITY_HYDRO_FINE
     enumerator::MDL_COOLING_FINE
     enumerator::MDL_NEWDT_PART
     enumerator::MDL_BROADCAST_DT
     enumerator::MDL_MAKE_INITIAL_PHI
     enumerator::MDL_RECURRENCE_ON_P
     enumerator::MDL_RECURRENCE_X_AND_R
     enumerator::MDL_CMP_RESIDUAL_CG
     enumerator::MDL_CMP_AP_CG
     enumerator::MDL_CMP_RHS_NORM
     enumerator::MDL_CMP_R2_CG
     enumerator::MDL_CMP_PAP_CG
     enumerator::MDL_INIT_MG
     enumerator::MDL_BUILD_MG
     enumerator::MDL_CLEANUP_MG
     enumerator::MDL_MAKE_MASK
     enumerator::MDL_MAKE_BC_RHS
     enumerator::MDL_RESTRICT_MASK
     enumerator::MDL_CMP_RESIDUAL_MG
     enumerator::MDL_GAUSS_SEIDEL_MG
     enumerator::MDL_RESET_CORRECTION
     enumerator::MDL_RESTRICT_RESIDUAL
     enumerator::MDL_INTERPOLATE_AND_CORRECT
     enumerator::MDL_SET_SCAN_FLAG
     enumerator::MDL_CMP_RESIDUAL_NORM2
     enumerator::MDL_OUTPUT_FRAME
  end enum
  
  ! Maximum number of cpus
  integer,parameter::MDL_MAX_CPU=262144
  
end module mdl_parameters

module mdl_commons
  
  use mdl_parameters
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_FUNPTR

  type mdl_t
     
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
     
  end type mdl_t

  contains
    subroutine mdl_add_service(mdl,sid,p1,service,nin,nout)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_FUNPTR
      type(mdl_t)::mdl
      integer::sid
      type(*)::p1
      type(c_funptr), intent(in), value :: service
      integer :: nin, nout ! NOTE: nin is the size IN BYTES!!
      mdl%callback(sid) = service
      mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,nin/4) ! Divide by four to get the number of Integers

    end subroutine mdl_add_service


!     callback(MDL_BCAST_PARAMS)%proc => r_broadcast_params
!     mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,storage_size(pst%s%r)/32)

end module mdl_commons
