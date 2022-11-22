!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine mdl_init
  use call_back
  use amr_parameters, only: flen
  use mdl_module
#ifdef MDL2
  use ramses_commons, only: pst_t
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_NULL_PTR, C_FUNLOC
  call mdl_launch( command_argument_count(), C_NULL_PTR, C_FUNLOC(master), C_FUNLOC(worker_init), C_FUNLOC(worker_done) )
#else
  use ramses_commons, only: pst_t, ramses_t
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_FUNLOC, C_SIZEOF
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  type(mdl_t),pointer::mdl
  type(pst_t),allocatable::pst
  integer,dimension(1)::ncpu

  allocate(mdl)
  call mdl_initialize(mdl)
!  associate(mdl=>pst%s%mdl)

  pst = worker_init(mdl)

#ifndef WITHOUTMPI
  ! Allocate input and output buffer sizes
  allocate(mdl%mpi_input_buffer(1:32+mdl%MDL_INPUT_MAXSIZE))

  ! Initialize software cache
  call init_cache(mdl)
#endif
  
  ! For slave workers, go into waiting loop
  if(mdl_self(mdl)>1)then
     call mdl_wait(pst)
  else
     call master(mdl,pst)
  endif

! #ifndef WITHOUTMPI
!   write(*,*)"MYID ",mdl_self(mdl)," TERMINATING AND EXITING"
! #else
!   write(*,*)"TERMINATING AND EXITING"
! #endif  
  
#ifndef WITHOUTMPI
  call MPI_FINALIZE(info)
#endif

!  end associate
#endif

contains
!##############################################################
!##############################################################
!##############################################################
!##############################################################
#ifndef MDL2
recursive subroutine r_clean_stop(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CLEAN_STOP,pst%iUpper+1)
     call r_clean_stop(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  endif
  
end subroutine r_clean_stop
#endif
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine master(mdl,pst)
  use init_amr_module, only: r_set_add
  implicit none
  type(*)::mdl
  type(pst_t)::pst
  call r_set_add(pst,mdl_threads(pst%s%mdl),1)
  call adaptive_loop(pst)
#ifndef MDL2
  call r_clean_stop(pst)
#endif
end subroutine master
!##############################################################
!##############################################################
!##############################################################
!##############################################################
function worker_init(mdl) result(pst)
  use ramses_commons, only: pst_t, ramses_t
  use call_back, only: ramses_function
  use mdl_parameters
  use flag_utils, only: r_init_flag, r_ensure_ref_rules, r_user_flag
  use init_amr_module, only: r_init_amr, r_set_add
  use params_module, only: r_broadcast_params,r_broadcast_global
  use init_time_module, only: r_init_time
  use init_hydro_module, only: r_init_hydro
  use init_part_module, only: r_init_part
  use input_part_grafic_module, only: r_input_part_grafic
  use input_part_ascii_module, only: r_input_part_ascii
  use input_part_restart_module, only: r_input_part_restart
  use input_part_module, only: r_npart_max, r_mass_min_part, r_broadcast_mp_min
  use update_time_module, only: r_broadcast_aexp
  use init_refine_basegrid_module, only:r_init_refine_basegrid,r_collect_noct,r_noct_tot,r_noct_min,r_noct_max,&
                                        r_noct_used_max,r_gather_noct_max
  use init_refine_restart_module, only: r_init_refine_restart
  use load_balance_module, only: r_load_balance,r_balance_part,r_broadcast_bound_key,r_collect_bound_key
  use refine_utils, only: r_refine_fine
  use smooth_module, only: r_smooth_fine
  use input_hydro_condinit_module, only: r_input_hydro_condinit
  use input_hydro_grafic_module, only: r_input_hydro_grafic
  use upload_module, only: r_upload_fine
  use move_fine_module, only: r_kick_drift_part
  use output_amr_module, only: r_output_amr
  use output_hydro_module, only: r_output_hydro
  use output_poisson_module, only: r_output_poisson
  use output_part_module, only: r_output_part
  use synchro_hydro_fine_module, only: r_synchro_hydro_fine, r_gravity_hydro_fine
  use nbors_utils_p, only: r_save_phi_old
  use courant_fine_module, only: r_courant_fine
  use godunov_fine_module, only: r_godunov_fine,r_set_unew,r_set_uold
  use cooling_fine_module, only: r_cooling_fine
  use newdt_fine_module, only: r_newdt_part,r_broadcast_dt
#ifdef GRAV
  use force_fine_module, only: r_force_analytic,r_compute_epot,r_compute_rhomax,r_gradient_phi
  use rho_fine_module, only: r_multipole_leaf_cells,r_multipole_split_cells,r_broadcast_multipole,r_collect_multipole,&
                            r_cic_multipole,r_cic_part,r_split_part,r_reset_rho
  use phi_fine_cg_module, only: r_cmp_pAp_cg,r_cmp_r2_cg,r_cmp_residual_cg,r_cmp_rhs_norm,&
                                r_make_initial_phi,r_recurrence_on_p,r_recurrence_x_and_r
  use multigrid_fine_commons, only: r_init_mg,r_build_mg,r_cleanup_mg,r_make_mask,r_make_bc_rhs
  use multigrid_fine_coarse, only: r_restrict_mask,r_cmp_residual_mg,r_cmp_residual_norm2,r_restrict_residual,&
                                r_reset_correction,r_set_scan_flag,r_gauss_seidel_mg,r_interpolate_and_correct
#endif
  use movie_module, only: r_output_frame
  use amr_parameters, only: nhilbert

  implicit none

  integer::dummy
  integer(kind=8)::dummy8
  integer::ncpu
#ifdef MDL2
  type(c_ptr),value::mdl
  type(mdl_t),pointer::mdl_wrapper
#else
  type(mdl_t),target::mdl
#endif
  type(pst_t),allocatable::pst

  allocate(pst)
  allocate(pst%s)
#ifdef MDL2
  allocate(mdl_wrapper)
  mdl_wrapper%mdl2 = mdl
  pst%s%mdl => mdl_wrapper
#else
  pst%s%mdl => mdl
#endif
  ncpu = mdl_threads(pst%s%mdl)

  ! Store cpu info as a global variable
  pst%s%g%myid=mdl_self(pst%s%mdl)
  pst%s%g%ncpu=mdl_threads(pst%s%mdl)

  ! Input and Output sizes here are in BYTES!
#ifndef MDL2
  call mdl_add_service(pst%s%mdl,MDL_CLEAN_STOP,             pst,C_FUNLOC(r_clean_stop),0,0,"clean_stop")
#endif
  call mdl_add_service(pst%s%mdl,MDL_SET_ADD,                pst,C_FUNLOC(r_set_add),storage_size(dummy)/8,0, "set_add")
  call mdl_add_service(pst%s%mdl,MDL_BCAST_PARAMS,           pst,C_FUNLOC(r_broadcast_params),storage_size(pst%s%r)/8,0,"broadcast_params")
  call mdl_add_service(pst%s%mdl,MDL_BCAST_GLOBAL,           pst,C_FUNLOC(r_broadcast_global),storage_size(pst%s%g)/8,0,"broadcast_global")
  call mdl_add_service(pst%s%mdl,MDL_INIT_AMR,               pst,C_FUNLOC(r_init_amr),0,0,"init_amr")
  call mdl_add_service(pst%s%mdl,MDL_INIT_TIME,              pst,C_FUNLOC(r_init_time),0,0,"init_time")
  call mdl_add_service(pst%s%mdl,MDL_INIT_HYDRO,             pst,C_FUNLOC(r_init_hydro),0,0,"init_hydro")
  call mdl_add_service(pst%s%mdl,MDL_INIT_PART,              pst,C_FUNLOC(r_init_part),0,0,"init_part")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_PART_GRAFIC,      pst,C_FUNLOC(r_input_part_grafic),storage_size(pst%s%p%npart_tot)/8,0,"input_part_grafic")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_PART_ASCII,       pst,C_FUNLOC(r_input_part_ascii),storage_size(pst%s%p%npart_tot)/8,0,"input_part_ascii")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_PART_RESTART,     pst,C_FUNLOC(r_input_part_restart),MDL_MAX_CPU*4,0,"input_part_restart") ! number of integers
  call mdl_add_service(pst%s%mdl,MDL_NPART_MAX,              pst,C_FUNLOC(r_npart_max),0,4,"npart_max")
  call mdl_add_service(pst%s%mdl,MDL_INIT_FLAG,              pst,C_FUNLOC(r_init_flag),4,4,"init_flag")
  call mdl_add_service(pst%s%mdl,MDL_USER_FLAG,              pst,C_FUNLOC(r_user_flag),4,4,"user_flag")
  call mdl_add_service(pst%s%mdl,MDL_ENSURE_REF_RULES,       pst,C_FUNLOC(r_ensure_ref_rules),4,0,"ensure_ref_rules")
  call mdl_add_service(pst%s%mdl,MDL_COLLECT_NOCT,           pst,C_FUNLOC(r_collect_noct),4,ncpu*storage_size(dummy)/8,"collect_noct")
  call mdl_add_service(pst%s%mdl,MDL_NOCT_TOT,               pst,C_FUNLOC(r_noct_tot),4,4,"noct_tot")
  call mdl_add_service(pst%s%mdl,MDL_NOCT_MIN,               pst,C_FUNLOC(r_noct_min),4,4,"noct_min")
  call mdl_add_service(pst%s%mdl,MDL_NOCT_MAX,               pst,C_FUNLOC(r_noct_max),4,4,"noct_max")
  call mdl_add_service(pst%s%mdl,MDL_NOCT_USED_MAX,          pst,C_FUNLOC(r_noct_used_max),4,4,"noct_used_max")
  call mdl_add_service(pst%s%mdl,MDL_GATHER_NOCT_MAX,        pst,C_FUNLOC(r_gather_noct_max),0,0,"gather_noct_max")
  call mdl_add_service(pst%s%mdl,MDL_INIT_REFINE_BASEGRID,   pst,C_FUNLOC(r_init_refine_basegrid),4,0,"init_refine_basegrid")
  call mdl_add_service(pst%s%mdl,MDL_INIT_REFINE_RESTART,    pst,C_FUNLOC(r_init_refine_restart),0,2*nhilbert*(pst%s%g%ncpu+1)*4,"init_refine_restart")
  call mdl_add_service(pst%s%mdl,MDL_COLLECT_BOUND_KEY,      pst,C_FUNLOC(r_collect_bound_key),(MDL_MAX_CPU+1)*4,nhilbert*(ncpu+1)*storage_size(dummy8)/8,"collect_bound_key")
  call mdl_add_service(pst%s%mdl,MDL_BROADCAST_BOUND_KEY,    pst,C_FUNLOC(r_broadcast_bound_key),nhilbert*(ncpu+1)*storage_size(dummy8)/8 + 4,0,"broadcast_bound_key")
  call mdl_add_service(pst%s%mdl,MDL_LOAD_BALANCE,           pst,C_FUNLOC(r_load_balance),4,0,"load_balance")
  call mdl_add_service(pst%s%mdl,MDL_BALANCE_PART,           pst,C_FUNLOC(r_balance_part),0,0,"balance_part")
  call mdl_add_service(pst%s%mdl,MDL_REFINE_FINE,            pst,C_FUNLOC(r_refine_fine),4,8,"refine_fine")
  call mdl_add_service(pst%s%mdl,MDL_SMOOTH_FINE,            pst,C_FUNLOC(r_smooth_fine),4,4,"smooth_fine")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_HYDRO_CONDINIT,   pst,C_FUNLOC(r_input_hydro_condinit),4,0,"input_hydro_condinit")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_HYDRO_GRAFIC,     pst,C_FUNLOC(r_input_hydro_grafic),0,0,"input_hydro_grafic")
  call mdl_add_service(pst%s%mdl,MDL_UPLOAD_FINE,            pst,C_FUNLOC(r_upload_fine),4,0,"upload_fine")
  call mdl_add_service(pst%s%mdl,MDL_KICK_DRIFT_PART,        pst,C_FUNLOC(r_kick_drift_part),0,0,"kick_drift_part")
  call mdl_add_service(pst%s%mdl,MDL_MASS_MIN_PART,          pst,C_FUNLOC(r_mass_min_part),0,8,"mass_min_part")
  call mdl_add_service(pst%s%mdl,MDL_BROADCAST_MP_MIN,       pst,C_FUNLOC(r_broadcast_mp_min),0,0,"broadcast_mp_min")
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_AMR,             pst,C_FUNLOC(r_output_amr),flen,0,"output_amr")
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_HYDRO,           pst,C_FUNLOC(r_output_hydro),flen,0,"output_hydro")
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_POISSON,         pst,C_FUNLOC(r_output_poisson),flen,0,"output_poisson")
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_PART,            pst,C_FUNLOC(r_output_part),flen,0,"output_part")
  call mdl_add_service(pst%s%mdl,MDL_BROADCAST_AEXP,         pst,C_FUNLOC(r_broadcast_aexp),16,0,"broadcast_aexp")
  call mdl_add_service(pst%s%mdl,MDL_COURANT_FINE,           pst,C_FUNLOC(r_courant_fine),4,8*4,"courant_fine")
  call mdl_add_service(pst%s%mdl,MDL_GODUNOV_FINE,           pst,C_FUNLOC(r_godunov_fine),4,0,"godunov_fine")
  call mdl_add_service(pst%s%mdl,MDL_SET_UNEW,               pst,C_FUNLOC(r_set_unew),4,0,"set_unew")
  call mdl_add_service(pst%s%mdl,MDL_SET_UOLD,               pst,C_FUNLOC(r_set_uold),4,0,"set_uold")
  call mdl_add_service(pst%s%mdl,MDL_COOLING_FINE,           pst,C_FUNLOC(r_cooling_fine),4,0,"cooling_fine")
  call mdl_add_service(pst%s%mdl,MDL_NEWDT_PART,             pst,C_FUNLOC(r_newdt_part),0,0,"newdt_part")
  call mdl_add_service(pst%s%mdl,MDL_BROADCAST_DT,           pst,C_FUNLOC(r_broadcast_dt),24*4,0,"broadcast_dt")
#ifdef GRAV
  call mdl_add_service(pst%s%mdl,MDL_MULTIPOLE_LEAF_CELLS,   pst,C_FUNLOC(r_multipole_leaf_cells),4,0,"multipole_leaf_cells")
  call mdl_add_service(pst%s%mdl,MDL_MULTIPOLE_SPLIT_CELLS,  pst,C_FUNLOC(r_multipole_split_cells),4,0,"multipole_split_cells")
  call mdl_add_service(pst%s%mdl,MDL_RESET_RHO,              pst,C_FUNLOC(r_reset_rho),4,0,"reset_rho")
  call mdl_add_service(pst%s%mdl,MDL_CIC_MULTIPOLE,          pst,C_FUNLOC(r_cic_multipole),4,0,"cic_multipole")
  call mdl_add_service(pst%s%mdl,MDL_CIC_PART,               pst,C_FUNLOC(r_cic_part),0,0,"cic_part")
  call mdl_add_service(pst%s%mdl,MDL_COLLECT_MULTIPOLE,      pst,C_FUNLOC(r_collect_multipole),4,storage_size(pst%s%g%multipole)/8,"collect_multipole")
  call mdl_add_service(pst%s%mdl,MDL_COMPUTE_RHOMAX,         pst,C_FUNLOC(r_compute_rhomax),1*4,8,"compute_rhomax")
  call mdl_add_service(pst%s%mdl,MDL_SPLIT_PART,             pst,C_FUNLOC(r_split_part),0,0,"split_part")
  call mdl_add_service(pst%s%mdl,MDL_BROADCAST_MULTIPOLE,    pst,C_FUNLOC(r_broadcast_multipole),storage_size(pst%s%g%multipole)/8,0,"broadcast_multipole")
  call mdl_add_service(pst%s%mdl,MDL_SYNCHRO_HYDRO_FINE,     pst,C_FUNLOC(r_synchro_hydro_fine),3*4,0,"synchro_hydro_fine")
  call mdl_add_service(pst%s%mdl,MDL_SAVE_PHI_OLD,           pst,C_FUNLOC(r_save_phi_old),1*4,0,"save_phi_old")
  call mdl_add_service(pst%s%mdl,MDL_FORCE_ANALYTIC,         pst,C_FUNLOC(r_force_analytic),1*4,0,"force_analytic")
  call mdl_add_service(pst%s%mdl,MDL_GRADIENT_PHI,           pst,C_FUNLOC(r_gradient_phi),2*4,0,"gradient_phi")
  call mdl_add_service(pst%s%mdl,MDL_COMPUTE_EPOT,           pst,C_FUNLOC(r_compute_epot),1*4,8,"compute_epot")
  call mdl_add_service(pst%s%mdl,MDL_GRAVITY_HYDRO_FINE,     pst,C_FUNLOC(r_gravity_hydro_fine),4,0,"gravity_hydro_fine")
  call mdl_add_service(pst%s%mdl,MDL_MAKE_INITIAL_PHI,       pst,C_FUNLOC(r_make_initial_phi),8,0,"make_initial_phi")
  call mdl_add_service(pst%s%mdl,MDL_RECURRENCE_ON_P,        pst,C_FUNLOC(r_recurrence_on_p),0,0,"recurrence_on_p")
  call mdl_add_service(pst%s%mdl,MDL_RECURRENCE_X_AND_R,     pst,C_FUNLOC(r_recurrence_x_and_r),0,0,"recurrence_x_and_r")
  call mdl_add_service(pst%s%mdl,MDL_CMP_RESIDUAL_CG,        pst,C_FUNLOC(r_cmp_residual_cg),0,0,"cmp_residual_cg")
  call mdl_add_service(pst%s%mdl,MDL_CMP_R2_CG,              pst,C_FUNLOC(r_cmp_r2_cg),0,8,"cmp_r2_cg")
  call mdl_add_service(pst%s%mdl,MDL_CMP_PAP_CG,             pst,C_FUNLOC(r_cmp_pAp_cg),0,8,"cmp_pAp_cg")
  call mdl_add_service(pst%s%mdl,MDL_CMP_RHS_NORM,           pst,C_FUNLOC(r_cmp_rhs_norm),0,8,"cmp_rhs_norm")
  call mdl_add_service(pst%s%mdl,MDL_INIT_MG,                pst,C_FUNLOC(r_init_mg),4,0,"init_mg")
  call mdl_add_service(pst%s%mdl,MDL_BUILD_MG,               pst,C_FUNLOC(r_build_mg),4,0,"build_mg")
  call mdl_add_service(pst%s%mdl,MDL_CLEANUP_MG,             pst,C_FUNLOC(r_cleanup_mg),4,0,"cleanup_mg")
  call mdl_add_service(pst%s%mdl,MDL_MAKE_MASK,              pst,C_FUNLOC(r_make_mask),4,0,"make_mask")
  call mdl_add_service(pst%s%mdl,MDL_MAKE_BC_RHS,            pst,C_FUNLOC(r_make_bc_rhs),8,0,"make_bc_rhs")
  call mdl_add_service(pst%s%mdl,MDL_RESTRICT_MASK,          pst,C_FUNLOC(r_restrict_mask),4,4,"restrict_mask")
  call mdl_add_service(pst%s%mdl,MDL_CMP_RESIDUAL_MG,        pst,C_FUNLOC(r_cmp_residual_mg),8,0,"cmp_residual_mg")
  call mdl_add_service(pst%s%mdl,MDL_GAUSS_SEIDEL_MG,        pst,C_FUNLOC(r_gauss_seidel_mg),16,0,"gauss_seidel_mg")
  call mdl_add_service(pst%s%mdl,MDL_RESET_CORRECTION,       pst,C_FUNLOC(r_reset_correction),4,0,"reset_correction")
  call mdl_add_service(pst%s%mdl,MDL_RESTRICT_RESIDUAL,      pst,C_FUNLOC(r_restrict_residual),4,0,"restrict_residual")
  call mdl_add_service(pst%s%mdl,MDL_INTERPOLATE_AND_CORRECT,pst,C_FUNLOC(r_interpolate_and_correct),4,0,"interpolate_and_correct")
  call mdl_add_service(pst%s%mdl,MDL_SET_SCAN_FLAG,          pst,C_FUNLOC(r_set_scan_flag),8,0,"set_scan_flag")
  call mdl_add_service(pst%s%mdl,MDL_CMP_RESIDUAL_NORM2,     pst,C_FUNLOC(r_cmp_residual_norm2),4,8,"cmp_residual_norm2")
#endif
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_FRAME,           pst,C_FUNLOC(r_output_frame),8,&
                      2*pst%s%r%nw_frame*pst%s%r%nh_frame*4,"output_frame")
end function worker_init
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine worker_done(mdl,pst)
  type(mdl_t)::mdl
  type(pst_t),allocatable::pst
  ! FIXME: These should be deallocatable, but they aren't ???
  !deallocate(pst%s)
  !deallocate(pst)
end subroutine worker_done

end subroutine mdl_init
!##############################################################
!##############################################################
!##############################################################
!##############################################################
#ifndef MDL2
subroutine mdl_wait(pst)
  use ramses_commons, only: pst_t
  use call_back, only: call_back_f, ramses_function
  use mdl_module
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_F_PROCPOINTER, C_PTR, C_LOC, C_F_POINTER
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  type(pst_t)::pst

  type generic
  end type generic

  procedure(ramses_function),pointer::mdl_function
  type(generic),pointer::ipar,opar

#ifndef WITHOUTMPI

  logical::stop_order_received,order_received
  integer::order_id,order_tag=101,output_tag=203,output_id
  integer,dimension(MPI_STATUS_SIZE)::order_status,output_status
  integer::input_size
  integer::output_size,source_cpu,function_id
  integer,dimension(:),allocatable,target::input_array,output_array
  integer,dimension(1)::dummy
  integer,dimension(1:32)::header

  associate(s=>pst%s,mdl=>pst%s%mdl)
  
  ! Post the first RECV for a launch order
  call MPI_IRECV(mdl%mpi_input_buffer,mdl%MDL_INPUT_MAXSIZE+32,MPI_INTEGER,MPI_ANY_SOURCE,order_tag,MPI_COMM_WORLD,order_id,info)

  stop_order_received=.false.
  do while(.NOT. stop_order_received)

     call MPI_Test(order_id,order_received,order_status,info)
     
     if(order_received)then

        ! Execute call-back functions
        header=mdl%mpi_input_buffer(1:32)
        function_id=header(1)
        if(function_id==0)stop_order_received=.true.
        ! Get source cpu
        source_cpu=order_status(MPI_SOURCE)
        
        ! Allocate input and output arrays
        input_size=header(2)
        output_size=header(3)
        
!        if(input_size>0)then
           allocate(input_array(1:input_size))
           input_array(1:input_size)=mdl%mpi_input_buffer(33:32+input_size)
!        endif
        
!        if(output_size>0)then
           allocate(output_array(1:output_size))
           output_array=0
!        endif
        
        ! Launch the corresponding call-back function
        call c_f_pointer(c_loc(input_array),ipar)
        call c_f_pointer(c_loc(output_array),opar)
        CALL C_F_PROCPOINTER (mdl%callback(function_id), mdl_function)
        call mdl_function(pst,ipar,input_size,opar,output_size)
        
        ! Deallocate input array
!        if(input_size>0)then
        deallocate(input_array)
!        endif
        
        ! Always send the output back to the source cpu (even if not allocated = handshake)
        if(output_size>0)then
           call MPI_ISEND(output_array,output_size,MPI_INTEGER,source_cpu,output_tag,MPI_COMM_WORLD,output_id,info)
        else
           dummy(1)=1
           call MPI_ISEND(dummy,1,MPI_INTEGER,source_cpu,output_tag,MPI_COMM_WORLD,output_id,info)
        endif
        call MPI_WAIT(output_id,output_status,info)

        ! Deallocate output array
!        if(output_size>0)then
        deallocate(output_array)
!        endif
        
        ! Post a new RECV for the next launch order
        if(.NOT. stop_order_received)then
           call MPI_IRECV(mdl%mpi_input_buffer,mdl%MDL_INPUT_MAXSIZE+32,MPI_INTEGER,MPI_ANY_SOURCE,order_tag,MPI_COMM_WORLD,order_id,info)
        endif

     endif
  end do

  end associate

#endif

end subroutine mdl_wait
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine init_cache(mdl)
  use mdl_module
  use cache_commons
  implicit none
  type(mdl_t)::mdl
  
  integer::ncpu,ibuf
  type(msg_large_realdp)::dummy_large_realdp

  ncpu=mdl_threads(mdl)
  
#ifndef WITHOUTMPI  

  ! Allocate all communication and cache-related variables
  allocate(mdl%reply_id(1:ncpu))

  ! Compute largest possible message size
  mdl%size_request_array=1+ndim
  mdl%size_msg_array=storage_size(dummy_large_realdp)/32
  mdl%size_flush_array=1+(1+ndim+mdl%size_msg_array)*nflushmax
  mdl%size_fetch_array=2+(1+ndim+mdl%size_msg_array)*ntilemax

  ! Allocate large enough message communication buffers
  allocate(mdl%recv_request_array(1:mdl%size_request_array))
  allocate(mdl%send_request_array(1:mdl%size_request_array))
  allocate(mdl%recv_fetch_array(1:mdl%size_fetch_array))
  allocate(mdl%recv_flush_array(1:mdl%size_flush_array))

  allocate(mdl%cpu2buf_fetch(1:ncpu))
  allocate(mdl%cpu2buf_flush(1:ncpu))
  allocate(mdl%send_fetch(1:ncpu))
  allocate(mdl%send_flush(1:ncpu))

  mdl%nbuffer_fetch=MIN(ncpu,nbuffermax)
  mdl%nbuffer_flush=MIN(ncpu,nbuffermax)
  do ibuf=1,mdl%nbuffer_fetch
     allocate(mdl%send_fetch(ibuf)%array(mdl%size_fetch_array))
  end do
  do ibuf=1,mdl%nbuffer_flush
     allocate(mdl%send_flush(ibuf)%array(mdl%size_flush_array))
  end do

  mdl%ibuffer_fetch=0
  mdl%ibuffer_flush=0
  mdl%cpu2buf_fetch=0
  mdl%cpu2buf_flush=0

#endif

end subroutine init_cache
#endif
!##############################################################
!##############################################################
!##############################################################
!##############################################################
