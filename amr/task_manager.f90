!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine mdl_init
  use amr_parameters, only: flen
  use ramses_commons, only: pst_t, ramses_t
  use mdl_commons
  use call_back
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  
  type(pst_t)::pst
  type(ramses_t),target::s
  type(call_back_f),dimension(0:100)::callback
  integer,dimension(1)::ncpu
  
  pst%s => s
  
  associate(mdl=>pst%s%mdl)

  ! MPI initialization
#ifndef WITHOUTMPI
  call MPI_INIT(info)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mdl%myid,info)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mdl%ncpu,info)
  mdl%myid=mdl%myid+1 ! Careful with this...
  if(mdl%myid==1)then
     write(*,'(" Launching MPI with nproc = ",I4)')mdl%ncpu
  endif
#else
  write(*,'(" Serial execution (no MPI).")')
  mdl%ncpu=1
  mdl%myid=1
#endif

  ! Store cpu info as a global variable
  pst%s%g%myid=mdl%myid
  pst%s%g%ncpu=mdl%ncpu

  ! Register call-back functions

#ifndef WITHOUTMPI

  callback(MDL_CLEAN_STOP)%proc => r_clean_stop

  callback(MDL_SET_ADD)%proc => r_set_add

  callback(MDL_BCAST_PARAMS)%proc => r_broadcast_params
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,storage_size(pst%s%r)/32)

  callback(MDL_BCAST_GLOBAL)%proc => r_broadcast_global
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,storage_size(pst%s%g)/32)

  callback(MDL_INIT_AMR)%proc => r_init_amr

  callback(MDL_INIT_TIME)%proc => r_init_time

  callback(MDL_INIT_HYDRO)%proc => r_init_hydro

  callback(MDL_INIT_PART)%proc => r_init_part

  callback(MDL_INPUT_PART_GRAFIC)%proc => r_input_part_grafic
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,storage_size(pst%s%p%npart_tot)/32)

  callback(MDL_INPUT_PART_ASCII)%proc => r_input_part_ascii
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,storage_size(pst%s%p%npart_tot)/32)

  callback(MDL_INPUT_PART_RESTART)%proc => r_input_part_restart
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,MDL_MAX_CPU)

  callback(MDL_NPART_MAX)%proc => r_npart_max

  callback(MDL_INIT_FLAG)%proc => r_init_flag

  callback(MDL_USER_FLAG)%proc => r_user_flag

  callback(MDL_ENSURE_REF_RULES)%proc => r_ensure_ref_rules

  callback(MDL_COLLECT_NOCT)%proc => r_collect_noct

  callback(MDL_NOCT_TOT)%proc => r_noct_tot

  callback(MDL_NOCT_MIN)%proc => r_noct_min

  callback(MDL_NOCT_MAX)%proc => r_noct_max

  callback(MDL_NOCT_USED_MAX)%proc => r_noct_used_max

  callback(MDL_GATHER_NOCT_MAX)%proc => r_gather_noct_max

  callback(MDL_INIT_REFINE_BASEGRID)%proc => r_init_refine_basegrid

  callback(MDL_INIT_REFINE_RESTART)%proc => r_init_refine_restart

  callback(MDL_COLLECT_BOUND_KEY)%proc => r_collect_bound_key
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,MDL_MAX_CPU+1)

  callback(MDL_BROADCAST_BOUND_KEY)%proc => r_broadcast_bound_key
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,(storage_size(pst%s%m%domain)+32)/32)

  callback(MDL_LOAD_BALANCE)%proc => r_load_balance

  callback(MDL_BALANCE_PART)%proc => r_balance_part

  callback(MDL_REFINE_FINE)%proc => r_refine_fine

  callback(MDL_SMOOTH_FINE)%proc => r_smooth_fine

  callback(MDL_INPUT_HYDRO_CONDINIT)%proc => r_input_hydro_condinit

  callback(MDL_INPUT_HYDRO_GRAFIC)%proc => r_input_hydro_grafic

  callback(MDL_UPLOAD_FINE)%proc => r_upload_fine

  callback(MDL_MULTIPOLE_LEAF_CELLS)%proc => r_multipole_leaf_cells

  callback(MDL_MULTIPOLE_SPLIT_CELLS)%proc => r_multipole_split_cells

  callback(MDL_RESET_RHO)%proc => r_reset_rho

  callback(MDL_CIC_MULTIPOLE)%proc => r_cic_multipole

  callback(MDL_CIC_PART)%proc => r_cic_part

  callback(MDL_SPLIT_PART)%proc => r_split_part

  callback(MDL_KICK_DRIFT_PART)%proc => r_kick_drift_part

  callback(MDL_MASS_MIN_PART)%proc => r_mass_min_part

  callback(MDL_BROADCAST_MP_MIN)%proc => r_broadcast_mp_min

  callback(MDL_COLLECT_MULTIPOLE)%proc => r_collect_multipole

  callback(MDL_BROADCAST_MULTIPOLE)%proc => r_broadcast_multipole
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,storage_size(pst%s%g%multipole)/32)

  callback(MDL_OUTPUT_AMR)%proc => r_output_amr
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,flen/4)

  callback(MDL_OUTPUT_HYDRO)%proc => r_output_hydro
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,flen/4)

  callback(MDL_OUTPUT_POISSON)%proc => r_output_poisson
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,flen/4)

  callback(MDL_OUTPUT_PART)%proc => r_output_part
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,flen/4)

  callback(MDL_SYNCHRO_HYDRO_FINE)%proc => r_synchro_hydro_fine
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,3)

  callback(MDL_SAVE_PHI_OLD)%proc => r_save_phi_old
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,1)

  callback(MDL_FORCE_ANALYTIC)%proc => r_force_analytic
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,1)

  callback(MDL_GRADIENT_PHI)%proc => r_gradient_phi
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,2)

  callback(MDL_COMPUTE_EPOT)%proc => r_compute_epot
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,1)

  callback(MDL_COMPUTE_RHOMAX)%proc => r_compute_rhomax
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,1)

  callback(MDL_BROADCAST_AEXP)%proc => r_broadcast_aexp

  callback(MDL_COURANT_FINE)%proc => r_courant_fine

  callback(MDL_GODUNOV_FINE)%proc => r_godunov_fine
  
  callback(MDL_SET_UNEW)%proc => r_set_unew
  
  callback(MDL_SET_UOLD)%proc => r_set_uold
  
  callback(MDL_GRAVITY_HYDRO_FINE)%proc => r_gravity_hydro_fine
  
  callback(MDL_COOLING_FINE)%proc => r_cooling_fine
  
  callback(MDL_NEWDT_PART)%proc => r_newdt_part
  
  callback(MDL_BROADCAST_DT)%proc => r_broadcast_dt
  
  callback(MDL_MAKE_INITIAL_PHI)%proc => r_make_initial_phi
  
  callback(MDL_RECURRENCE_ON_P)%proc => r_recurrence_on_p
  
  callback(MDL_RECURRENCE_X_AND_R)%proc => r_recurrence_x_and_r
  
  callback(MDL_CMP_RESIDUAL_CG)%proc => r_cmp_residual_cg
  
  callback(MDL_CMP_R2_CG)%proc => r_cmp_r2_cg
  
  callback(MDL_CMP_PAP_CG)%proc => r_cmp_pAp_cg
  
  callback(MDL_CMP_RHS_NORM)%proc => r_cmp_rhs_norm
  
  callback(MDL_INIT_MG)%proc => r_init_mg

  callback(MDL_BUILD_MG)%proc => r_build_mg
  
  callback(MDL_CLEANUP_MG)%proc => r_cleanup_mg
  
  callback(MDL_MAKE_MASK)%proc => r_make_mask
  
  callback(MDL_MAKE_BC_RHS)%proc => r_make_bc_rhs
  
  callback(MDL_RESTRICT_MASK)%proc => r_restrict_mask

  callback(MDL_CMP_RESIDUAL_MG)%proc => r_cmp_residual_mg

  callback(MDL_GAUSS_SEIDEL_MG)%proc => r_gauss_seidel_mg

  callback(MDL_RESET_CORRECTION)%proc => r_reset_correction

  callback(MDL_RESTRICT_RESIDUAL)%proc => r_restrict_residual

  callback(MDL_INTERPOLATE_AND_CORRECT)%proc => r_interpolate_and_correct

  callback(MDL_SET_SCAN_FLAG)%proc => r_set_scan_flag

  callback(MDL_CMP_RESIDUAL_NORM2)%proc => r_cmp_residual_norm2

  callback(MDL_OUTPUT_FRAME)%proc => r_output_frame

  ! Allocate input and output buffer sizes
  allocate(mdl%mpi_input_buffer(1:32+mdl%MDL_INPUT_MAXSIZE))

  ! Initialize software cache
  call init_cache(mdl)

#endif
  
  ! For slave workers, go into waiting loop
  if(mdl%myid>1)then
     call mdl_wait(pst,callback)
  else
     ncpu(1)=mdl%ncpu
     call r_set_add(pst,1,0,ncpu)
     call adaptive_loop(pst)
  endif

#ifndef WITHOUTMPI
  write(*,*)"MYID ",mdl%myid," TERMINATING AND EXITING"
#else
  write(*,*)"TERMINATING AND EXITING"
#endif  
  
#ifndef WITHOUTMPI
  call MPI_FINALIZE(info)
#endif

  end associate
  
end subroutine mdl_init
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine mdl_wait(pst,callback)
  use ramses_commons, only: pst_t
  use call_back, only: call_back_f
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  type(pst_t)::pst
  type(call_back_f),dimension(0:100)::callback

#ifndef WITHOUTMPI

  logical::stop_order_received,order_received
  integer::order_id,order_tag=101,output_tag=203,output_id
  integer,dimension(MPI_STATUS_SIZE)::order_status,output_status
  integer::input_size,output_size,source_cpu,function_id
  integer,dimension(:),allocatable::input_array,output_array
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
        
        if(input_size>0)then
           allocate(input_array(1:input_size))
           input_array(1:input_size)=mdl%mpi_input_buffer(33:32+input_size)
        endif
        
        if(output_size>0)then
           allocate(output_array(1:output_size))
           output_array=0
        endif
        
        ! Launch the corresponding call-back function
        call callback(function_id)%proc(pst,input_size,output_size,input_array,output_array)
        
        ! Deallocate input array
        if(input_size>0)then
           deallocate(input_array)
        endif
        
        ! Send the output back to the source cpu
        if(output_size>0)then
           call MPI_ISEND(output_array,output_size,MPI_INTEGER,source_cpu,output_tag,MPI_COMM_WORLD,output_id,info)
           call MPI_WAIT(output_id,output_status,info)
           deallocate(output_array)
        endif
        
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
subroutine mdl_send_request(mdl,mdl_function_id,target_cpu,input_size,output_size,input_array)
  use mdl_commons, only: mdl_t
  implicit none
  type(mdl_t)::mdl
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  integer::mdl_function_id
  integer::target_cpu,input_size,output_size
  integer,dimension(1:input_size),optional::input_array

#ifndef WITHOUTMPI
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

end subroutine mdl_send_request
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine mdl_get_reply(mdl,target_cpu,output_size,output_array)
  use mdl_commons, only: mdl_t
  implicit none
  type(mdl_t)::mdl
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  integer::target_cpu,output_size
  integer,dimension(1:output_size)::output_array

#ifndef WITHOUTMPI  
  integer::output_tag=203,output_id
  integer,dimension(MPI_STATUS_SIZE)::output_status  
  
  ! Post a RECV for the output back from target_cpu
  call MPI_IRECV(output_array,output_size,MPI_INTEGER,target_cpu-1,output_tag,MPI_COMM_WORLD,output_id,info)

  ! Wait for ISEND completion to free memory in corresponding MPI buffer
  call MPI_WAIT(output_id,output_status,info)

#endif  

end subroutine mdl_get_reply
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine mdl_abort
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
  call MPI_ABORT(MPI_COMM_WORLD,info)
#else
  stop
#endif  
end subroutine mdl_abort
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine init_cache(mdl)
  use mdl_commons, only: mdl_t
  use cache_commons
  implicit none
  type(mdl_t)::mdl
  
  integer::ncpu
  type(msg_large_realdp)::dummy_large_realdp

  ncpu=mdl%ncpu
  
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
  allocate(mdl%send_fetch_array(1:ncpu*mdl%size_fetch_array))
  allocate(mdl%send_flush_array(1:ncpu*mdl%size_flush_array))

#endif

end subroutine init_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################












