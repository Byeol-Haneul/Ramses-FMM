!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine mdl_init
  use amr_parameters, only: flen
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons
  use call_back
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  type(call_back_f),dimension(0:100)::callback
  
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
  g%myid=mdl%myid
  g%ncpu=mdl%ncpu

  ! Register call-back functions

#ifndef WITHOUTMPI

  callback(MDL_CLEAN_STOP)%proc => r_clean_stop

  callback(MDL_BCAST_PARAMS)%proc => r_broadcast_params
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,storage_size(r)/32)

  callback(MDL_INIT_AMR)%proc => r_init_amr

  callback(MDL_INIT_TIME)%proc => r_init_time

  callback(MDL_INIT_HYDRO)%proc => r_init_hydro

  callback(MDL_INIT_PART)%proc => r_init_part

  callback(MDL_INPUT_PART_GRAFIC)%proc => r_input_part_grafic
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,storage_size(p%npart_tot)/32)

  callback(MDL_INPUT_PART_ASCII)%proc => r_input_part_ascii
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,storage_size(p%npart_tot)/32)

  callback(MDL_INPUT_PART_RESTART)%proc => r_input_part_restart
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,MDL_MAX_CPU)

  callback(MDL_INIT_FLAG)%proc => r_init_flag

  callback(MDL_USER_FLAG)%proc => r_user_flag

  callback(MDL_ENSURE_REF_RULES)%proc => r_ensure_ref_rules

  callback(MDL_COLLECT_NOCT)%proc => r_collect_noct

  callback(MDL_NOCT_TOT)%proc => r_noct_tot

  callback(MDL_NOCT_MIN)%proc => r_noct_min

  callback(MDL_NOCT_MAX)%proc => r_noct_max

  callback(MDL_GATHER_NOCT_MAX)%proc => r_gather_noct_max

  callback(MDL_INIT_REFINE_BASEGRID)%proc => r_init_refine_basegrid

  callback(MDL_INIT_REFINE_RESTART)%proc => r_init_refine_restart

  callback(MDL_COLLECT_BOUND_KEY)%proc => r_collect_bound_key
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,MDL_MAX_CPU+1)

  callback(MDL_BROADCAST_BOUND_KEY)%proc => r_broadcast_bound_key
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,(storage_size(m%domain)+32)/32)

  callback(MDL_LOAD_BALANCE)%proc => r_load_balance

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

  callback(MDL_COLLECT_MULTIPOLE)%proc => r_collect_multipole

  callback(MDL_BROADCAST_MULTIPOLE)%proc => r_broadcast_multipole
  mdl%MDL_INPUT_MAXSIZE=MAX(mdl%MDL_INPUT_MAXSIZE,storage_size(g%multipole)/32)

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
     call mdl_wait(r,g,m,p,mdl,callback)
  else
     call adaptive_loop(r,g,m,p,mdl)
  endif

#ifndef WITHOUTMPI
  write(*,*)"MYID ",mdl%myid," TERMINATING AND EXITING"
#else
  write(*,*)"TERMINATING AND EXITING"
#endif  
  
#ifndef WITHOUTMPI
  call MPI_FINALIZE(info)
#endif
  
end subroutine mdl_init
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine mdl_wait(r,g,m,p,mdl,callback)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use call_back, only: call_back_f
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  type(call_back_f),dimension(0:100)::callback

#ifndef WITHOUTMPI

  logical::stop_order_received,order_received
  integer::order_id,order_tag=101,output_tag=203,output_id
  integer,dimension(MPI_STATUS_SIZE)::order_status,output_status
  integer::input_size,output_size,cpu_range,source_cpu,function_id
  integer,dimension(:),allocatable::input_array,output_array
  integer,dimension(1:32)::header
  
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
        cpu_range=header(2)
        input_size=header(3)
        output_size=header(4)
        
        if(input_size>0)then
           allocate(input_array(1:input_size))
           input_array(1:input_size)=mdl%mpi_input_buffer(33:32+input_size)
        endif
        
        if(output_size>0)then
           allocate(output_array(1:output_size))
           output_array=0
        endif
        
        ! Launch the corresponding call-back function
        call callback(function_id)%proc(r,g,m,p,mdl,cpu_range,input_size,output_size,input_array,output_array)
        
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

#endif

end subroutine mdl_wait
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine mdl_send_request(mdl,mdl_function_id,target_cpu,cpu_range,input_size,output_size,input_array)
  use mdl_commons, only: mdl_t
  implicit none
  type(mdl_t)::mdl
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  integer::mdl_function_id
  integer::target_cpu,cpu_range,input_size,output_size
  integer,dimension(1:input_size),optional::input_array

#ifndef WITHOUTMPI
  integer::launch_id,launch_tag=101
  integer,dimension(MPI_STATUS_SIZE)::launch_status
  integer,dimension(1:32)::header=0

  ! Assemble MPI message
  header(1)=mdl_function_id
  header(2)=cpu_range
  header(3)=input_size
  header(4)=output_size
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
#ifndef WITHOUTMPI
  include 'mpif.h'  
  integer::icpu,info
  integer::intex,realdpex,msg_size
  integer,dimension(1:10)::new_type_disp,new_type_type,new_type_length
#endif
  type(mdl_t)::mdl
  
  integer::ncpu
  type(msg_large_realdp)::dummy_large_realdp

  ncpu=mdl%ncpu
  
#ifndef WITHOUTMPI  

  ! Allocate all communication and cache-related variables
  allocate(reply_id(1:ncpu))

#ifdef OLDCACHE
  allocate(reply_interpol(1:ncpu))
  allocate(reply_mg(1:ncpu))
  allocate(reply_flag(1:ncpu))
  allocate(reply_hydro(1:ncpu))
  allocate(reply_poisson(1:ncpu))
  allocate(reply_refine(1:ncpu))
  allocate(send_flush_interpol(1:ncpu))
  allocate(send_flush_mg(1:ncpu))
  allocate(send_flush_flag(1:ncpu))
  allocate(send_flush_hydro(1:ncpu))
  allocate(send_flush_poisson(1:ncpu))
  allocate(send_flush_refine(1:ncpu))
  do icpu=1,ncpu
     send_flush_interpol(icpu)%nflush=0
     send_flush_mg(icpu)%nflush=0
     send_flush_flag(icpu)%nflush=0
     send_flush_hydro(icpu)%nflush=0
     send_flush_poisson(icpu)%nflush=0
     send_flush_refine(icpu)%nflush=0
  end do

  ! Create and commit MPI derived types
  call MPI_TYPE_EXTENT(MPI_INTEGER,intex,info)
  call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,realdpex,info)

  ! New type for int4_msg
  new_type_disp(1)=0
  new_type_disp(2)=2*intex
  new_type_disp(3)=(2+ntilemax*(1+ndim))*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_INTEGER
  new_type_length(1)=2
  new_type_length(2)=ntilemax*(1+ndim)
  new_type_length(3)=ntilemax*(twotondim)
  call MPI_TYPE_STRUCT(3,new_type_length,new_type_disp,new_type_type,new_mpi_int4_msg,info)
  call MPI_TYPE_COMMIT(new_mpi_int4_msg,info)

  ! New type for int4_flush
  new_type_disp(1)=0
  new_type_disp(2)=intex
  new_type_disp(3)=(1+nflushmax*(1+ndim))*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_INTEGER
  new_type_length(1)=1
  new_type_length(2)=nflushmax*(1+ndim)
  new_type_length(3)=nflushmax*(twotondim)
  call MPI_TYPE_STRUCT(3,new_type_length,new_type_disp,new_type_type,new_mpi_int4_flush,info)
  call MPI_TYPE_COMMIT(new_mpi_int4_flush,info)

  ! New type for realdp_msg
  new_type_disp(1)=0
  new_type_disp(2)=2*intex
  new_type_disp(3)=(2+ntilemax*(1+ndim))*intex
  new_type_disp(4)=(2+ntilemax*(1+ndim)+ntilemax*(twotondim))*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_INTEGER
  new_type_type(4)=MPI_DOUBLE_PRECISION
  new_type_length(1)=2
  new_type_length(2)=ntilemax*(1+ndim)
  new_type_length(3)=ntilemax*(twotondim)
  new_type_length(4)=ntilemax*(twotondim*nvar)
  call MPI_TYPE_STRUCT(4,new_type_length,new_type_disp,new_type_type,new_mpi_realdp_msg,info)
  call MPI_TYPE_COMMIT(new_mpi_realdp_msg,info)

  ! New type for realdp_flush
  new_type_disp(1)=0
  new_type_disp(2)=intex
  new_type_disp(3)=(1+nflushmax*(1+ndim))*intex
  new_type_disp(4)=(1+nflushmax*(1+ndim)+nflushmax*(twotondim))*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_INTEGER
  new_type_type(4)=MPI_DOUBLE_PRECISION
  new_type_length(1)=1
  new_type_length(2)=nflushmax*(1+ndim)
  new_type_length(3)=nflushmax*(twotondim)
  new_type_length(4)=nflushmax*(twotondim*nvar)
  call MPI_TYPE_STRUCT(4,new_type_length,new_type_disp,new_type_type,new_mpi_realdp_flush,info)
  call MPI_TYPE_COMMIT(new_mpi_realdp_flush,info)

  ! New type for small_realdp_msg
  new_type_disp(1)=0
  new_type_disp(2)=2*intex
  new_type_disp(3)=(2+ntilemax*(1+ndim))*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_DOUBLE_PRECISION
  new_type_length(1)=2
  new_type_length(2)=ntilemax*(1+ndim)
  new_type_length(3)=ntilemax*(twotondim)
  call MPI_TYPE_STRUCT(3,new_type_length,new_type_disp,new_type_type,new_mpi_small_realdp_msg,info)
  call MPI_TYPE_COMMIT(new_mpi_small_realdp_msg,info)

  ! New type for twin_realdp_msg
  new_type_disp(1)=0
  new_type_disp(2)=2*intex
  new_type_disp(3)=(2+ntilemax*(1+ndim))*intex
  new_type_disp(4)=(2+ntilemax*(1+ndim))*intex+ntilemax*twotondim*realdpex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_DOUBLE_PRECISION
  new_type_type(4)=MPI_DOUBLE_PRECISION
  new_type_length(1)=2
  new_type_length(2)=ntilemax*(1+ndim)
  new_type_length(3)=ntilemax*(twotondim)
  new_type_length(4)=ntilemax*(twotondim)
  call MPI_TYPE_STRUCT(4,new_type_length,new_type_disp,new_type_type,new_mpi_twin_realdp_msg,info)
  call MPI_TYPE_COMMIT(new_mpi_twin_realdp_msg,info)

  ! New type for three_realdp_msg
  new_type_disp(1)=0
  new_type_disp(2)=2*intex
  new_type_disp(3)=(2+ntilemax*(1+ndim))*intex
  new_type_disp(4)=(2+ntilemax*(1+ndim))*intex+ntilemax*twotondim*realdpex
  new_type_disp(5)=(2+ntilemax*(1+ndim))*intex+2*ntilemax*twotondim*realdpex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_DOUBLE_PRECISION
  new_type_type(4)=MPI_DOUBLE_PRECISION
  new_type_type(5)=MPI_DOUBLE_PRECISION
  new_type_length(1)=2
  new_type_length(2)=ntilemax*(1+ndim)
  new_type_length(3)=ntilemax*(twotondim)
  new_type_length(4)=ntilemax*(twotondim)
  new_type_length(5)=ntilemax*(twotondim)
  call MPI_TYPE_STRUCT(5,new_type_length,new_type_disp,new_type_type,new_mpi_three_realdp_msg,info)
  call MPI_TYPE_COMMIT(new_mpi_three_realdp_msg,info)

  ! New type for large_realdp_msg
  new_type_disp(1)=0
  new_type_disp(2)=2*intex
  new_type_disp(3)=new_type_disp(2)+ntilemax*(1+ndim)*intex
  new_type_disp(4)=new_type_disp(3)+ntilemax*twotondim*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_INTEGER
  new_type_length(1)=2
  new_type_length(2)=ntilemax*(1+ndim)
  new_type_length(3)=ntilemax*twotondim
  msg_size=3

#ifdef HYDRO
  msg_size=msg_size+1
  new_type_length(msg_size)=ntilemax*twotondim*nvar
  new_type_type(msg_size)=MPI_DOUBLE_PRECISION
  new_type_disp(msg_size+1)=new_type_disp(msg_size)+ntilemax*twotondim*nvar*realdpex
#endif

#ifdef GRAV
  msg_size=msg_size+1
  new_type_length(msg_size)=ntilemax*twotondim*(ndim+2)
  new_type_type(msg_size)=MPI_DOUBLE_PRECISION
  new_type_disp(msg_size+1)=new_type_disp(msg_size)+ntilemax*twotondim*(ndim+2)*realdpex
#endif

  call MPI_TYPE_STRUCT(msg_size,new_type_length,new_type_disp,new_type_type,new_mpi_large_realdp_msg,info)
  call MPI_TYPE_COMMIT(new_mpi_large_realdp_msg,info)

  ! New type for small_realdp_flush
  new_type_disp(1)=0
  new_type_disp(2)=intex
  new_type_disp(3)=(1+nflushmax*(1+ndim))*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_DOUBLE_PRECISION
  new_type_length(1)=1
  new_type_length(2)=nflushmax*(1+ndim)
  new_type_length(3)=nflushmax*twotondim
  call MPI_TYPE_STRUCT(3,new_type_length,new_type_disp,new_type_type,new_mpi_small_realdp_flush,info)
  call MPI_TYPE_COMMIT(new_mpi_small_realdp_flush,info)

  ! New type for twin_realdp_flush
  new_type_disp(1)=0
  new_type_disp(2)=intex
  new_type_disp(3)=(1+nflushmax*(1+ndim))*intex
  new_type_disp(4)=(1+nflushmax*(1+ndim))*intex+nflushmax*twotondim*realdpex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_DOUBLE_PRECISION
  new_type_type(4)=MPI_DOUBLE_PRECISION
  new_type_length(1)=1
  new_type_length(2)=nflushmax*(1+ndim)
  new_type_length(3)=nflushmax*twotondim
  new_type_length(4)=nflushmax*twotondim
  call MPI_TYPE_STRUCT(4,new_type_length,new_type_disp,new_type_type,new_mpi_twin_realdp_flush,info)
  call MPI_TYPE_COMMIT(new_mpi_twin_realdp_flush,info)

  ! New type for three_realdp_flush
  new_type_disp(1)=0
  new_type_disp(2)=intex
  new_type_disp(3)=(1+nflushmax*(1+ndim))*intex
  new_type_disp(4)=(1+nflushmax*(1+ndim))*intex+nflushmax*twotondim*realdpex
  new_type_disp(5)=(1+nflushmax*(1+ndim))*intex+2*nflushmax*twotondim*realdpex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_DOUBLE_PRECISION
  new_type_type(4)=MPI_DOUBLE_PRECISION
  new_type_type(5)=MPI_DOUBLE_PRECISION
  new_type_length(1)=1
  new_type_length(2)=nflushmax*(1+ndim)
  new_type_length(3)=nflushmax*twotondim
  new_type_length(4)=nflushmax*twotondim
  new_type_length(5)=nflushmax*twotondim
  call MPI_TYPE_STRUCT(5,new_type_length,new_type_disp,new_type_type,new_mpi_three_realdp_flush,info)
  call MPI_TYPE_COMMIT(new_mpi_three_realdp_flush,info)

  ! New type for large_realdp_flush
  new_type_disp(1)=0
  new_type_disp(2)=intex
  new_type_disp(3)=new_type_disp(2)+nflushmax*(1+ndim)*intex
  new_type_disp(4)=new_type_disp(3)+nflushmax*twotondim*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_INTEGER
  new_type_length(1)=1
  new_type_length(2)=nflushmax*(1+ndim)
  new_type_length(3)=nflushmax*twotondim
  msg_size=3

#ifdef HYDRO
  msg_size=msg_size+1
  new_type_type(msg_size)=MPI_DOUBLE_PRECISION
  new_type_length(msg_size)=nflushmax*twotondim*nvar
  new_type_disp(msg_size+1)=new_type_disp(msg_size)+nflushmax*twotondim*nvar*realdpex
#endif

#ifdef GRAV
  msg_size=msg_size+1
  new_type_type(msg_size)=MPI_DOUBLE_PRECISION
  new_type_length(msg_size)=nflushmax*twotondim*(ndim+2)
  new_type_disp(msg_size+1)=new_type_disp(msg_size)+nflushmax*twotondim*(ndim+2)*realdpex
#endif

  call MPI_TYPE_STRUCT(msg_size,new_type_length,new_type_disp,new_type_type,new_mpi_large_realdp_flush,info)
  call MPI_TYPE_COMMIT(new_mpi_large_realdp_flush,info)

  ! New type for request
  new_type_disp(1)=0
  new_type_disp(2)=intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_length(1)=1
  new_type_length(2)=ndim
  call MPI_TYPE_STRUCT(2,new_type_length,new_type_disp,new_type_type,new_mpi_request,info)
  call MPI_TYPE_COMMIT(new_mpi_request,info)

#else
  size_request_array=1+ndim
  size_msg_array=storage_size(dummy_large_realdp)/32
  size_flush_array=1+(1+ndim+size_msg_array)*nflushmax
  size_fetch_array=2+(1+ndim+size_msg_array)*ntilemax

  allocate(recv_request_array(1:size_request_array))
  allocate(send_request_array(1:size_request_array))
  allocate(recv_fetch_array(1:size_fetch_array))
  allocate(send_fetch_array(1:ncpu*size_fetch_array))
  allocate(recv_flush_array(1:size_flush_array))
  allocate(send_flush_array(1:ncpu*size_flush_array))
  
#endif
#endif

end subroutine init_cache










