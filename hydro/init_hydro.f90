!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_init_hydro(r,g,m,p,mdl,cpu_range,input_size,output_size)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size

  integer::next_range,next_cpu

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_INIT_HYDRO,next_cpu,next_range,input_size,output_size)
     call r_init_hydro(r,g,m,p,mdl,next_range,input_size,output_size)
  else
     call init_hydro(r,m)
  endif

end subroutine r_init_hydro
!###############################################
!###############################################
!###############################################
!###############################################
subroutine init_hydro(r,m)
  use amr_commons, only: run_t,mesh_t
  implicit none
  type(run_t)::r
  type(mesh_t)::m
  
  ! Initialise workspace for hydro kernels
  associate(h => m%hydro_w)
    call h%kernel_1%init(2)
    if(r%nsuperoct>0) call h%kernel_2%init(4)
    if(r%nsuperoct>1) call h%kernel_4%init(8)
    if(r%nsuperoct>2) call h%kernel_8%init(16)
    if(r%nsuperoct>3) call h%kernel_16%init(32)
    if(r%nsuperoct>4) call h%kernel_32%init(64)
  end associate

end subroutine init_hydro
!###############################################
!###############################################
!###############################################
!###############################################




