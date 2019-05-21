module init_hydro_module

contains
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_init_hydro(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INIT_HYDRO,pst%iUpper+1,input_size,output_size,input_array)
     call r_init_hydro(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     call init_hydro(pst%s%r,pst%s%m)
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
end module init_hydro_module