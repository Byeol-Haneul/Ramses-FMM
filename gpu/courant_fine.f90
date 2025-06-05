module courant_fine_module

type :: out_courant_fine_t
  real(kind=8)::mass,ekin,eint,emag,dt
end type out_courant_fine_t

contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_courant_fine(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use gpu_runner, only: gpu_cmpdt
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(out_courant_fine_t)::output, next_output

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_COURANT_FINE,pst%iUpper+1,input_size,output_size,ilevel)
     call r_courant_fine(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
     output%ekin=output%ekin+next_output%ekin
     output%eint=output%eint+next_output%eint
     output%emag=output%emag+next_output%emag
     output%dt=MIN(output%dt,next_output%dt)
  else
     call gpu_cmpdt(pst%s,ilevel,output%mass,output%ekin,output%eint,output%emag,output%dt)
  endif

end subroutine r_courant_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module courant_fine_module
