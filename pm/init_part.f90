!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_init_part(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::rID

  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to allocate
  ! particle-based arrays.
  !--------------------------------------------------------------------
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INIT_PART,pst%iUpper+1,input_size,output_size,input_array)
     call r_init_part(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     call init_part(pst%s%r,pst%s%g,pst%s%p)
  endif

end subroutine r_init_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_part(r,g,p)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  !------------------------------
  ! Allocate particle variables
  !------------------------------
  allocate(p%xp    (r%npartmax,ndim))
  allocate(p%vp    (r%npartmax,ndim))
  allocate(p%mp    (r%npartmax))
  allocate(p%levelp(r%npartmax))
  allocate(p%idp   (r%npartmax))
  allocate(p%sortp (r%npartmax))
  allocate(p%workp (r%npartmax))
#ifdef OUTPUT_PARTICLE_POTENTIAL
  allocate(p%phip  (r%npartmax))
#endif
  ! Allocate pointers to particle levels
  allocate(p%headp(r%levelmin:r%nlevelmax))
  allocate(p%tailp(r%levelmin:r%nlevelmax))
  ! No particle just yet
  p%headp=1
  p%tailp=0
end subroutine init_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

