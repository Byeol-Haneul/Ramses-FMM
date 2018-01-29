!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_init_part(s,cpu_range,input_size,output_size)
  use ramses_commons, only: ramses_t
  use mdl_parameters
  implicit none
  type(ramses_t)::s
  integer::cpu_range,input_size,output_size
  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to allocate
  ! particle-based arrays.
  !--------------------------------------------------------------------
  integer::next_range,next_cpu
  
  next_range=cpu_range/2
  next_cpu=s%g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(s%mdl,MDL_INIT_PART,next_cpu,next_range,input_size,output_size)
     call r_init_part(s,next_range,input_size,output_size)
  else
     call init_part(s%r,s%g,s%p)
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

