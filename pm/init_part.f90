module init_part_module

contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_init_part(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID

  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to allocate
  ! particle-based arrays.
  !--------------------------------------------------------------------
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INIT_PART,pst%iUpper+1)
     call r_init_part(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call init_part(pst%s%r,pst%s%g,pst%s%p)
     if(pst%s%r%star)then
        call init_star(pst%s%r,pst%s%g,pst%s%s)
     end if
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
  !---------------------------------
  ! Allocate DM particle variables
  !---------------------------------
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
subroutine init_star(r,g,s)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::s
  !-----------------------------------
  ! Allocate star particle variables
  !------------------------------------
  allocate(s%xp    (r%nstarmax,ndim))
  allocate(s%vp    (r%nstarmax,ndim))
  allocate(s%mp    (r%nstarmax))
  allocate(s%zp    (r%nstarmax))
  allocate(s%tp    (r%nstarmax))
  allocate(s%levelp(r%nstarmax))
  allocate(s%idp   (r%nstarmax))
  allocate(s%sortp (r%nstarmax))
  allocate(s%workp (r%nstarmax))
#ifdef OUTPUT_PARTICLE_POTENTIAL
  allocate(s%phip  (r%nstarmax))
#endif
  ! Allocate pointers to particle levels
  allocate(s%headp(r%levelmin:r%nlevelmax))
  allocate(s%tailp(r%levelmin:r%nlevelmax))
  ! No particle just yet
  s%headp=1
  s%tailp=0
end subroutine init_star
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_gas_part(r,g,p,ngas)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  !------------------------------------------
  ! Allocate gadget gas particle variables
  !------------------------------------------
  allocate(p%xp    (r%ngasmax,ndim))
  allocate(p%vp    (r%ngasmax,ndim))
  allocate(p%mp    (r%ngasmax))
  allocate(p%zp    (r%ngasmax))
  allocate(p%up    (r%ngasmax))
  allocate(p%levelp(r%ngasmax))
  allocate(p%idp   (r%ngasmax))
  allocate(p%sortp (r%ngasmax))
  allocate(p%workp (r%ngasmax))

  ! Allocate pointers to particle levels
  allocate(p%headp(r%levelmin:r%nlevelmax))
  allocate(p%tailp(r%levelmin:r%nlevelmax))
  ! No particle just yet
  p%headp=1
  p%tailp=0
end subroutine init_gas_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module init_part_module
