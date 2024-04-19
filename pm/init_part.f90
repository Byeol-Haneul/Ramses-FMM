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
  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to allocate
  ! particle-based arrays.
  !--------------------------------------------------------------------
  integer::rID
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
  use pm_parameters, only: DM_TYPE
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  !---------------------------------
  ! Allocate DM particle variables
  !---------------------------------
  p%type=DM_TYPE
  allocate(p%xp    (r%npartmax,ndim))
  allocate(p%vp    (r%npartmax,ndim))
  allocate(p%mp    (r%npartmax))
  allocate(p%levelp(r%npartmax))
  allocate(p%idp   (r%npartmax))
  p%nvaralloc=2*ndim+3
#ifdef OUTPUT_PARTICLE_POTENTIAL
  allocate(p%phip  (r%npartmax))
  p%nvaralloc=p%nvaralloc+1
#endif
  ! ALlocate workspace variables
  allocate(p%sortp (r%npartmax))
  allocate(p%workp (r%npartmax))
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
  use pm_parameters, only: STAR_TYPE
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::s
  !-----------------------------------
  ! Allocate star particle variables
  !------------------------------------
  s%type=STAR_TYPE
  allocate(s%xp    (r%nstarmax,ndim))
  allocate(s%vp    (r%nstarmax,ndim))
  allocate(s%mp    (r%nstarmax))
  allocate(s%zp    (r%nstarmax))
  allocate(s%tp    (r%nstarmax))
  allocate(s%levelp(r%nstarmax))
  allocate(s%idp   (r%nstarmax))
  s%nvaralloc=2*ndim+5
#ifdef OUTPUT_PARTICLE_POTENTIAL
  allocate(s%phip  (r%nstarmax))
  s%nvaralloc=p%nvaralloc+1
#endif
  ! Allocate workspace variables
  allocate(s%sortp (r%nstarmax))
  allocate(s%workp (r%nstarmax))
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
subroutine allocate_gas(r,g,gas)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::gas
  !-----------------------------------
  ! Allocate gas sph particle variables
  !------------------------------------
  allocate(gas%xp    (gas%npart,ndim))
  allocate(gas%vp    (gas%npart,ndim))
  allocate(gas%mp    (gas%npart))
  allocate(gas%zp    (gas%npart))
  allocate(gas%up    (gas%npart))
  allocate(gas%levelp(gas%npart))
  gas%nvaralloc=2*ndim+5
  ! Allocate workspace variables
  allocate(gas%sortp (gas%npart))
  allocate(gas%workp (gas%npart))
  ! Allocate pointers to particle levels
  allocate(gas%headp(r%levelmin:r%nlevelmax))
  allocate(gas%tailp(r%levelmin:r%nlevelmax))
  ! No particle just yet
  gas%headp=1
  gas%tailp=0
end subroutine allocate_gas
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_deallocate_gas(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to deallocate
  ! gas sph particle used for Gadget initial conditions.
  !--------------------------------------------------------------------
  integer::rID
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_DEALLOCATE_GAS,pst%iUpper+1)
     call r_deallocate_gas(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call deallocate_gas(pst%s%r,pst%s%g,pst%s%gas)
  endif

end subroutine r_deallocate_gas
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine deallocate_gas(r,g,gas)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::gas
  !-----------------------------------------
  ! Deallocate gas sph particle variables
  !-----------------------------------------
  deallocate(gas%xp)
  deallocate(gas%vp)
  deallocate(gas%mp)
  deallocate(gas%zp)
  deallocate(gas%up)
  deallocate(gas%levelp)
  deallocate(gas%sortp)
  deallocate(gas%workp)
  deallocate(gas%headp)
  deallocate(gas%tailp)
  gas%nvaralloc=0
  gas%npart=0
end subroutine deallocate_gas
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module init_part_module
