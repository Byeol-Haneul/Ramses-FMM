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
        call init_star(pst%s%r,pst%s%g,pst%s%star)
     end if
     if(pst%s%r%sink)then
        call init_sink(pst%s%r,pst%s%g,pst%s%sink)
     end if
     if(pst%s%r%tree)then
        call init_tree(pst%s%r,pst%s%g,pst%s%tree)
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
subroutine init_star(r,g,p)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t,global_t
  use pm_parameters, only: STAR_TYPE
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  !-----------------------------------
  ! Allocate star particle variables
  !------------------------------------
  p%type=STAR_TYPE
  allocate(p%xp    (r%nstarmax,ndim))
  allocate(p%vp    (r%nstarmax,ndim))
  allocate(p%mp    (r%nstarmax))
  allocate(p%zp    (r%nstarmax))
  allocate(p%tp    (r%nstarmax))
  allocate(p%levelp(r%nstarmax))
  allocate(p%idp   (r%nstarmax))
  p%nvaralloc=2*ndim+5
#ifdef OUTPUT_PARTICLE_POTENTIAL
  allocate(p%phip  (r%nstarmax))
  p%nvaralloc=p%nvaralloc+1
#endif
  ! Allocate workspace variables
  allocate(p%sortp (r%nstarmax))
  allocate(p%workp (r%nstarmax))
  ! Allocate pointers to particle levels
  allocate(p%headp(r%levelmin:r%nlevelmax))
  allocate(p%tailp(r%levelmin:r%nlevelmax))
  ! No particle just yet
  p%headp=1
  p%tailp=0
end subroutine init_star
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_sink(r,g,p)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t,global_t
  use pm_parameters, only: SINK_TYPE
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  !-----------------------------------
  ! Allocate sink particle variables
  !------------------------------------
  p%type=SINK_TYPE
  allocate(p%xp    (r%nsinkmax,ndim))
  allocate(p%vp    (r%nsinkmax,ndim))
  allocate(p%fp    (r%nsinkmax,ndim))
  allocate(p%jp    (r%nsinkmax,ndim))
  allocate(p%mp    (r%nsinkmax))
  allocate(p%tp    (r%nsinkmax))
  allocate(p%levelp(r%nsinkmax))
  allocate(p%idp   (r%nsinkmax))
  p%nvaralloc=3*ndim+4
#ifdef OUTPUT_PARTICLE_POTENTIAL
  allocate(p%phip  (r%nsinkmax))
  p%nvaralloc=p%nvaralloc+1
#endif
  ! Allocate workspace variables
  allocate(p%sortp (r%nsinkmax))
  allocate(p%workp (r%nsinkmax))
  ! Allocate pointers to particle levels
  allocate(p%headp(r%levelmin:r%nlevelmax))
  allocate(p%tailp(r%levelmin:r%nlevelmax))
  ! No particle just yet
  p%headp=1
  p%tailp=0
end subroutine init_sink
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_tree(r,g,p)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t,global_t
  use pm_parameters, only: TREE_TYPE
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  !-----------------------------------
  ! Allocate tree particle variables
  !------------------------------------
  p%type=TREE_TYPE
  allocate(p%xp    (r%ntreemax,ndim))
  allocate(p%vp    (r%ntreemax,ndim))
  allocate(p%mp    (r%ntreemax))
  allocate(p%tp    (r%ntreemax))
  allocate(p%tm    (r%ntreemax))
  allocate(p%levelp(r%ntreemax))
  allocate(p%idp   (r%ntreemax))
  allocate(p%idm   (r%ntreemax))
  p%nvaralloc=2*ndim+6
#ifdef OUTPUT_PARTICLE_POTENTIAL
  allocate(p%phip  (r%ntreemax))
  p%nvaralloc=p%nvaralloc+1
#endif
  ! Allocate workspace variables
  allocate(p%sortp (r%ntreemax))
  allocate(p%workp (r%ntreemax))
  ! Allocate pointers to particle levels
  allocate(p%headp(r%levelmin:r%nlevelmax))
  allocate(p%tailp(r%levelmin:r%nlevelmax))
  ! No particle just yet
  p%headp=1
  p%tailp=0
end subroutine init_tree
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine allocate_gas(r,g,p)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  !-----------------------------------
  ! Allocate gas sph particle variables
  !------------------------------------
  allocate(p%xp    (p%npart,ndim))
  allocate(p%vp    (p%npart,ndim))
  allocate(p%mp    (p%npart))
  allocate(p%zp    (p%npart))
  allocate(p%up    (p%npart))
  allocate(p%levelp(p%npart))
  p%nvaralloc=2*ndim+5
  ! Allocate workspace variables
  allocate(p%sortp (p%npart))
  allocate(p%workp (p%npart))
  ! Allocate pointers to particle levels
  allocate(p%headp(r%levelmin:r%nlevelmax))
  allocate(p%tailp(r%levelmin:r%nlevelmax))
  ! No particle just yet
  p%headp=1
  p%tailp=0
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
subroutine deallocate_gas(r,g,p)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  !-----------------------------------------
  ! Deallocate gas sph particle variables
  !-----------------------------------------
  deallocate(p%xp)
  deallocate(p%vp)
  deallocate(p%mp)
  deallocate(p%zp)
  deallocate(p%up)
  deallocate(p%levelp)
  deallocate(p%sortp)
  deallocate(p%workp)
  deallocate(p%headp)
  deallocate(p%tailp)
  p%nvaralloc=0
  p%npart=0
end subroutine deallocate_gas
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module init_part_module
