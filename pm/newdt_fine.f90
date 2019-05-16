!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine m_newdt_fine(pst,ilevel)
  use amr_parameters, only: dp,nvector
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !-----------------------------------------------------------
  ! This routine compute the time step using 3 constraints:
  ! 1- a Courant-type condition using particle velocity
  ! 2- the gravity free-fall time
  ! 3- 10% maximum variation for aexp 
  ! This routine also compute the particle kinetic energy.
  !-----------------------------------------------------------
  real(dp)::dx,tff,fourpi,threepi2
  real(kind=8)::mass,ekin,eint,dt,vmax
  integer,dimension(1:4)::part_array,dummy
  integer,dimension(1:5)::input_array
  integer,dimension(1:8)::output_array
  
  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)

  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,'("   Entering newdt_fine for level ",I2)')ilevel

  ! Save old time step
  g%dtold(ilevel)=g%dtnew(ilevel)

  ! Compute local cell spacing
  dx=r%boxlen/2.0d0**ilevel

  ! Maximum time step
  g%dtnew(ilevel)=dx/r%smallc

  ! Gravity-based Courant condition
  if(r%poisson.and.r%gravity_type<=0)then
     fourpi=4.0d0*ACOS(-1.0d0)
     if(r%cosmo)fourpi=1.5d0*g%omega_m*g%aexp
     threepi2=3.0d0*ACOS(-1.0d0)**2
     tff=sqrt(threepi2/8./fourpi/g%rho_max(ilevel))
     g%dtnew(ilevel)=MIN(g%dtnew(ilevel),r%courant_factor*tff)
  end if

  ! Cosmic-expansion-based Courant condition
  if(r%cosmo)then
     g%dtnew(ilevel)=MIN(g%dtnew(ilevel),0.1/g%hexp)
  end if

  ! Particle-based Courant condition
  if(r%pic)then
     call r_newdt_part(pst,ilevel,1,part_array,4)
     ekin=transfer(part_array(1:2),ekin)
     vmax=transfer(part_array(3:4),vmax)
     dt=r%courant_factor * dx/vmax
     g%ekin_tot=g%ekin_tot+ekin
     g%dtnew(ilevel)=MIN(g%dtnew(ilevel),dt)
  endif

  ! Hydro-based Courant condition
  if(r%hydro)then
     call r_courant_fine(pst,ilevel,1,output_array,8)
     mass=transfer(output_array(1:2),mass)
     ekin=transfer(output_array(3:4),ekin)
     eint=transfer(output_array(5:6),eint)
     dt=transfer(output_array(7:8),dt)
     g%mass_tot=g%mass_tot+mass
     g%ekin_tot=g%ekin_tot+ekin
     g%eint_tot=g%eint_tot+eint
     g%dtnew(ilevel)=MIN(g%dtnew(ilevel),dt)
  endif  
  
  ! Adaptive time step condition
  if(ilevel>r%levelmin)then
     g%dtnew(ilevel)=MIN(g%dtnew(ilevel-1)/real(r%nsubcycle(ilevel-1)),g%dtnew(ilevel))
  end if
  if(r%verbose)write(*,'("   New time step for level ",I2," is ",1PE12.5)')ilevel,g%dtnew(ilevel)

  ! Broadcast new and old time steps to all CPUs
  input_array(1)=ilevel
  input_array(2:3)=transfer(g%dtnew(ilevel),input_array)
  input_array(4:5)=transfer(g%dtold(ilevel),input_array)
  call r_broadcast_dt(pst,input_array,5,dummy,0)

  end associate
  
end subroutine m_newdt_fine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_newdt_part(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer,dimension(1:output_size)::next_output_array
  integer::ilevel
  real(kind=8)::ekin,vmax
  real(kind=8)::next_ekin,next_vmax

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_NEWDT_PART,pst%iUpper+1,input_size,output_size,input_array)
     call r_newdt_part(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size,next_output_array)
     ekin=transfer(output_array(1:2),ekin)
     vmax=transfer(output_array(3:4),vmax)
     next_ekin=transfer(next_output_array(1:2),next_ekin)
     next_vmax=transfer(next_output_array(3:4),next_vmax)
     ekin=ekin+next_ekin
     vmax=MAX(vmax,next_vmax)
     output_array(1:2)=transfer(ekin,output_array)
     output_array(3:4)=transfer(vmax,output_array)
  else
     ilevel=input_array(1)
     call newdt_part(pst%s%r,pst%s%g,pst%s%p,ilevel,ekin,vmax)
     output_array(1:2)=transfer(ekin,output_array)
     output_array(3:4)=transfer(vmax,output_array)
  endif

end subroutine r_newdt_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine newdt_part(r,g,p,ilevel,ekin,vmax)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  integer::ilevel
  real(kind=8)::ekin,vmax

  integer::ipart,idim

  ! Compute maximum particle velocity
  vmax = 0.0d0
  do idim = 1, ndim
     do ipart = p%headp(ilevel), p%tailp(ilevel)
        vmax = MAX(vmax, ABS(p%vp(ipart, idim)))
     end do
  end do
  
  ! Compute kinetic energy
  ekin = 0.0d0
  do idim = 1, ndim
     do ipart = p%headp(ilevel), p%tailp(ilevel)
        ekin = ekin + 0.5D0 * p%mp(ipart) * p%vp(ipart, idim)**2
     end do
  end do
    
end subroutine newdt_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_broadcast_dt(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel
  real(kind=8)::dt

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_BROADCAST_DT,pst%iUpper+1,input_size,output_size,input_array)
     call r_broadcast_dt(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     ilevel=input_array(1)
     pst%s%g%dtnew(ilevel)=transfer(input_array(2:3),dt)
     pst%s%g%dtold(ilevel)=transfer(input_array(4:5),dt)
  endif

end subroutine r_broadcast_dt
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################



