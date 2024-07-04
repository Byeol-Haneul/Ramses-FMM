module newdt_fine_module

type :: out_newdt_part_t
  real(kind=8)::ekin,vmax
end type out_newdt_part_t

type :: in_broadcast_dt_t
  integer::ilevel
  real(kind=8)::dtnew,dtold
end type in_broadcast_dt_t

contains
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine m_newdt_fine(pst,ilevel)
  use amr_parameters, only: dp,nvector
  use ramses_commons, only: pst_t
  use courant_fine_module, only: r_courant_fine, out_courant_fine_t
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
  real(kind=8)::ekin,vmax
  type(out_courant_fine_t)::out_courant_fine
  type(out_newdt_part_t)::out_newdt_part
  type(in_broadcast_dt_t)::in_broadcast_dt
  
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
     if(r%verbose)write(*,'("   Entering newdt_part for level ",I2)')ilevel
     call r_newdt_part(pst,ilevel,1,out_newdt_part,storage_size(out_newdt_part)/32)
     ekin=out_newdt_part%ekin
     g%ekin_tot=g%ekin_tot+ekin
     vmax=out_newdt_part%vmax
     if(vmax>0.0d0)then
        g%dtnew(ilevel)=MIN(g%dtnew(ilevel),r%courant_factor*dx/vmax)
     endif
  endif

  ! Hydro-based Courant condition
  if(r%hydro)then
     if(r%verbose)write(*,'("   Entering newdt_hydro for level ",I2)')ilevel
     call r_courant_fine(pst,ilevel,1,out_courant_fine,storage_size(out_courant_fine)/32)
     g%mass_tot=g%mass_tot+out_courant_fine%mass
     g%ekin_tot=g%ekin_tot+out_courant_fine%ekin
     g%eint_tot=g%eint_tot+out_courant_fine%eint
     g%emag_tot=g%emag_tot+out_courant_fine%emag
     g%dtnew(ilevel)=MIN(g%dtnew(ilevel),out_courant_fine%dt)
  endif  
  
  ! Adaptive time step condition
  if(ilevel>r%levelmin)then
     g%dtnew(ilevel)=MIN(g%dtnew(ilevel-1)/real(r%nsubcycle(ilevel-1)),g%dtnew(ilevel))
  end if
  if(r%verbose)write(*,'("   New time step for level ",I2," is ",1PE12.5)')ilevel,g%dtnew(ilevel)

  ! Broadcast new and old time steps to all CPUs
  in_broadcast_dt%ilevel=ilevel
  in_broadcast_dt%dtnew=g%dtnew(ilevel)
  in_broadcast_dt%dtold=g%dtold(ilevel)
  call r_broadcast_dt(pst,in_broadcast_dt,storage_size(in_broadcast_dt)/32)

  end associate
  
end subroutine m_newdt_fine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_newdt_part(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(out_newdt_part_t)::output,next_output
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_NEWDT_PART,pst%iUpper+1,input_size,output_size,ilevel)
     call r_newdt_part(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%ekin=output%ekin+next_output%ekin
     output%vmax=MAX(output%vmax,next_output%vmax)
  else
     output%vmax=0.0d0
     output%ekin=0.0d0
     call newdt_part(pst%s%r,pst%s%g,pst%s%p,ilevel,output%ekin,output%vmax)
     if(pst%s%r%star)then
        call newdt_part(pst%s%r,pst%s%g,pst%s%star,ilevel,output%ekin,output%vmax)
     endif
     if(pst%s%r%sink)then
        call newdt_part(pst%s%r,pst%s%g,pst%s%sink,ilevel,output%ekin,output%vmax)
     endif
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
  do idim = 1, ndim
     do ipart = p%headp(ilevel), p%tailp(ilevel)
        vmax = MAX(vmax, ABS(p%vp(ipart, idim)))
     end do
  end do
  
  ! Compute kinetic energy
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
recursive subroutine r_broadcast_dt(pst,input,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(in_broadcast_dt_t)::input

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BROADCAST_DT,pst%iUpper+1,input_size,0,input)
     call r_broadcast_dt(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     pst%s%g%dtnew(input%ilevel)=input%dtnew
     pst%s%g%dtold(input%ilevel)=input%dtold
  endif

end subroutine r_broadcast_dt
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
end module newdt_fine_module
