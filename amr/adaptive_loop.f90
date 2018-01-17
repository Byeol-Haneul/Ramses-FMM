subroutine adaptive_loop(r,g,m,p,mdl)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only:mdl_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl

  ! Local variables
  integer::ilevel
  real::tt1,tt2
  
  call cpu_time(tt1)

  ! Read run parameters
  call m_read_params(r,g,m,p,mdl)

  ! Initialize grid variables
  call r_init_amr(r,g,m,p,mdl,mdl%ncpu,0,0)

  ! Initialize time variables
  call r_init_time(r,g,m,p,mdl,mdl%ncpu,0,0)

  ! Initialize hydro kernel workspace
  if(r%hydro)call r_init_hydro(r,g,m,p,mdl,mdl%ncpu,0,0)

  ! Initialize particle variables
  if(r%pic)call r_init_part(r,g,m,p,mdl,mdl%ncpu,0,0)

  ! Read initial particle properties from files
  if(r%pic)call m_input_part(r,g,m,p,mdl)
  
  ! Build initial AMR grid
  if(r%nrestart==0)then
     call m_init_refine_basegrid(r,g,m,p,mdl) ! Build coarse grid
     call m_init_refine_adaptive(r,g,m,p,mdl) ! Build adaptive grid
  else
     call m_init_refine_restart(r,g,m,p,mdl) ! Build AMR grid from restart file
  endif

  ! Timing since startup
  call cpu_time(tt2)
  write(*,*)'Time elapsed since startup:',tt2-tt1

  ! Output mesh structure
  write(*,*)'Initial mesh structure'
  do ilevel=r%levelmin,r%nlevelmax
     if(m%noct_tot(ilevel)>0)write(*,999)&
          & ilevel,m%noct_tot(ilevel),m%noct_min(ilevel),m%noct_max(ilevel),m%noct_tot(ilevel)/g%ncpu
  end do
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

  g%nstep_coarse_old=g%nstep_coarse

  write(*,*)'Starting time integration' 

  do ! Main time loop

     call cpu_time(tt1)

     if(r%verbose)write(*,*)'Entering amr_step_coarse'

     g%epot_tot=0.0D0  ! Reset total potential energy
     g%ekin_tot=0.0D0  ! Reset total kinetic energy
     g%mass_tot=0.0D0  ! Reset total mass
     g%eint_tot=0.0D0  ! Reset total internal energy

     ! Call base level
     call m_amr_step(r,g,m,p,mdl,r%levelmin,1)

     ! New coarse time-step
     g%nstep_coarse=g%nstep_coarse+1

     call cpu_time(tt2)
     write(*,*)'Time elapsed since last coarse step:',tt2-tt1
     
  end do

  call r_clean_stop(r,g,m,p,mdl,mdl%ncpu,0,0)

  return

end subroutine adaptive_loop
