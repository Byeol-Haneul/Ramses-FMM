subroutine adaptive_loop(s)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s

  ! Local variables
  integer::ilevel,ncpu
  real::tt1,tt2
  
  call cpu_time(tt1)

  ncpu=s%mdl%ncpu
  
  ! Read run parameters
  call m_read_params(s)

  ! Initialize grid variables
  call r_init_amr(s,ncpu,0,0)

  ! Initialize time variables
  call r_init_time(s,ncpu,0,0)

  ! Initialize hydro kernel workspace
  if(s%r%hydro)call r_init_hydro(s,ncpu,0,0)

  ! Initialize particle variables
  if(s%r%pic)call r_init_part(s,ncpu,0,0)

  ! Read initial particle properties from files
  if(s%r%pic)call m_input_part(s)
  
  ! Build initial AMR grid
  if(s%r%nrestart==0)then
     call m_init_refine_basegrid(s) ! Build coarse grid
     call m_init_refine_adaptive(s) ! Build adaptive grid
  else
     call m_init_refine_restart(s) ! Build AMR grid from restart file
  endif

  ! Timing since startup
  call cpu_time(tt2)
  write(*,*)'Time elapsed since startup:',tt2-tt1

  ! Output mesh structure
  write(*,*)'Initial mesh structure'
  do ilevel=s%r%levelmin,s%r%nlevelmax
     if(s%m%noct_tot(ilevel)>0)write(*,999)&
          & ilevel,s%m%noct_tot(ilevel),s%m%noct_min(ilevel),s%m%noct_max(ilevel),s%m%noct_tot(ilevel)/ncpu
  end do
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

  s%g%nstep_coarse_old=s%g%nstep_coarse

  write(*,*)'Starting time integration' 

  do ! Main time loop

     call cpu_time(tt1)

     if(s%r%verbose)write(*,*)'Entering amr_step_coarse'

     s%g%epot_tot=0.0D0  ! Reset total potential energy
     s%g%ekin_tot=0.0D0  ! Reset total kinetic energy
     s%g%mass_tot=0.0D0  ! Reset total mass
     s%g%eint_tot=0.0D0  ! Reset total internal energy

     ! Call base level
     call m_amr_step(s,s%r%levelmin,1)

     ! New coarse time-step
     s%g%nstep_coarse=s%g%nstep_coarse+1

     call cpu_time(tt2)
     write(*,*)'Time elapsed since last coarse step:',tt2-tt1
     
  end do

  call r_clean_stop(s,ncpu,0,0)

  return

end subroutine adaptive_loop
