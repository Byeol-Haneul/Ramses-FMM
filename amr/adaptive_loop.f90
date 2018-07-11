subroutine adaptive_loop(pst)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst

  ! Local variables
  integer::ilevel,dummy
  real::tt1,tt2

  associate(mdl=>pst%s%mdl,r=>pst%s%r,m=>pst%s%m,g=>pst%s%g)
  
  call cpu_time(tt1)

  ! Read run parameters
  call m_read_params(pst)

  ! Initialize grid variables
  call r_init_amr(pst,dummy,0,dummy,0)

  ! Initialize time variables
  call r_init_time(pst,dummy,0,dummy,0)

  ! Initialize hydro kernel workspace
  if(r%hydro)call r_init_hydro(pst,dummy,0,dummy,0)

  ! Initialize particle variables
  if(r%pic)call r_init_part(pst,dummy,0,dummy,0)

  ! Read initial particle properties from files
  if(r%pic)call m_input_part(pst)
  
  ! Build initial AMR grid
  if(r%nrestart==0)then
     call m_init_refine_basegrid(pst) ! Build coarse grid
     call m_init_refine_adaptive(pst) ! Build adaptive grid
  else
     call m_init_refine_restart(pst) ! Build AMR grid from restart file
  endif

  ! Timing since startup
  call cpu_time(tt2)
  write(*,*)'Time elapsed since startup:',tt2-tt1

  ! Output mesh structure
  write(*,*)'Initial mesh structure'
  do ilevel=r%levelmin,r%nlevelmax
     if(m%noct_tot(ilevel)>0)write(*,999)&
          & ilevel,m%noct_tot(ilevel),m%noct_min(ilevel),m%noct_max(ilevel),m%noct_tot(ilevel)/mdl%ncpu
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
     call m_amr_step(pst,r%levelmin,1)

     ! New coarse time-step
     g%nstep_coarse=g%nstep_coarse+1

     call cpu_time(tt2)
     write(*,*)'Time elapsed since last coarse step:',tt2-tt1
     
  end do

  call r_clean_stop(pst,dummy,0,dummy,0)

  return

  end associate
  
end subroutine adaptive_loop
