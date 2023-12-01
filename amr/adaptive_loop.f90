subroutine adaptive_loop(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use init_amr_module, only: r_init_amr
  use params_module, only: m_read_params
  use init_time_module, only: r_init_time
  use init_hydro_module, only: r_init_hydro
  use init_part_module, only: r_init_part
  use input_part_module, only: m_input_part
  use init_refine_basegrid_module, only: m_init_refine_basegrid
  use init_refine_restart_module, only: m_init_refine_restart
  use amr_step, only: m_amr_step
  use update_time_module, only: getmem, writemem
  implicit none
  type(pst_t)::pst
  logical::done

  ! Local variables
  integer::ilevel
  double precision::tt1,tt2
  real(kind=4)::core_mem

  associate(mdl=>pst%s%mdl,r=>pst%s%r,m=>pst%s%m,g=>pst%s%g)

  tt1 = mdl_wtime(mdl)

  ! Read run parameters
  call m_read_params(pst)

  ! Initialize grid variables
  call r_init_amr(pst)

  ! Initialize time variables
  call r_init_time(pst)

  ! Initialize hydro kernel workspace
  if(r%hydro)call r_init_hydro(pst)

  ! Initialize particle variables
  if(r%pic)call r_init_part(pst)

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
  tt2 = mdl_wtime(mdl)
  print '(A,F14.7)',' Time elapsed since startup:',tt2-tt1

  ! Output mesh structure
  write(*,*)'Initial mesh structure'
  do ilevel=r%levelmin,r%nlevelmax
     if(m%noct_tot(ilevel)>0)write(*,999)&
          & ilevel,m%noct_tot(ilevel),m%noct_min(ilevel),m%noct_max(ilevel),m%noct_tot(ilevel)/mdl_threads(mdl)
  end do
999 format(' Level ',I2,' has ',I11,' grids (',3(I8,','),')')

  g%nstep_coarse_old=g%nstep_coarse

  write(*,*)'Starting time integration' 

  done = .false.
  do while(.not.done) ! Main time loop

     tt1 = mdl_wtime(mdl)

     if(r%verbose)write(*,*)'Entering amr_step_coarse'

     g%epot_tot=0.0D0  ! Reset total potential energy
     g%ekin_tot=0.0D0  ! Reset total kinetic energy
     g%mass_tot=0.0D0  ! Reset total mass
     g%eint_tot=0.0D0  ! Reset total internal energy

     ! Call base level
     call m_amr_step(pst,r%levelmin,1,done)

     ! New coarse time-step
     g%nstep_coarse=g%nstep_coarse+1

     tt2 = mdl_wtime(mdl)
     if(mod(g%nstep_coarse,r%ncontrol)==0)then
        if(.not. done)print '(A,F14.7)',' Time elapsed since last coarse step:',tt2-tt1
     endif

     call getmem(core_mem)
     call writemem(core_mem)
     
  end do

  call m_output_timer(pst,.false.,'dummy')
  
  return

  end associate
  
end subroutine adaptive_loop
