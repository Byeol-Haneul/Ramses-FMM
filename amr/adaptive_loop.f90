subroutine adaptive_loop(r,g,m,p)
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p

  ! Local variables
  integer::ilevel,info
  real(kind=8)::tt1,tt2
  real(kind=4)::real_mem,real_mem_tot

#ifndef WITHOUTMPI
  tt1=MPI_WTIME(info)
#endif

  ! Initialize grid variables
  call init_amr
  call init_amr_2(r,g,m) ! Thread safe version

  ! Initialize software cache
  call init_cache

  ! Initialize time variables
  call init_time
  call init_time_2(r,g) ! Thread safe version

  ! Initialize hydro kernel workspace
  if(hydro)call init_hydro

  ! Initialize particle variables from files
  if(pic)call init_part_file
  if(r%pic)call init_part_file_2(r,g,p) ! Thread safe version

  ! Build initial AMR grid
  if(nrestart==0)then
     call init_refine_basegrid ! Build coarsest grid
     call init_refine_adaptive ! Build adaptive grid
  else
     call init_refine_restart ! Build AMR grid from restart file
  endif
  if(r%nrestart==0)then ! Thread safe version
     call init_refine_basegrid_2(r,g,m,p) ! Build coarsest grid
     call init_refine_adaptive_2(r,g,m,p) ! Build adaptive grid
  else
     call init_refine_restart_2(r,g,m) ! Build AMR grid from restart file
  endif

  ! Initialize particle from the AMR grid
  if(pic)call init_part_grid
  if(r%pic)call init_part_grid_2(r,g,m,p) ! Thread safe version

  ! Timing since startup
#ifndef WITHOUTMPI
  tt2=MPI_WTIME(info)
  call getmem(real_mem)
  call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
  if(myid==1)then
     write(*,*)'Time elapsed since startup:',tt2-tt1
     call writemem(real_mem_tot)
  endif
#endif

  ! Output mesh structure
  if(myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=levelmin,nlevelmax
        if(noct_tot(ilevel)>0)write(*,999)ilevel,noct_tot(ilevel),noct_min(ilevel),noct_max(ilevel),noct_tot(ilevel)/ncpu
     end do
  end if
  ! Thread safe version
  if(g%myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=r%levelmin,r%nlevelmax
        if(m%noct_tot(ilevel)>0)write(*,999)ilevel,m%noct_tot(ilevel),m%noct_min(ilevel),m%noct_max(ilevel),m%noct_tot(ilevel)/g%ncpu
     end do
  end if
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

  nstep_coarse_old=nstep_coarse
  g%nstep_coarse_old=g%nstep_coarse

  if(myid==1)write(*,*)'Starting time integration' 

  do ! Main time loop

#ifndef WITHOUTMPI
     tt1=MPI_WTIME(info)
#endif

     if(verbose)write(*,*)'Entering amr_step_coarse'

     epot_tot=0.0D0  ! Reset total potential energy
     ekin_tot=0.0D0  ! Reset total kinetic energy
     mass_tot=0.0D0  ! Reset total mass
     eint_tot=0.0D0  ! Reset total internal energy

     g%epot_tot=0.0D0  ! Reset total potential energy
     g%ekin_tot=0.0D0  ! Reset total kinetic energy
     g%mass_tot=0.0D0  ! Reset total mass
     g%eint_tot=0.0D0  ! Reset total internal energy

     ! Call base level
     call amr_step(levelmin,1)
     call amr_step_2(r,g,m,p,r%levelmin,1)

     ! New coarse time-step
     nstep_coarse=nstep_coarse+1
     g%nstep_coarse=g%nstep_coarse+1
     
#ifndef WITHOUTMPI
     tt2=MPI_WTIME(info)
     if(mod(nstep_coarse,ncontrol)==0)then
        call getmem(real_mem)
        call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        if(myid==1)then
           write(*,*)'Time elapsed since last coarse step:',tt2-tt1
           call writemem(real_mem_tot)
        endif
     endif
#else
     if(mod(nstep_coarse,ncontrol)==0)then
        call getmem(real_mem)
        call writemem(real_mem)
     endif
#endif

  end do

end subroutine adaptive_loop
