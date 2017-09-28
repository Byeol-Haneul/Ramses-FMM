subroutine adaptive_loop(r,g,m,p)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  real(kind=4)::real_mem_tot
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p

  ! Local variables
  integer::ilevel
  real(kind=4)::real_mem

#ifndef WITHOUTMPI
  tt1=MPI_WTIME(info)
#endif

  ! Initialize grid variables
  call init_amr(r,g,m)

  ! Initialize software cache
  call init_cache(r,g)

  ! Initialize time variables
  call init_time(r,g)

  ! Initialize hydro kernel workspace
  if(r%hydro)call init_hydro(r)

  ! Initialize particle variables from files
  if(r%pic)call init_part_file(r,g,p)

  ! Build initial AMR grid
  if(r%nrestart==0)then
     call init_refine_basegrid(r,g,m,p) ! Build coarsest grid
     call init_refine_adaptive(r,g,m,p) ! Build adaptive grid
  else
     call init_refine_restart(r,g,m) ! Build AMR grid from restart file
  endif

  ! Initialize particle from the AMR grid
  if(r%pic)call init_part_grid(r,g,m,p)

  ! Timing since startup
#ifndef WITHOUTMPI
  tt2=MPI_WTIME(info)
  call getmem(real_mem)
  call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
  if(g%myid==1)then
     write(*,*)'Time elapsed since startup:',tt2-tt1
     call writemem(real_mem_tot)
  endif
#endif

  ! Output mesh structure
  if(g%myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=r%levelmin,r%nlevelmax
        if(m%noct_tot(ilevel)>0)write(*,999)ilevel,m%noct_tot(ilevel),m%noct_min(ilevel),m%noct_max(ilevel),m%noct_tot(ilevel)/g%ncpu
     end do
  end if
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

  g%nstep_coarse_old=g%nstep_coarse

  if(g%myid==1)write(*,*)'Starting time integration' 

  do ! Main time loop

#ifndef WITHOUTMPI
     tt1=MPI_WTIME(info)
#endif

     if(r%verbose)write(*,*)'Entering amr_step_coarse'

     g%epot_tot=0.0D0  ! Reset total potential energy
     g%ekin_tot=0.0D0  ! Reset total kinetic energy
     g%mass_tot=0.0D0  ! Reset total mass
     g%eint_tot=0.0D0  ! Reset total internal energy

     ! Call base level
     call amr_step(r,g,m,p,r%levelmin,1)

     ! New coarse time-step
     g%nstep_coarse=g%nstep_coarse+1
     
#ifndef WITHOUTMPI
     tt2=MPI_WTIME(info)
     if(mod(g%nstep_coarse,r%ncontrol)==0)then
        call getmem(real_mem)
        call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        if(g%myid==1)then
           write(*,*)'Time elapsed since last coarse step:',tt2-tt1
           call writemem(real_mem_tot)
        endif
     endif
#else
     if(mod(g%nstep_coarse,r%ncontrol)==0)then
        call getmem(real_mem)
        call writemem(real_mem)
     endif
#endif

  end do

end subroutine adaptive_loop
