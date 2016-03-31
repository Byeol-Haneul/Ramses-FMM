subroutine adaptive_loop
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,info
  real(kind=8)::tt1,tt2
  real(kind=4)::real_mem,real_mem_tot

#ifndef WITHOUTMPI
  tt1=MPI_WTIME(info)
#endif

  call init_amr                ! Initialize grid variables
  call init_time               ! Initialize time variables
  if(hydro)call init_hydro     ! Initialize hydro variables
  if(poisson)call init_poisson ! Initialize poisson variables
  if(pic)call init_part        ! Initialize particle variables
  if(nrestart==0)then
     call init_refine_basegrid ! Build initial coarsest grid
     call init_refine_adaptive ! Build initial adaptive grid
  else
     call init_refine_restart  ! Build AMR grid from restart file
  endif

#ifndef WITHOUTMPI
  tt2=MPI_WTIME(info)
  call getmem(real_mem)
  call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
  if(myid==1)then
     write(*,*)'Time elapsed since startup:',tt2-tt1
     call writemem(real_mem_tot)
  endif
#endif

  if(myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=levelmin,nlevelmax
        if(noct_tot(ilevel)>0)write(*,999)ilevel,noct_tot(ilevel),noct_min(ilevel),noct_max(ilevel),noct_tot(ilevel)/ncpu
     end do
  end if
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

  nstep_coarse_old=nstep_coarse

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

     ! Call base level
     call amr_step(levelmin,1)

     ! New coarse time-step
     nstep_coarse=nstep_coarse+1

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
