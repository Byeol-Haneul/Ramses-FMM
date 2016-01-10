subroutine adaptive_loop
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,idim,ivar,info
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
  endif

#ifndef WITHOUTMPI
  tt2=MPI_WTIME(info)
  if(myid==1)write(*,*)'Time elapsed since startup:',tt2-tt1
#endif

  if(myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=levelmin,nlevelmax
        if(noct(ilevel)>0)write(*,999)ilevel,noct(ilevel),head(ilevel),tail(ilevel)
     end do
  end if
999 format(' Level ',I2,' has ',I10,' grids (',2(I8,','),')')

  nstep_coarse_old=nstep_coarse

  if(myid==1)write(*,*)'Starting time integration' 

#ifdef TOTO

  do ! Main time loop

#ifndef WITHOUTMPI
     tt1=MPI_WTIME(info)
#endif

     if(verbose)write(*,*)'Entering amr_step_coarse'

     epot_tot=0.0D0  ! Reset total potential energy
     ekin_tot=0.0D0  ! Reset total kinetic energy
     mass_tot=0.0D0  ! Reset total mass
     eint_tot=0.0D0  ! Reset total internal energy

     ! Make new refinements
     if(levelmin.lt.nlevelmax .and..not.static)then
        call refine_coarse
        do ilevel=1,levelmin
           call build_comm(ilevel)
           call make_virtual_fine_int(cpu_map(1),ilevel)
           if(hydro)then
              do ivar=1,nvar
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
              end do
              if(simple_boundary)call make_boundary_hydro(ilevel)
           endif
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              do idim=1,ndim
                 call make_virtual_fine_dp(f(1,idim),ilevel)
              end do
           end if
           if(ilevel<levelmin)call refine_fine(ilevel)
        end do
     endif

     ! Call base level
     call amr_step(levelmin,1)

     if(levelmin.lt.nlevelmax .and..not. static)then
        do ilevel=levelmin-1,1,-1
           ! Hydro book-keeping
           if(hydro)then
              call upload_fine(ilevel)
              do ivar=1,nvar
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
              end do
              if(simple_boundary)call make_boundary_hydro(ilevel)
           end if
           ! Gravity book-keeping
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              do idim=1,ndim
                 call make_virtual_fine_dp(f(1,idim),ilevel)
              end do
           end if
        end do
        
        ! Build refinement map
        do ilevel=levelmin-1,1,-1
           call flag_fine(ilevel,2)
        end do
        call flag_coarse
     endif

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
#endif

  end do

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

#endif

end subroutine adaptive_loop
