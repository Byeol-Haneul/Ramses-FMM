subroutine adaptive_loop(test_ok)
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::test_ok
  integer::ilevel,idim,ivar,info
  real(kind=8)::tt1,tt2
  real(kind=4)::real_mem,real_mem_tot

#ifndef WITHOUTMPI
  tt1=MPI_WTIME(info)
#endif

  call init_amr                      ! Initialize AMR variables
  call init_time                     ! Initialize time variables
  if(hydro)call init_hydro           ! Initialize hydro variables
  if(poisson)call init_poisson       ! Initialize poisson variables
  if(nrestart==0)call init_refine    ! Build initial AMR grid
  if(pic)call init_part              ! Initialize particle variables
  if(pic)call init_part_andreas      ! Initialize Andreas particle variables
  if(pic)call init_tree              ! Initialize particle tree
  if(nrestart==0)call init_refine_2  ! Build initial AMR grid again

#ifndef WITHOUTMPI
  tt2=MPI_WTIME(info)
!  if(myid==1)write(*,*)'Time elapsed since startup:',tt2-tt1
#endif

!  if(myid==1)then
!     write(*,*)'Initial mesh structure'
!     do ilevel=1,nlevelmax
!        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
!     end do
!  end if

  nstep_coarse_old=nstep_coarse

!  if(myid==1)write(*,*)'Starting time integration' 


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
  call amr_step(levelmin,1,test_ok)
  

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine adaptive_loop
