recursive subroutine amr_step(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  integer::i,idim,ivar,info, ilev, j
  logical::ok_defrag
  logical,save::first_step=.true., use_histograms
  real(dp)::told,tnew,dthilbert,dtrho
  real(dp), dimension(1:nvector, 1:3)::pos
  integer, dimension(1:nvector)::cell_level, cell_index, cc
  integer::k, nmax
  if(numbtot(1,ilevel)==0)return

  if(verbose)write(*,999)icount,ilevel

  !-------------------------------------------
  ! Make new refinements and update boundaries
  !-------------------------------------------
  if(levelmin.lt.nlevelmax .and..not. static)then
     if(ilevel==levelmin.or.icount>1)then
        do i=ilevel,nlevelmax
           if(i>levelmin)then

              !--------------------------
              ! Build communicators
              !--------------------------
              call build_comm(i)

               !--------------------------
              ! Update boundaries
              !--------------------------
              call make_virtual_fine_int(cpu_map(1),i)
              if(hydro)then
                 do ivar=1,nvar
                    call make_virtual_fine_dp(uold(1,ivar),i)
                 end do
                 if(simple_boundary)call make_boundary_hydro(i)
              end if
              if(poisson)then
                 call make_virtual_fine_dp(phi(1),i)
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),i)
                 end do
                 if(simple_boundary)call make_boundary_force(i)
              end if
           end if

           !--------------------------
           ! Refine grids
           !--------------------------
           call refine_fine(i)
        end do
     end if
  end if

  !--------------------------
  ! Load balance
  !--------------------------
  ok_defrag=.false.
  if(levelmin.lt.nlevelmax)then
     if(ilevel==levelmin)then
        if(nremap>0)then
           ! Skip first load balance because 
           ! it has been performed before file dump
           if(nrestart>0.and.first_step)then
              first_step=.false.
           else
              if(MOD(nstep_coarse,nremap)==0)then
                 call load_balance
                 call defrag
                 ok_defrag=.true.
              endif
           end if
        end if
     endif
  end if

  call cmp_particle_boundary_key
    
  !------------------------
  ! Output results to files
  !------------------------
  if(ilevel==levelmin)then
     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout))then
        if(.not.ok_defrag)then
           call defrag
        endif
        call dump_all
     endif
  endif
  
  !----------------------------
  ! Output frame to movie dump (without synced levels)
  !----------------------------
  if(movie) then
     if(imov.le.imovout)then 
        if(aexp>=amovout(imov).or.t>=tmovout(imov))then
           call output_frame()
        endif
     endif
  end if

  ! particles must be re-sorted before density is computed
  if(pic .and. (ilevel == levelmin .or. icount > 1))then
     use_histograms = .true.
     do ilev=ilevel, nlevelmax
        call sort_particles(ilev, use_histograms)
     end do
  end if

  !--------------------
  ! Poisson source term
  !--------------------
  if(poisson)then
     !save old potential for time-extrapolation at level boundaries
     call save_phi_old(ilevel)

     ! Rho is now updated for ALL levels >= ilevel at one call
     if(ilevel==levelmin.or.icount>1)then  
        call rho_fine(ilevel)
     end if
  endif
           
  !---------------
  ! Gravity update
  !---------------
  if(poisson)then
!     if (nstep_coarse == 10 .and. ilevel == nlevelmax)then
!        do i=1,size(rho)
!        print*,"phi", i,t, phi(i)
!           if (rho(i) > 0.d0)print*,"rho", i,nstep_coarse, rho(i)
!        end do
!     end if
        ! Remove gravity source term with half time step and old force
     if(hydro)then
        call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel))
     endif
     
     if (myid==1)print*,'here comes rho', ilevel
     nmax = 2**ilevel
     do i=1,nmax
        do j=1,nmax
           do k=1,nmax
              pos(1,1)=(i-0.5)/nmax
              pos(1,2)=(j-0.5)/nmax
              pos(1,3)=(k-0.5)/nmax
              call cmp_cpumap(pos,cc,1)
              if (cc(1)==myid)then
                 call get_cell_index(cell_index, cell_level, pos, ilevel, 1 )
                 print*,'rho', i,j,k, rho(cell_index(1)), ilevel, cell_level(1)
                 write(*,'(A,3(X,I5),X,E16.8E2,2(X,I5))'),'rho', i,j,k, rho(cell_index(1)), ilevel, cell_level(1) 
              end if
           end do
        end do
     end do
     if (myid==1)print*,'here comes phi - 1', ilevel
     nmax = 2**(ilevel-1)
     do i=1,nmax
        do j=1,nmax
           do k=1,nmax
              pos(1,1)=(i-0.5)/nmax
              pos(1,2)=(j-0.5)/nmax
              pos(1,3)=(k-0.5)/nmax
              call cmp_cpumap(pos,cc,1)
              if (cc(1)==myid)then
                 call get_cell_index(cell_index, cell_level, pos, ilevel-1, 1 )
                 write(*,'(A,3(X,I5),X,E16.8E2,2(X,I5))'),'phi', i,j,k, phi(cell_index(1)), ilevel-1, cell_level(1) 
              end if
           end do
        end do
     end do

     ! Compute gravitational potential
     if(ilevel>levelmin)then
        if(ilevel .ge. cg_levelmin) then
           call phi_fine_cg(ilevel,icount)
        else
           call multigrid_fine(ilevel,icount)
        end if
     else
        call multigrid_fine(levelmin,icount)
     end if

     !probe field
     if (myid==1)print*,'here comes phi', ilevel
     nmax = 2**ilevel
     do i=1,nmax
        do j=1,nmax
           do k=1,nmax
              pos(1,1)=(i-0.5)/nmax
              pos(1,2)=(j-0.5)/nmax
              pos(1,3)=(k-0.5)/nmax
              call cmp_cpumap(pos,cc,1)
              if (cc(1)==myid)then
                 call get_cell_index(cell_index, cell_level, pos, ilevel, 1 )
                 write(*,'(A,3(X,I5),X,E16.8E2,2(X,I5))'),'phi', i,j,k, phi(cell_index(1)), ilevel, cell_level(1)
              end if
           end do
        end do
     end do
     

     !when there is no old potential...
     if (nstep==0)call save_phi_old(ilevel)

     ! Compute gravitational acceleration
     call force_fine(ilevel,icount)

     ! Synchronize remaining particles for gravity
     if(pic)then
        call second_kick(ilevel)
     end if
     if (ilevel==levelmin)then
     do i=1,npartmax
        do j=1,npart
           if (idp(j)==i)then
              write(*,'(A,X,I3,X,I8,3(X,F10.6))'),"xp:",nstep_coarse,i,xp(j,:)
              write(*,'(A,X,I3,X,I8,3(X,F10.6))'),"vp:",nstep_coarse,i,vp(j,:)
           end if
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,info)
     end do
     endif

     if (myid==1)print*,'here comes f', ilevel
     nmax = 2**ilevel
     do i=1,nmax
        do j=1,nmax
           do k=1,nmax
              pos(1,1)=(i-0.5)/nmax
              pos(1,2)=(j-0.5)/nmax
              pos(1,3)=(k-0.5)/nmax
              call cmp_cpumap(pos,cc,1)
              if (cc(1)==myid)then
                 call get_cell_index(cell_index, cell_level, pos, ilevel, 1 )
                 write(*,'(A,3(X,I5),X,3(E18.10E2),2(X,I5))'),'f', i,j,k, f(cell_index(1),1:3), ilevel, cell_level(1)
              end if
           end do
        end do
     end do
     
     
     if(hydro)then

        ! Add gravity source term with half time step and new force
        call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))

        ! Update boundaries
        do ivar=1,nvar
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
        end do
        if(simple_boundary)call make_boundary_hydro(ilevel)
     end if
  end if

  !----------------------
  ! Compute new time step
  !----------------------
  call newdt_fine(ilevel)
  if(ilevel>levelmin)then
     dtnew(ilevel)=MIN(dtnew(ilevel-1)/real(nsubcycle(ilevel-1)),dtnew(ilevel))
  end if


  ! Set unew equal to uold
  if(hydro)call set_unew(ilevel)

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
  if(ilevel<nlevelmax)then
     if(numbtot(1,ilevel+1)>0)then
        if(nsubcycle(ilevel)==2)then
           call amr_step(ilevel+1,1)
           call amr_step(ilevel+1,2)
        else
           call amr_step(ilevel+1,1)
        endif
     else 
        ! Otherwise, update time and finer level time-step
        dtold(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        dtnew(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        call update_time(ilevel)
     end if
  else
     call update_time(ilevel)
  end if

  !---------------
  ! Move particles
  !---------------
  if(pic)then
     call kick_drift(ilevel) ! Only remaining particles
  end if

  if (ilevel==levelmin)then
     do i=1,npartmax
        do j=1,npart
           if (idp(j)==i)then
              write(*,'(A,X,I3,X,I8,3(X,F10.6))'),"xp2:",nstep_coarse,i,xp(j,:)
              write(*,'(A,X,I3,X,I8,3(X,F10.6))'),"vp2:",nstep_coarse,i,vp(j,:)
           end if
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,info)
     end do
  endif
  !-----------
  ! Hydro step
  !-----------
  if(hydro)then

     ! Hyperbolic solver
     call godunov_fine(ilevel)

     ! Reverse update boundaries
     do ivar=1,nvar
        call make_virtual_reverse_dp(unew(1,ivar),ilevel)
     end do
     if(pressure_fix)then
        call make_virtual_reverse_dp(enew(1),ilevel)
        call make_virtual_reverse_dp(divu(1),ilevel)
     endif

     ! Set uold equal to unew
     call set_uold(ilevel)

     ! Add gravity source term with half time step and old force
     ! in order to complete the time step 
     if(poisson)call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))

     ! Restriction operator
     call upload_fine(ilevel)

  endif

  !---------------------------------------
  ! Update physical and virtual boundaries
  !---------------------------------------
  if(hydro)then
     do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     end do
     if(simple_boundary)call make_boundary_hydro(ilevel)
  endif

  !-----------------------
  ! Compute refinement map
  !-----------------------
  if(.not.static) call flag_fine(ilevel,icount)

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>levelmin)then
     if(nsubcycle(ilevel-1)==1)dtnew(ilevel-1)=dtnew(ilevel)
     if(icount==2)dtnew(ilevel-1)=dtold(ilevel)+dtnew(ilevel)
  end if

999 format(' Entering amr_step',i1,' for level',i2)

end subroutine amr_step




