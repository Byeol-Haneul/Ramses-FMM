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
  integer::i,idim,ivar,info, idomain, ilev
  logical::ok_defrag
  logical,save::first_step=.true.
  real(dp)::told,tnew,dthilbert,dtrho
  logical::use_histograms

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


  do idomain=1,ndomain - 1
     bound_key_level(idomain,nlevelmax) = nint(bound_key(idomain) / 8.)
  end do
  bound_key_level(0,nlevelmax) = floor(bound_key(0) / 8.)
  bound_key_level(ndomain,nlevelmax) = ceiling(bound_key(ndomain) / 8.)
  
  do ilev=nlevelmax-1, levelmin, - 1
     do idomain=1,ndomain - 1
        bound_key_level(idomain,ilev) = nint(bound_key_level(idomain, ilev +1) / 8.)
     end do
     bound_key_level(0,ilev) = floor(bound_key_level(0, ilev +1) / 8.)
     bound_key_level(ndomain,ilev) = ceiling(bound_key_level(ndomain, ilev +1) / 8.)
  end do
  
  call MPI_BARRIER(MPI_COMM_WORLD,info)
  if (myid==1)print*, 'starting tests now', npart_andreas


  
  !-----------------
  ! Particle leakage
  !-----------------
  if(pic)call make_tree_fine(ilevel)
  
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

  !-----------------------------------------------------------
  ! Put here all stuffs that are done only at coarse time step
  !-----------------------------------------------------------

  !--------------------
  ! Poisson source term
  !--------------------
  if(poisson)then
     !save old potential for time-extrapolation at level boundaries
     call save_phi_old(ilevel)

! #ifndef WITHOUTMPI
!      told=MPI_WTIME(info)
! #endif

     call rho_fine(ilevel,icount,.false.)     

! #ifndef WITHOUTMPI
!      tnew=MPI_WTIME(info)
!      dtrho=tnew-told
!      told=tnew
! #endif

!      if(pic)call hilbert_allparts(ilevel)

! #ifndef WITHOUTMPI
!      tnew=MPI_WTIME(info)
!      dthilbert=tnew-told
!      told=tnew
! #endif

!      if (dthilbert>0. .and. dtrho > 0.)then
!         print*,'dt(hilbert+sort)/dt(rho): ',dthilbert/dtrho,myid,ilevel
!      end if
  endif
!  stop
  !-------------------------------------------
  ! Sort particles between ilevel and ilevel+1
  !-------------------------------------------
  if(pic)then
     ! Remove particles to finer levels
     call kill_tree_fine(ilevel)
     ! Update boundary conditions for remaining particles
     call virtual_tree_fine(ilevel)
  end if

  !---------------
  ! Gravity update
  !---------------
  if(poisson)then
 
     ! Remove gravity source term with half time step and old force
     if(hydro)then
        call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel))
     endif
     
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
     !when there is no old potential...
     if (nstep==0)call save_phi_old(ilevel)

     ! Compute gravitational acceleration
     call force_fine(ilevel,icount)


     use_histograms = .false.
     call sort_particles(ilevel, use_histograms)
     
     ! Synchronize remaining particles for gravity
     if(pic)then
        call synchro_fine(ilevel)
        call second_kick(ilevel)
     end if

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
     call move_fine(ilevel) ! Only remaining particles
     call kick_drift(ilevel) ! Only remaining particles
  end if

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

  !----------------------------
  ! Merge finer level particles
  !----------------------------
  if(pic)call merge_tree_fine(ilevel)

  if (ilevel==levelmin)then
     call write_ascii_parts
     print*,'written ascii', t
     if(t>0.2)stop
  endif
  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>levelmin)then
     if(nsubcycle(ilevel-1)==1)dtnew(ilevel-1)=dtnew(ilevel)
     if(icount==2)dtnew(ilevel-1)=dtold(ilevel)+dtnew(ilevel)
  end if

999 format(' Entering amr_step',i1,' for level',i2)

end subroutine amr_step




