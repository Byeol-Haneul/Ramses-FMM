!#####################################################
!#####################################################
!#####################################################
!#####################################################
recursive subroutine amr_step_2(r,g,m,p,ilevel,icount)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use pm_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  integer::ilevel,icount
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  logical,save::first_step=.true.

  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,999)icount,ilevel

  !---------------------
  ! Make new refinements
  !---------------------
  if(ilevel==r%levelmin.or.icount>1)then
                               call timer('refine','start')
     call refine_fine_2(r,g,m,ilevel)
                               call timer('load balance','start')
     call load_balance_2(r,g,m,ilevel)
  endif

  !------------------------------
  ! Balance particles across cpus
  ! Careful rho must have been called once !
  !------------------------------
  if(first_step)then
     first_step=.false.
  else
     if(ilevel==r%levelmin)then
                               call timer('particles','start')
        if(r%pic)call balance_part_2(r,g,m,p,ilevel)
     endif
  endif

  !------------------------
  ! Output results to files
  !------------------------
                               call timer('output','start')
  if(ilevel==r%levelmin)then
     if(mod(g%nstep_coarse,r%foutput)==0.or.g%aexp>=r%aout(g%iout).or.g%t>=r%tout(g%iout))then
        call dump_all_2(r,g,m,p)
     endif
  endif
  
  !----------------------------
  ! Output frame to movie dump
  !----------------------------
  if(r%movie) then
     if(r%imov.le.r%imovout)then 
        if(g%aexp>=r%amovout(r%imov).or.g%t>=r%tmovout(r%imov))then
           call output_frame_2(r,g,m,p)
        endif
     endif
  end if

  !--------------------
  ! Poisson source term
  !--------------------
  if(r%poisson)then
     if(ilevel==r%levelmin.or.icount>1)then
                               call timer('rho','start')
        call rho_fine_2(r,g,m,p,ilevel)
     endif
  endif

  !---------------
  ! Gravity solver
  !---------------
#ifdef GRAV
  if(r%poisson)then
                               call timer('poisson','start')
     ! Remove gravity source term with half time step and old force
     if(r%hydro)then
        call synchro_hydro_fine_2(r,m,ilevel,-0.5*g%dtnew(ilevel))
     endif

     ! Save old potential for time-extrapolation at level boundaries
     call save_phi_old_2(m,ilevel)

     ! Compute new gravitational potential
     if(ilevel > r%levelmin)then
        if(ilevel >= r%cg_levelmin) then
           call phi_fine_cg_2(r,g,m,ilevel,icount)
        else
           call multigrid_2(r,g,m,ilevel,icount)
        end if
     else
        call multigrid_2(r,g,m,r%levelmin,icount)
     end if

     ! Initial old potential
     if (g%nstep==0)call save_phi_old_2(m,ilevel)

     ! Compute gravitational acceleration
     call force_fine_2(r,g,m,ilevel,icount)

     ! Perform second kick for particles
                               call timer('particles','start')
     if(r%pic)call kick_drift_part_2(r,g,m,p,ilevel,action_kick_only)

     ! Add gravity source term with half time step and new force
     if(r%hydro)then
                               call timer('poisson','start')
        call synchro_hydro_fine_2(r,m,ilevel,+0.5*g%dtnew(ilevel))
     end if

  end if
#endif

  !----------------------
  ! Compute new time step
  !----------------------
                               call timer('courant','start')
  call newdt_fine_2(r,g,m,p,ilevel)
  if(ilevel>r%levelmin)then
     g%dtnew(ilevel)=MIN(g%dtnew(ilevel-1)/real(r%nsubcycle(ilevel-1)),g%dtnew(ilevel))
  end if
  
  !-----------------------
  ! Set unew equal to uold
  !-----------------------
                               call timer('hydro - set unew','start')
  if(r%hydro)call set_unew_2(r,g,m,ilevel)

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
                               call timer('recursive call','start')
  if(ilevel<r%nlevelmax)then
     if(m%noct_tot(ilevel+1)>0)then
        if(r%nsubcycle(ilevel)==2)then
           call amr_step_2(r,g,m,p,ilevel+1,1)
           call amr_step_2(r,g,m,p,ilevel+1,2)
        else
           call amr_step_2(r,g,m,p,ilevel+1,1)
        endif
     else 
        ! Otherwise, update time and finer level time-step
        g%dtold(ilevel+1)=g%dtnew(ilevel)/dble(r%nsubcycle(ilevel))
        g%dtnew(ilevel+1)=g%dtnew(ilevel)/dble(r%nsubcycle(ilevel))
        call update_time_2(r,g,m,p,ilevel)
     end if
  else
     call update_time_2(r,g,m,p,ilevel)
  end if

  !-----------
  ! Hydro step
  !-----------
  if(r%hydro)then
     ! Hyperbolic solver
                               call timer('hydro - godunov','start')
     call godunov_fine_2(r,g,m,ilevel)
     ! Add gravity source terms to unew with half time step
                               call timer('poisson - synchro','start')
     if(r%poisson)call add_gravity_source_terms_2(r,g,m,ilevel)

     ! Set uold equal to unew
                               call timer('hydro - set uold','start')
     call set_uold_2(r,g,m,ilevel)
     ! Add gravity source terms to uold with half time step
     ! to complete the time step (will be removed later)
                               call timer('poisson - synchro','start')
     if(r%poisson)call synchro_hydro_fine_2(r,m,ilevel,+0.5*g%dtnew(ilevel))
     ! Restriction operator
                               call timer('hydro - upload','start')
     call upload_fine_2(r,g,m,ilevel)
  endif

  !----------------------------
  ! Compute cooling/heating
  !----------------------------
                               call timer('cooling','start')
  if(r%cooling)call cooling_fine_2(r,g,m,ilevel)

  !-------------------------------------------
  ! Perform first kick and drift for particles
  !-------------------------------------------
                               call timer('particles','start')
  if(r%pic)call kick_drift_part_2(r,g,m,p,ilevel,action_kick_drift)

  !-----------------------
  ! Compute refinement map
  !-----------------------
                               call timer('flag','start')
  if(.not.r%static) call flag_fine_2(r,g,m,ilevel,icount)

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
                               call timer('recursive call','start')
  if(ilevel>r%levelmin)then
     if(r%nsubcycle(ilevel-1)==1)g%dtnew(ilevel-1)=g%dtnew(ilevel)
     if(icount==2)g%dtnew(ilevel-1)=g%dtold(ilevel)+g%dtnew(ilevel)
  end if

999 format(' Entering amr_step',i1,' for level',i2)

end subroutine amr_step_2




