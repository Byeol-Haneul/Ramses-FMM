recursive subroutine amr_step(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount,ilev
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!

  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,999)icount,ilevel

  !---------------------
  ! Make new refinements
  !---------------------
  if(ilevel==levelmin.or.icount>1)then
                               call timer('refine','start')
     call refine_fine(ilevel)
                               call timer('load balance','start')
     call load_balance(ilevel)
  endif

  !------------------------
  ! Output results to files
  !------------------------
  if(ilevel==levelmin)then
     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout))then
        call dump_all
     endif
  endif
  
  !----------------------------
  ! Output frame to movie dump
  !----------------------------
  if(movie) then
     if(imov.le.imovout)then 
        if(aexp>=amovout(imov).or.t>=tmovout(imov))then
           call output_frame()
        endif
     endif
  end if

  !--------------------
  ! Poisson source term
  !--------------------
  if(poisson)then
     if(ilevel==levelmin.or.icount>1)then
                               call timer('poisson - rho','start')
        call rho_fine(ilevel)
     endif
  endif

  !---------------
  ! Gravity solver
  !---------------
#ifdef GRAV
  if(poisson)then

     ! Remove gravity source term with half time step and old force
     if(hydro)then
                               call timer('poisson - synchro','start')
        call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel))
     endif

     ! Compute gravitational potential
                               call timer('poisson - solver','start')
     ! Save old potential for time-extrapolation at level boundaries
     call save_phi_old(ilevel)
     if(ilevel > levelmin)then
        if(ilevel >= cg_levelmin) then
           call phi_fine_cg(ilevel,icount)
        else
           call multigrid(ilevel,icount)
        end if
     else
        call multigrid(levelmin,icount)
     end if
     ! Initial old potential
     if (nstep==0)call save_phi_old(ilevel)

     ! Compute gravitational acceleration
     call force_fine(ilevel,icount)

     ! Perform second kick for particles
                               call timer('particles','start')
     call kick_drift_part(ilevel,action_kick_only)

     ! Add gravity source term with half time step and new force
     if(hydro)then
                               call timer('poisson - synchro','start')
        call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))
     end if

  end if
#endif

  !----------------------
  ! Compute new time step
  !----------------------
                               call timer('courant','start')
  call newdt_fine(ilevel)
  if(ilevel>levelmin)then
     dtnew(ilevel)=MIN(dtnew(ilevel-1)/real(nsubcycle(ilevel-1)),dtnew(ilevel))
  end if
  
  !-----------------------
  ! Set unew equal to uold
  !-----------------------
                               call timer('hydro - set unew','start')
  if(hydro)call set_unew(ilevel)

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
                               call timer('recursive call','start')
  if(ilevel<nlevelmax)then
     if(noct_tot(ilevel+1)>0)then
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

  !-----------
  ! Hydro step
  !-----------
  if(hydro)then
     ! Hyperbolic solver
                               call timer('hydro - godunov','start')
     call godunov_fine(ilevel)
     ! Add gravity source terms to unew with half time step
                               call timer('poisson - synchro','start')
     if(poisson)call add_gravity_source_terms(ilevel)

     ! Set uold equal to unew
                               call timer('hydro - set uold','start')
     call set_uold(ilevel)
     ! Add gravity source terms to uold with half time step
     ! to complete the time step (will be removed later)
                               call timer('poisson - synchro','start')
     if(poisson)call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))
     ! Restriction operator
                               call timer('hydro - upload','start')
     call upload_fine(ilevel)
  endif

  !----------------------------
  ! Compute cooling/heating
  !----------------------------
                               call timer('cooling','start')
  if(cooling)call cooling_fine(ilevel)

  !-------------------------------------------
  ! Perform first kick and drift for particles
  !-------------------------------------------
                               call timer('particles','start')
  call kick_drift_part(ilevel,action_kick_drift)

  !-----------------------
  ! Compute refinement map
  !-----------------------
                               call timer('flag','start')
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




