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
  ! Output frame to movie dump (without synced levels)
  !----------------------------
!!$  if(movie) then
!!$     if(imov.le.imovout)then 
!!$        if(aexp>=amovout(imov).or.t>=tmovout(imov))then
!!$           call output_frame()
!!$        endif
!!$     endif
!!$  end if

  !--------------------
  ! Poisson source term
  !--------------------
  if(poisson)then
                               call timer('poisson - rho','start')
     ! Compute total mass density at level ilevel
     call rho_fine(ilevel,icount)
  endif

  !---------------
  ! Gravity solver
  !---------------
  if(poisson)then
                               call timer('poisson - synchro','start')
     ! Remove gravity source term with half time step and old force
     if(hydro)then
        call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel))
     endif

                               call timer('poisson - solver','start')
     ! Save old potential for time-extrapolation at level boundaries
     call save_phi_old(ilevel)

     ! Compute gravitational potential
     if(ilevel > levelmin)then
        if(ilevel >= cg_levelmin) then
           call phi_fine_cg(ilevel,icount)
        else
!           call multigrid(ilevel,icount)
        end if
     else
        call phi_fine_cg(ilevel,icount)
!        call multigrid(levelmin,icount)
     end if

     ! Initial old potential
     if (nstep==0)call save_phi_old(ilevel)

     ! Compute gravitational acceleration
     call force_fine(ilevel,icount)

     if(hydro)then
                               call timer('poisson - synchro','start')
        ! Add gravity source term with half time step and new force
        call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))
     end if

  end if
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




