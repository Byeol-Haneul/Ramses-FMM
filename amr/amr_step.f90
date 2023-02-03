module amr_step
contains
!#####################################################
!#####################################################
!#####################################################
!#####################################################
recursive subroutine m_amr_step(pst,ilevel,icount,done)
  use ramses_commons, only: pst_t
  use pm_parameters
  use flag_utils, only: m_flag_fine
  use update_time_module, only: m_update_time
  use refine_utils, only: m_refine_fine
  use upload_module, only: m_upload_fine
#ifdef GRAV
  use rho_fine_module, only: m_rho_fine
  use phi_fine_cg_module, only: m_phi_fine_cg
  use multigrid_fine_commons, only: multigrid
  use force_fine_module, only: m_force_fine
#endif
  use move_fine_module, only: m_kick_drift_part
  use output_amr_module, only: m_dump_all
  use synchro_hydro_fine_module, only: m_synchro_hydro_fine, r_gravity_hydro_fine
  use nbors_utils, only: r_save_phi_old
  use godunov_fine_module, only: r_godunov_fine,r_set_unew,r_set_uold
  use cooling_fine_module, only: r_cooling_fine
  use newdt_fine_module, only: m_newdt_fine,r_broadcast_dt,in_broadcast_dt_t
  use movie_module, only: m_output_frame

  implicit none

  type(pst_t)::pst
  integer::ilevel,icount
  logical::done
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  type(in_broadcast_dt_t)::in_broadcast_dt

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,mdl=>pst%s%mdl)
  
  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,'(" Entering amr_step",i1," for level",i2)')icount,ilevel

  !------------------------------
  ! Make new refinements and load
  ! balance grids and particles.
  !------------------------------
  if(ilevel==r%levelmin.or.icount>1)then
                                    call m_timer(pst,'refine','start')
     call m_refine_fine(pst,ilevel)
  endif
  
  !------------------------
  ! Output results to files
  !------------------------
  if(ilevel==r%levelmin)then
     if(r%foutput>0)then
        if(mod(g%nstep_coarse,r%foutput)==0.or.g%aexp>=r%aout(g%iout).or.g%t>=r%tout(g%iout))then
                                    call m_timer(pst,'output','start')
           call m_dump_all(pst)
        endif
     endif
  endif
  
  !----------------------------
  ! Output frame to movie dump
  !----------------------------
  if(r%movie) then
     if(r%imov.le.r%imovout)then 
        if(g%aexp>=r%aendmov*dble(r%imov)/dble(r%imovout).or.g%t>=r%tendmov*dble(r%imov)/dble(r%imovout))then
           call m_output_frame(pst)
        endif
     endif
  end if

  !--------------------
  ! Poisson source term
  !--------------------
#ifdef GRAV
  if(r%poisson)then
     if(ilevel==r%levelmin.or.icount>1)then
                                    call m_timer(pst,'rho','start')
        call m_rho_fine(pst,ilevel)
     endif
  endif
#endif
  !---------------
  ! Gravity solver
  !---------------
#ifdef GRAV
  if(r%poisson)then

                                    call m_timer(pst,'poisson','start')
     ! Remove gravity source term with half time step and old force
     if(r%hydro)then
        call m_synchro_hydro_fine(pst,ilevel,-0.5d0*g%dtnew(ilevel))
     endif

     ! Save old potential for time-extrapolation at level boundaries
     call r_save_phi_old(pst,ilevel,1)

     ! Compute new gravitational potential
     if(ilevel > r%levelmin)then
        if(ilevel >= r%cg_levelmin) then
           call m_phi_fine_cg(pst,ilevel,icount)
        else
           call multigrid(pst,ilevel,icount)
        end if
     else
        call multigrid(pst,r%levelmin,icount)
     end if

     ! Initial old potential
     if (g%nstep==0)call r_save_phi_old(pst,ilevel,1)

     ! Compute gravitational acceleration
     call m_force_fine(pst,ilevel,icount)

     ! Perform second kick for particles
     if(r%pic)then
                                    call m_timer(pst,'particle - kickdrift','start')
        call m_kick_drift_part(pst,ilevel,action_kick_only)
     endif
     ! Add gravity source term with half time step and new force
     if(r%hydro)then
                                    call m_timer(pst,'poisson','start')
        call m_synchro_hydro_fine(pst,ilevel,+0.5d0*g%dtnew(ilevel))
     end if

  end if
#endif

  !----------------------
  ! Compute new time step
  !----------------------
                                    call m_timer(pst,'courant','start')
  call m_newdt_fine(pst,ilevel)
  
  !-----------------------
  ! Set unew equal to uold
  !-----------------------
  if(r%hydro)then
                                    call m_timer(pst,'hydro - set unew','start')
     call r_set_unew(pst,ilevel,1)
  endif

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
                                    call m_timer(pst,'recursive call','start')
  if(ilevel<r%nlevelmax)then
     if(m%noct_tot(ilevel+1)>0)then
        if(r%nsubcycle(ilevel)==2)then
           call m_amr_step(pst,ilevel+1,1,done)
           if (done)return
           call m_amr_step(pst,ilevel+1,2,done)
        else
           call m_amr_step(pst,ilevel+1,1,done)
        endif
     else 
        ! Otherwise, modify finer level time-step
        g%dtold(ilevel+1)=g%dtnew(ilevel)/dble(r%nsubcycle(ilevel))
        g%dtnew(ilevel+1)=g%dtnew(ilevel)/dble(r%nsubcycle(ilevel))

        ! Broadcast modified time step to all CPUs
        in_broadcast_dt%ilevel=ilevel+1
        in_broadcast_dt%dtnew=g%dtnew(ilevel+1)
        in_broadcast_dt%dtold=g%dtold(ilevel+1)
        call r_broadcast_dt(pst,in_broadcast_dt,storage_size(in_broadcast_dt)/32)

        ! Update time variable
        call m_update_time(pst,ilevel,done)
     end if
  else
     call m_update_time(pst,ilevel,done)
  end if
  if (done)return

  !-----------
  ! Hydro step
  !-----------
  if(r%hydro)then

     ! Hyperbolic solver
                                    call m_timer(pst,'hydro - godunov','start')
     if(.not.r%static)call r_godunov_fine(pst,ilevel,1)

     ! Add gravity source terms to unew with half time step
     if(r%poisson)then
                                    call m_timer(pst,'hydro - gravity','start')
        call r_gravity_hydro_fine(pst,ilevel,1)
     endif

     ! Set uold equal to unew
                                    call m_timer(pst,'hydro - set uold','start')
     call r_set_uold(pst,ilevel,1)

     ! Add gravity source terms to uold with half time step
     ! to complete the time step with old force (will be removed later)
     if(r%poisson)then
                                    call m_timer(pst,'poisson - synchro','start')
        call m_synchro_hydro_fine(pst,ilevel,+0.5d0*g%dtnew(ilevel))
     endif

     ! Restriction operator
                                    call m_timer(pst,'hydro - upload','start')
     call m_upload_fine(pst,ilevel)
  endif

  !----------------------------
  ! Compute cooling/heating
  !----------------------------
  if(r%cooling)then
                                    call m_timer(pst,'cooling','start')
     call r_cooling_fine(pst,ilevel,1)
  endif

  !-------------------------------------------
  ! Perform first kick and drift for particles
  !-------------------------------------------
  if(r%pic)then
                                    call m_timer(pst,'particle - kickdrift','start')
     call m_kick_drift_part(pst,ilevel,action_kick_drift)
  endif

  !-----------------------
  ! Compute refinement map
  !-----------------------
                                    call m_timer(pst,'flag','start')
  if(.not.r%static)call m_flag_fine(pst,ilevel,icount)

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>r%levelmin)then
     ! Impose adaptive time step constraints
     if(r%nsubcycle(ilevel-1)==1)g%dtnew(ilevel-1)=g%dtnew(ilevel)
     if(icount==2)g%dtnew(ilevel-1)=g%dtold(ilevel)+g%dtnew(ilevel)

     ! Broadcast updated time step to all CPUs
                                    call m_timer(pst,'recursive call','start')
     in_broadcast_dt%ilevel=ilevel-1
     in_broadcast_dt%dtnew=g%dtnew(ilevel-1)
     in_broadcast_dt%dtold=g%dtold(ilevel-1)
     call r_broadcast_dt(pst,in_broadcast_dt,storage_size(in_broadcast_dt)/32)
  end if

  end associate

end subroutine m_amr_step
end module amr_step
