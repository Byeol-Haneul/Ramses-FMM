!#####################################################
!#####################################################
!#####################################################
!#####################################################
recursive subroutine m_amr_step(pst,ilevel,icount)
  use ramses_commons, only: pst_t
  use pm_parameters
  implicit none
  type(pst_t)::pst
  integer::ilevel,icount
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  integer,dimension(1:5)::input_array

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,mdl=>pst%s%mdl)
  
  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,'(" Entering amr_step",i1," for level",i2)')icount,ilevel

  !------------------------------
  ! Make new refinements and load
  ! balance grids and particles.
  !------------------------------
  if(ilevel==r%levelmin.or.icount>1)then
     call m_refine_fine(pst,ilevel)
  endif
  
  !------------------------
  ! Output results to files
  !------------------------
  if(ilevel==r%levelmin)then
     if(mod(g%nstep_coarse,r%foutput)==0.or.g%aexp>=r%aout(g%iout).or.g%t>=r%tout(g%iout))then
        call m_dump_all(pst)
     endif
  endif
  
  !----------------------------
  ! Output frame to movie dump
  !----------------------------
  if(r%movie) then
     if(r%imov.le.r%imovout)then 
        if(g%aexp>=r%amovout(r%imov).or.g%t>=r%tmovout(r%imov))then
           call m_output_frame(pst)
        endif
     endif
  end if

  !--------------------
  ! Poisson source term
  !--------------------
  if(r%poisson)then
     if(ilevel==r%levelmin.or.icount>1)then
        call m_rho_fine(pst,ilevel)
     endif
  endif

  !---------------
  ! Gravity solver
  !---------------
#ifdef GRAV
  if(r%poisson)then

     ! Remove gravity source term with half time step and old force
     if(r%hydro)then
        call m_synchro_hydro_fine(pst,ilevel,-0.5d0*g%dtnew(ilevel))
     endif

     ! Save old potential for time-extrapolation at level boundaries
     call r_save_phi_old(pst,mdl%ncpu,1,0,ilevel)

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
     if (g%nstep==0)call r_save_phi_old(pst,mdl%ncpu,1,0,ilevel)

     ! Compute gravitational acceleration
     call m_force_fine(pst,ilevel,icount)

     ! Perform second kick for particles
     if(r%pic)call m_kick_drift_part(pst,ilevel,action_kick_only)

     ! Add gravity source term with half time step and new force
     if(r%hydro)then
        call m_synchro_hydro_fine(pst,ilevel,+0.5d0*g%dtnew(ilevel))
     end if

  end if
#endif

  !----------------------
  ! Compute new time step
  !----------------------
  call m_newdt_fine(pst,ilevel)
  
  !-----------------------
  ! Set unew equal to uold
  !-----------------------
  if(r%hydro)call r_set_unew(pst,mdl%ncpu,1,0,ilevel)

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
  if(ilevel<r%nlevelmax)then
     if(m%noct_tot(ilevel+1)>0)then
        if(r%nsubcycle(ilevel)==2)then
           call m_amr_step(pst,ilevel+1,1)
           call m_amr_step(pst,ilevel+1,2)
        else
           call m_amr_step(pst,ilevel+1,1)
        endif
     else 
        ! Otherwise, modify finer level time-step
        g%dtold(ilevel+1)=g%dtnew(ilevel)/dble(r%nsubcycle(ilevel))
        g%dtnew(ilevel+1)=g%dtnew(ilevel)/dble(r%nsubcycle(ilevel))

        ! Broadcast modified time step to all CPUs
        input_array(1)=ilevel+1
        input_array(2:3)=transfer(g%dtnew(ilevel+1),input_array)
        input_array(4:5)=transfer(g%dtold(ilevel+1),input_array)
        call r_broadcast_dt(pst,mdl%ncpu,5,0,input_array)

        ! Update time variable
        call m_update_time(pst,ilevel)
     end if
  else
     call m_update_time(pst,ilevel)
  end if

  !-----------
  ! Hydro step
  !-----------
  if(r%hydro)then

     ! Hyperbolic solver
     if(.not.r%static)call r_godunov_fine(pst,mdl%ncpu,1,0,ilevel)

     ! Add gravity source terms to unew with half time step
     if(r%poisson)call r_gravity_hydro_fine(pst,mdl%ncpu,1,0,ilevel)

     ! Set uold equal to unew
     call r_set_uold(pst,mdl%ncpu,1,0,ilevel)

     ! Add gravity source terms to uold with half time step
     ! to complete the time step with old force (will be removed later)
     if(r%poisson)call m_synchro_hydro_fine(pst,ilevel,+0.5d0*g%dtnew(ilevel))

     ! Restriction operator
     call m_upload_fine(pst,ilevel)
  endif

  !----------------------------
  ! Compute cooling/heating
  !----------------------------
  if(r%cooling)call r_cooling_fine(pst,mdl%ncpu,1,0,ilevel)

  !-------------------------------------------
  ! Perform first kick and drift for particles
  !-------------------------------------------
  if(r%pic)call m_kick_drift_part(pst,ilevel,action_kick_drift)

  !-----------------------
  ! Compute refinement map
  !-----------------------
  if(.not.r%static)call m_flag_fine(pst,ilevel,icount)

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>r%levelmin)then
     ! Impose adaptive time step constraints
     if(r%nsubcycle(ilevel-1)==1)g%dtnew(ilevel-1)=g%dtnew(ilevel)
     if(icount==2)g%dtnew(ilevel-1)=g%dtold(ilevel)+g%dtnew(ilevel)

     ! Broadcast updated time step to all CPUs
     input_array(1)=ilevel-1
     input_array(2:3)=transfer(g%dtnew(ilevel-1),input_array)
     input_array(4:5)=transfer(g%dtold(ilevel-1),input_array)
     call r_broadcast_dt(pst,mdl%ncpu,5,0,input_array)
  end if

  end associate

end subroutine m_amr_step
