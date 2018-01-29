!#####################################################
!#####################################################
!#####################################################
!#####################################################
recursive subroutine m_amr_step(s,ilevel,icount)
  use ramses_commons, only: ramses_t
  use pm_parameters
  implicit none
  type(ramses_t)::s
  integer::ilevel,icount
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  integer,dimension(1:5)::input_array
  
  if(s%m%noct_tot(ilevel)==0)return
  if(s%r%verbose)write(*,'(" Entering amr_step",i1," for level",i2)')icount,ilevel

  !------------------------------
  ! Make new refinements and load
  ! balance grids and particles.
  !------------------------------
  if(ilevel==s%r%levelmin.or.icount>1)then
     call m_refine_fine(s,ilevel)
  endif
  
  !------------------------
  ! Output results to files
  !------------------------
  if(ilevel==s%r%levelmin)then
     if(mod(s%g%nstep_coarse,s%r%foutput)==0.or.s%g%aexp>=s%r%aout(s%g%iout).or.s%g%t>=s%r%tout(s%g%iout))then
        call m_dump_all(s)
     endif
  endif
  
  !----------------------------
  ! Output frame to movie dump
  !----------------------------
  if(s%r%movie) then
     if(s%r%imov.le.s%r%imovout)then 
        if(s%g%aexp>=s%r%amovout(s%r%imov).or.s%g%t>=s%r%tmovout(s%r%imov))then
           call m_output_frame(s)
        endif
     endif
  end if

  !--------------------
  ! Poisson source term
  !--------------------
  if(s%r%poisson)then
     if(ilevel==s%r%levelmin.or.icount>1)then
        call m_rho_fine(s,ilevel)
     endif
  endif

  !---------------
  ! Gravity solver
  !---------------
#ifdef GRAV
  if(s%r%poisson)then

     ! Remove gravity source term with half time step and old force
     if(s%r%hydro)then
        call m_synchro_hydro_fine(s,ilevel,-0.5d0*s%g%dtnew(ilevel))
     endif

     ! Save old potential for time-extrapolation at level boundaries
     call r_save_phi_old(s,s%mdl%ncpu,1,0,ilevel)

     ! Compute new gravitational potential
     if(ilevel > s%r%levelmin)then
        if(ilevel >= s%r%cg_levelmin) then
           call m_phi_fine_cg(s,ilevel,icount)
        else
           call multigrid(s,ilevel,icount)
        end if
     else
        call multigrid(s,s%r%levelmin,icount)
     end if

     ! Initial old potential
     if (s%g%nstep==0)call r_save_phi_old(s,s%mdl%ncpu,1,0,ilevel)

     ! Compute gravitational acceleration
     call m_force_fine(s,ilevel,icount)

     ! Perform second kick for particles
     if(s%r%pic)call m_kick_drift_part(s,ilevel,action_kick_only)

     ! Add gravity source term with half time step and new force
     if(s%r%hydro)then
        call m_synchro_hydro_fine(s,ilevel,+0.5d0*s%g%dtnew(ilevel))
     end if

  end if
#endif

  !----------------------
  ! Compute new time step
  !----------------------
  call m_newdt_fine(s,ilevel)
  
  !-----------------------
  ! Set unew equal to uold
  !-----------------------
  if(s%r%hydro)call r_set_unew(s,s%mdl%ncpu,1,0,ilevel)

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
  if(ilevel<s%r%nlevelmax)then
     if(s%m%noct_tot(ilevel+1)>0)then
        if(s%r%nsubcycle(ilevel)==2)then
           call m_amr_step(s,ilevel+1,1)
           call m_amr_step(s,ilevel+1,2)
        else
           call m_amr_step(s,ilevel+1,1)
        endif
     else 
        ! Otherwise, modify finer level time-step
        s%g%dtold(ilevel+1)=s%g%dtnew(ilevel)/dble(s%r%nsubcycle(ilevel))
        s%g%dtnew(ilevel+1)=s%g%dtnew(ilevel)/dble(s%r%nsubcycle(ilevel))

        ! Broadcast modified time step to all CPUs
        input_array(1)=ilevel+1
        input_array(2:3)=transfer(s%g%dtnew(ilevel+1),input_array)
        input_array(4:5)=transfer(s%g%dtold(ilevel+1),input_array)
        call r_broadcast_dt(s,s%mdl%ncpu,5,0,input_array)

        ! Update time variable
        call m_update_time(s,ilevel)
     end if
  else
     call m_update_time(s,ilevel)
  end if

  !-----------
  ! Hydro step
  !-----------
  if(s%r%hydro)then

     ! Hyperbolic solver
     if(.not.s%r%static)call r_godunov_fine(s,s%mdl%ncpu,1,0,ilevel)

     ! Add gravity source terms to unew with half time step
     if(s%r%poisson)call r_gravity_hydro_fine(s,s%mdl%ncpu,1,0,ilevel)

     ! Set uold equal to unew
     call r_set_uold(s,s%mdl%ncpu,1,0,ilevel)

     ! Add gravity source terms to uold with half time step
     ! to complete the time step with old force (will be removed later)
     if(s%r%poisson)call m_synchro_hydro_fine(s,ilevel,+0.5d0*s%g%dtnew(ilevel))

     ! Restriction operator
     call m_upload_fine(s,ilevel)
  endif

  !----------------------------
  ! Compute cooling/heating
  !----------------------------
  if(s%r%cooling)call r_cooling_fine(s,s%mdl%ncpu,1,0,ilevel)

  !-------------------------------------------
  ! Perform first kick and drift for particles
  !-------------------------------------------
  if(s%r%pic)call m_kick_drift_part(s,ilevel,action_kick_drift)

  !-----------------------
  ! Compute refinement map
  !-----------------------
  if(.not.s%r%static)call m_flag_fine(s,ilevel,icount)

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>s%r%levelmin)then
     ! Impose adaptive time step constraints
     if(s%r%nsubcycle(ilevel-1)==1)s%g%dtnew(ilevel-1)=s%g%dtnew(ilevel)
     if(icount==2)s%g%dtnew(ilevel-1)=s%g%dtold(ilevel)+s%g%dtnew(ilevel)

     ! Broadcast updated time step to all CPUs
     input_array(1)=ilevel-1
     input_array(2:3)=transfer(s%g%dtnew(ilevel-1),input_array)
     input_array(4:5)=transfer(s%g%dtold(ilevel-1),input_array)
     call r_broadcast_dt(s,s%mdl%ncpu,5,0,input_array)
  end if

end subroutine m_amr_step
