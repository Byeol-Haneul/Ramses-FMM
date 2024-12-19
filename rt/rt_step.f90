subroutine m_rt_step(pst,ilevel)
  use ramses_commons, only: pst_t

  real(dp) :: dt_hydro, t_left, dt_rt, t_save
  integer  :: i_substep, ivar

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,mdl=>pst%s%mdl)

  ! Store hydro timestep length
  dt_hydro = g%dtnew(ilevel)

  ! Time left until end of hydro time step
  t_left = dt_hydro

  ! We shift the time backwards one hydro-dt, to get evolution of stellar
  ! ages within the hydro timestep, in the case of rt subcycling:
  t_save=g%t
  g%t=g%t-t_left

  ! RT sub-cycle loop
  i_substep = 0
  do while (t_left > 0)

     ! Get RT courant time step at coarse level
     call rt_courant_coarse(dt_rt)

     ! Temporarily change timestep length to rt step
     g%dtnew(ilevel) = MIN(t_left, dt_rt/2**(ilevel-r%levelmin))

     ! Shift the time forwards one dt_rt
     g%t = g%t + g%dtnew(ilevel)

     ! Set rtunew equal to rtuold
     call rt_set_unew(ilevel)

     ! Radiative feedback from stars and sinkss
     if(rt_star) call star_RT_feedback(ilevel,g%dtnew(ilevel))
     if(rt_sink) call sink_RT_feedback(ilevel,g%dtnew(ilevel))

     ! Hyperbolic RT solver
     if(r%rt_advect) call r_rt_godunov_fine(pst,ilevel,1)

     ! Injection from radiation sources
     call add_rt_sources(ilevel,g%dtnew(ilevel))

     ! Set rtuold equal to rtunew
     call rt_set_uold(ilevel)

     ! Source terms for photo-chemistry
     if(neq_chem.or.cooling.or.T2_star>0.0.or.barotropic_eos)call cooling_fine(ilevel)

     ! Update time left until end of hydro step
     t_left = t_left - g%dtnew(ilevel)

  end do
  ! End RT subcycle loop

  ! Restore hydro timestep length
  g%dtnew(ilevel) = dt_hydro 
  g%t = t_save ! Restore original time (otherwise tiny roundoff error)

  ! Restriction operator to update coarser level split cells
  call rt_upload_fine(ilevel)

  if (myid==1 .and. rt_nsubcycle .gt. 1) write(*,901) ilevel, i_substep
901 format (' Performed level', I3, ' RT-step with ', I5, ' subcycles')

  end associate

end subroutine m_rt_step
