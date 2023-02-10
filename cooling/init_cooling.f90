module init_cooling_module
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_init_cooling(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INIT_TIME,pst%iUpper+1)
     call r_init_cooling(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call init_cooling(pst%s%mdl,pst%s%r,pst%s%g,pst%s%c)
  endif

end subroutine r_init_cooling
!###########################################################
!###########################################################
!###########################################################
!###########################################################
  subroutine init_cooling(mdl,r,g,c)
  use mdl_module
  use amr_commons, only: run_t,global_t
  use cooling_module, only: cooling_t
  implicit none
  type(mdl_t)::mdl
  type(run_t)::r
  type(global_t)::g
  type(cooling_t)::c

  if(r%cooling.and..not.r%cooling_ism) then
     if(myid==1)write(*,*)'Computing cooling model'
     Nmodel=-1
     if(.not. r%haardt_madau)then
        Nmodel=2
     endif
     if(r%cosmo)then
        ! Reonization redshift has to be later than starting redshift
        r%z_reion=min(1d0/(1.1d0*g%aexp_ini)-1d0,r%z_reion)
        call set_model(Nmodel,dble(J21*1d-21),-1.0d0,dble(a_spec),-1.0d0,dble(z_reion), &
             & -1,2, &
             & dble(h0/100.),dble(omega_b),dble(omega_m),dble(omega_l), &
             & dble(aexp_ini),T2_sim)
        T2_start=T2_sim
        if(nrestart==0)then
           if(myid==1)write(*,*)'Starting with T/mu (K) = ',T2_start
        end if
     else
        call set_model(Nmodel,dble(J21*1d-21),-1.0d0,dble(a_spec),-1.0d0,dble(z_reion), &
             & -1,2, &
             & dble(70./100.),dble(0.04),dble(0.3),dble(0.7), &
             & dble(aexp_ini),T2_sim)
     endif
  end if

end subroutine init_cooling
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module init_cooling_module
