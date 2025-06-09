module init_cooling_module
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
  subroutine init_cooling(r,g,c)
  use amr_commons, only: run_t,global_t
  use cooling_module, only: cooling_t, set_model
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(cooling_t)::c

  integer::Nmodel
  real(kind=8)::T2_sim

  if(g%myid==1)write(*,*)'Computing cooling model'
  Nmodel=-1
  if(.not. r%haardt_madau)then
     Nmodel=2
  endif
  if(r%cosmo)then
     ! Reonization redshift has to be later than starting redshift
     r%z_reion=min(1d0/(1.1d0*dble(g%aexp_ini))-1d0,dble(r%z_reion))
     call set_model(c,Nmodel,dble(r%J21*1d-21),-1.0d0,dble(r%a_spec),-1.0d0,dble(r%z_reion), &
          & -1,2, &
          & dble(g%h0/100.),dble(g%omega_b),dble(g%omega_m),dble(g%omega_l), &
          & dble(g%aexp_ini),T2_sim, r%mu_mol)
     g%T2_start=T2_sim
     if(r%nrestart==0)then
        if(g%myid==1)write(*,*)'Starting with T/mu (K) = ',g%T2_start
     end if
  else
     call set_model(c,Nmodel,dble(r%J21*1d-21),-1.0d0,dble(r%a_spec),-1.0d0,dble(r%z_reion), &
          & -1,2, &
          & dble(70./100.),dble(0.04),dble(0.3),dble(0.7), &
          & dble(r%aexp_ini),T2_sim, r%mu_mol)
  endif

end subroutine init_cooling
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module init_cooling_module
