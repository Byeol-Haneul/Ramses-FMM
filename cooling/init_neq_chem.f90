module init_neq_chem_module
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_neq_chem(r,g,tables)
  use amr_commons, only: run_t, global_t
  use neq_cooling_module, only: neq_set_model
  use coolrates_module, only: neq_cooling_t, update_rt_c
  use rtz_cooling_module, only: initialize_elements
  use cosmic_ray_ionization_module, only: initialize_cr_rates
  use photoionization_UVB_module, only: load_UVB_data, update_UVB
  use charge_exchange_module, only: load_ct_rates
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(neq_cooling_t)::tables

  real(kind=8)::T2_sim

  if(g%myid==1)write(*,*)'Computing non-equilibrium chemistry model'

  call update_rt_c(r, g, tables)

  if(r%cosmo)then
     call neq_set_model(r,tables,dble(g%h0/100.),dble(g%omega_b),dble(g%omega_m),dble(g%omega_l),dble(g%aexp_ini),T2_sim)
     g%T2_start=T2_sim
     if(r%nrestart==0)then
        if(g%myid==1)write(*,*)'Starting with T/mu (K) = ',g%T2_start
     end if
  else
     call neq_set_model(r,tables,dble(70./100.),dble(0.049),dble(0.3),dble(0.7),dble(r%aexp_ini),T2_sim)
  endif


#ifdef RTZ
  ! If we are running with RTZ than we have numerous other initializations

  ! First initialize the element structure
  call initialize_elements()

  ! Initialize cosmic ray data
  call initialize_cr_rates()

  ! Initialize the UV background data
  call load_UVB_data() 

  if(r%cosmo)then
     call update_UVB((1.d0/dble(g%aexp_ini))-1.d0)
  else
     call update_UVB((1.d0/dble(r%aexp_ini))-1.d0)
  end if

  !Initialize the charge transfer rates
  call load_ct_rates()

  ! TODO(code): WE NEED THE FOLLOWING INITIALIZATIONS

  !// Initialize high temperature cooling tables
  !initialize_high_temperature_metal_cooling();

  !// Initialize the low temperature cooling tables
  !init_fine_structure_tables();

#endif

end subroutine init_neq_chem
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module init_neq_chem_module
