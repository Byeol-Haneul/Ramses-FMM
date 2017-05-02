program ramses
  use amr_commons
  use pm_commons
  implicit none  

  type(run_t)::run_p
  type(global_t)::global_v
  type(mesh_t)::mesh_v
  type(part_t)::part_v
  
  ! Read run parameters
  call read_params(run_p,global_v)

!  call run_p%print

  ! Start time integration
  call adaptive_loop(run_p,global_v,mesh_v,part_v)

end program ramses

