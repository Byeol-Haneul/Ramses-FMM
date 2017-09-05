program ramses
  use amr_commons
  use pm_commons
  implicit none  

  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  
  ! Read run parameters
  call read_params(r,g)

  ! Start time integration
  call adaptive_loop(r,g,m,p)

end program ramses

