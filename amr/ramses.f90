program ramses
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
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

