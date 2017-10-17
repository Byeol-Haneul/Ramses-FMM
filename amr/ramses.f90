program ramses
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons
  implicit none  

  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  
  integer,dimension(:),allocatable::input
  type(run_t)::r2

  ! Initialisation of MDL
  call mdl_init(r,g,m,p)

  ! Read run parameters
  call read_params(r,g,m,p)

!!$  call r%print
!!$
!!$  write(*,*)storage_size(r)/32,MDL_INPUT_MAXSIZE
!!$
!!$  allocate(input(1:MDL_INPUT_MAXSIZE))
!!$  
!!$  input=transfer(r,input)
!!$
!!$  write(*,*)input
!!$
!!$  r2=transfer(input,r)
!!$
!!$  deallocate(input)
!!$
!!$  call r2%print
  stop


  ! Start time integration
  call adaptive_loop(r,g,m,p)

end program ramses

