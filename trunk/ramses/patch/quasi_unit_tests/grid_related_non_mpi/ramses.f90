program ramses
  implicit none  
  logical::test_ok
  ! Read run parameters
  call read_params

  ! Start time integration
  call adaptive_loop(test_ok)

  if (test_ok)then
     print*, 'passed'
  else
     print*, 'failed'
  end if
end program ramses

