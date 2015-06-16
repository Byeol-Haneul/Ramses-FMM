program ramses_unit_tests
  implicit none  
  logical::all_ok=.true.
  integer::i


  call amr_tests(all_ok)
!  call hydro_tests(all_ok)
  call pm_tests(all_ok)
!  call poisson_tests(all_ok)

  
  if (all_ok)then
     write(*,*)'ALL UNIT TESTS PASSED SUCCESSFULLY'
  else
     write(*,*)'ATTENTION: TEST FAILURES OCCURED'
  end if

end program ramses_unit_tests

