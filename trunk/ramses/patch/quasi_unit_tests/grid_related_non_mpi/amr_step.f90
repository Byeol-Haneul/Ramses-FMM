subroutine amr_step(ilevel,icount,test_ok)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  logical::test_ok
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  integer::i,idim,ivar,info
  logical::ok_defrag
  logical,save::first_step=.true.
  real(dp)::told,tnew,dthilbert,dtrho

  test_ok=.true.

  if(numbtot(1,ilevel)==0)return


  !-------------------------------------------
  ! Make new refinements and update boundaries
  !-------------------------------------------
  if(levelmin.lt.nlevelmax .and..not. static)then
     if(ilevel==levelmin.or.icount>1)then
        do i=ilevel,nlevelmax
           if(i>levelmin)then

              !--------------------------
              ! Build communicators
              !--------------------------
              call build_comm(i)

              !--------------------------
              ! Update boundaries
              !--------------------------
              call make_virtual_fine_int(cpu_map(1),i)
              if(hydro)then
                 do ivar=1,nvar
                    call make_virtual_fine_dp(uold(1,ivar),i)
                 end do
                 if(simple_boundary)call make_boundary_hydro(i)
              end if
              if(poisson)then
                 call make_virtual_fine_dp(phi(1),i)
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),i)
                 end do
                 if(simple_boundary)call make_boundary_force(i)
              end if
           end if

           !--------------------------
           ! Refine grids
           !--------------------------
           call refine_fine(i)
        end do
     end if
  end if

  !--------------------------
  ! Load balance
  !--------------------------
  ok_defrag=.false.
  if(levelmin.lt.nlevelmax)then
     if(ilevel==levelmin)then
        if(nremap>0)then
           ! Skip first load balance because 
           ! it has been performed before file dump
           if(nrestart>0.and.first_step)then
              first_step=.false.
           else
              if(MOD(nstep_coarse,nremap)==0)then
                 call load_balance
                 call defrag
                 ok_defrag=.true.
              endif
           end if
        end if
     endif
  end if



  ! HERE, the grid should be fully built. Do wild test things here...
  call hilbert_to_cellindex_test(test_ok)


end subroutine amr_step


subroutine hilbert_to_cellindex_test(test_ok)
  use amr_parameters, only:dp, nvector
  use amr_commons
  implicit none
  
  logical::test_ok
  real(dp),dimension(1:nvector,1:3)::xx
  integer,dimension(1:nvector)::cell_index1,cell_index2,cell_levl1,cell_levl2
  integer(kind=8),dimension(1:nvector)::hkey2,hkey1,hkey0
  integer::i, j, ilevel

  do j=1,100
     do ilevel=levelmin,nlevelmax
        call random_number(xx(1:nvector,1))
        call random_number(xx(1:nvector,2))
        call random_number(xx(1:nvector,3))
        
        call get_cell_index(cell_index1,cell_levl1,xx,ilevel,nvector)
        call cmp_ordering_int(xx,hkey2,hkey1,hkey0,nvector)     
        call get_cell_index_from_hilbertkey(cell_index2,cell_levl2,hkey2,hkey1,hkey0,nvector,ilevel)
        
        
        do i=1,nvector
           test_ok=test_ok .and. cell_index2(i)==cell_index1(i)
           test_ok=test_ok .and. cell_levl2(i)==cell_levl1(i)
        end do
     end do
  end do


end subroutine hilbert_to_cellindex_test
