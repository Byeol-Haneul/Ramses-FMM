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

  integer :: idomain, ilev
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

  !------------------------
  ! Output results to files
  !------------------------
  if(ilevel==levelmin)then
     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout))then
        if(.not.ok_defrag)then
           call defrag
        endif
        call dump_all
     endif
  endif




  if(pic)call make_tree_fine(ilevel)
  
!   if(poisson)then
!      !save old potential for time-extrapolation at level boundaries
!      call save_phi_old(ilevel)

! ! #ifndef WITHOUTMPI
! !      told=MPI_WTIME(info)
! ! #endif

!      call rho_fine(ilevel,icount)     

! ! #ifndef WITHOUTMPI
! !      tnew=MPI_WTIME(info)
! !      dtrho=tnew-told
! !      told=tnew
! ! #endif

! !      if(pic)call hilbert_allparts(ilevel)

! ! #ifndef WITHOUTMPI
! !      tnew=MPI_WTIME(info)
! !      dthilbert=tnew-told
! !      told=tnew
! ! #endif

! !      if (dthilbert>0. .and. dtrho > 0.)then
! !         print*,'dt(hilbert+sort)/dt(rho): ',dthilbert/dtrho,myid,ilevel
! !      end if
!   endif




  call kill_tree_fine(ilevel)



  do idomain=1,ndomain - 1
     bound_key_level(idomain,nlevelmax) = nint(bound_key(idomain) / 8.)
  end do
  bound_key_level(0,nlevelmax) = floor(bound_key(0) / 8.)
  bound_key_level(ndomain,nlevelmax) = ceiling(bound_key(ndomain) / 8.)
  
  do ilev=nlevelmax-1, levelmin, - 1
     do idomain=1,ndomain - 1
        bound_key_level(idomain,ilev) = nint(bound_key_level(idomain, ilev +1) / 8.)
     end do
     bound_key_level(0,ilev) = floor(bound_key_level(0, ilev +1) / 8.)
     bound_key_level(ndomain,ilev) = ceiling(bound_key_level(ndomain, ilev +1) / 8.)
  end do
  


  call MPI_BARRIER(MPI_COMM_WORLD,info)
  if (myid==1)print*, 'starting tests now', npart_andreas

  ! HERE, the grid should be fully built. Do wild test things here...
!  call check_refined_test(test_ok)
!  call check_new_hilbert(test_ok)
  
  do ilev=levelmin, nlevelmax
     told=MPI_WTIME(info)
     call sort_particles(ilev)
     print*, 'ilevel sort: ', ilev, MPI_WTIME(info)-told
     told=MPI_WTIME(info)
     call rho_fine(ilev,icount) 
     call kill_tree_fine(ilev)
     call virtual_tree_fine(ilev)
     print*, 'ilevel rho, kill, virt: ', ilev, MPI_WTIME(info)-told

  end do

end subroutine amr_step


! subroutine check_refined_test(test_ok)
!   use pm_commons
!   use amr_commons, only:levelmin,myid, son, nvector

!   implicit none
! #ifndef WITHOUTMPI
!   include 'mpif.h'
! #endif
!   check whether particles sit in split cells via mpi communication

!   logical::test_ok
!   integer::i, info, nkeys, key_offset,j,dint
!   integer,dimension(1:nvector)::cell_index, cell_levl
!   integer(kind=8),dimension(1:nvector)::hkey2,hkey1,hkey0
!   real(dp)::dreal

!   print*, 'entering'

!   levelp_andreas(1:npart_andreas)=levelmin

!   call hilbert_allparts_andreas(levelmin)

! !  do i=1,npart_andreas-1
! !     if (big_hkey(i,2)==0 .and. big_hkey(i,1)==0 .and. big_hkey(i,0)==0)then
! !        print*,'found zero hilbert key for particle ',i,' on domain ',myid
! !        print*,xp_andreas(i,1:3)
! !     end if
! !  end do


!   ! check if properly sorted in memory
!   do i=1,npart_andreas-1
!      if (big_hkey(i,2) > big_hkey(i+1,2))then
!         test_ok=.false.
! !        print*, big_hkey(i,2), big_hkey(i+1,2)
!      end if
!      if (big_hkey(i,2) == big_hkey(i+1,2) .and. &
!          big_hkey(i,1) > big_hkey(i+1,1))then
!         test_ok=.false.
! !        print*, big_hkey(i,2), big_hkey(i+1,2)
! !        print*, big_hkey(i,1), big_hkey(i+1,1)
!      end if
!      if (big_hkey(i,2) == big_hkey(i+1,2) .and. &
!          big_hkey(i,1) == big_hkey(i+1,1) .and. &
!          big_hkey(i,0) > big_hkey(i+1,0))then
!         test_ok=.false.
! !        print*, big_hkey(i,2), big_hkey(i+1,2)
! !        print*, big_hkey(i,1), big_hkey(i+1,1)
! !        print*, big_hkey(i,0), big_hkey(i+1,0)
!      end if
!   end do

  


!   call build_particle_communicator

!   print*,myid,' send cnt',part_send_cnt
!   print*,myid,' send oft',part_send_oft
!   print*,myid,' recv cnt',part_recv_cnt
!   print*,myid,' recv oft',part_recv_oft

!   received_keys=-1

!   do ipart

! #ifndef WITHOUTMPI                                                                           
!   do i=0,2
!      call MPI_ALLTOALLV(big_hkey(1,i),part_send_cnt,part_send_oft,MPI_INTEGER8, &           
!           &             received_keys(1,i),part_recv_cnt,part_recv_oft,&
!           &             MPI_INTEGER8,MPI_COMM_WORLD,info) 
!   end do
! #endif

!   print*,MINVAL(received_keys)
  
!   do key_offset=0,part_recv_tot-1,nvector
!      nkeys=min(part_recv_tot-key_offset,nvector)
!      hkey2(1:nkeys)=received_keys(key_offset+1:key_offset+nkeys,2)
!      hkey1(1:nkeys)=received_keys(key_offset+1:key_offset+nkeys,1)
!      hkey0(1:nkeys)=received_keys(key_offset+1:key_offset+nkeys,0)
!      call get_cell_index_from_hilbertkey(cell_index,cell_levl,hkey2,hkey1,hkey0,nkeys,levelmin)
!      received_keys(key_offset+1:key_offset+nkeys,0)=son(cell_index(1:nkeys))
!   end do
  
! #ifndef WITHOUTMPI                                                                           
!   call MPI_ALLTOALLV(received_keys(1,0),part_recv_cnt,part_recv_oft,MPI_INTEGER8, &           
!        &             part_ref_mask,part_send_cnt,part_send_oft,&
!        &             MPI_INTEGER8,MPI_COMM_WORLD,info) 
! #endif

!   do i=1,npart_andreas
!      if(part_ref_mask(i)>0)then
!         part_ref_mask(i)=1        
!         levelp_andreas(i)=levelp_andreas(i)+1
!      end if
!   end do

!   do j=1,10
!      call random_number(dreal)
!      dint= int(dreal*99999)+1
!      do i=1,npart_andreas
!         if (idp_andreas(i)==dint)then
!            print*,dint,'++++++++++++++++andreas'
!            print*,levelp_andreas(i)
!            print*,dint,'++++++++++++++++andreas'
!         end if
!      end do
     
!      do i=1,npart
!         if (idp(i)==dint)then
!            print*,dint,'+++++++++++++++++++++++++'
!            print*,levelp(i)
!            print*,dint,'+++++++++++++++++++++++++'
!         end if
!      end do
!   end do


! end subroutine check_refined_test



subroutine check_new_hilbert(test_ok)
  use pm_commons
  use amr_commons
  use sort, only: msd_radix_sort_particles, lsd_radix_sort_particles
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  logical::test_ok
  integer::i, info, nkeys, key_offset,j,dint
  integer(kind=8), dimension(1:npart_andreas) :: store_key
!  real(dp), dimension(1:npart_andreas,1:3) :: store_pos
  real(dp)::told
  integer ::  offset, np, ipart

  print*, 'entering'

  levelp_andreas(1:npart_andreas)=levelmin

  told=MPI_WTIME(info)
  call hilbert3d_for_particle(0,npart_andreas, 0, levelmin)
  print*, 'new hilbert: ',MPI_WTIME(info)-told


  store_key(1:npart_andreas) = part_hkey(1:npart_andreas, 0)

  
  told=MPI_WTIME(info)
  part_hkey(1:npart_andreas, 0) = store_key(1:npart_andreas)
  call lsd_radix_sort_particles(0, npart_andreas, levelmin, levelmin)
  print*, 'lsd  sort: ', MPI_WTIME(info)-told

  do i=1,npart_andreas-1
     if (part_hkey(i,0) > part_hkey(i+1,0))then
        print*,part_hkey(i,0),part_hkey(i+1,0)
        test_ok=.false.
     end if
  end do




  
  call MPI_BARRIER(MPI_COMM_WORLD,info)

end subroutine check_new_hilbert

subroutine sort_particles(ilevel)
  use pm_commons
  use amr_commons
  use sort, only:msd_radix_sort_particles, lsd_radix_sort_particles
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in) :: ilevel

  logical::test_ok
  integer::i, info, nkeys, key_offset,j,dint
  real(dp)::told
  integer ::  offset, np, ipart
  
  offset = part_level_offset(ilevel)
  if (ilevel < nlevelmax - 1)then
     np = part_level_offset(ilevel+2) - part_level_offset(ilevel)  
  else
     np = npart_andreas - part_level_offset(ilevel)
  end if

  if (ilevel == levelmin)then
     told=MPI_WTIME(info)
     call hilbert3d_for_particle(offset, npart_andreas-offset, 0, levelmin)


     told=MPI_WTIME(info)
     call lsd_radix_sort_particles(offset, np, levelmin, levelmin)


  else

     told=MPI_WTIME(info)
     call hilbert3d_for_particle(offset, npart_andreas-offset, ilevel-1, ilevel)


     told=MPI_WTIME(info)
     call msd_radix_sort_particles(offset, np, ilevel-1, ilevel, ilevel)

  end if


  told=MPI_WTIME(info)
  call compute_particle_histogram(offset, np)


  told=MPI_WTIME(info)
  call build_histogram_communicator(levelmin)
  call send_histogram_bins(offset, np, levelmin)


end subroutine sort_particles
