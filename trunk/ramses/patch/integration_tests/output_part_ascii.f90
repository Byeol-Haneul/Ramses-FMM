subroutine part_to_ascii
  use amr_commons,   only: myid
  use pm_commons,    only: npart, xp, vp, idp
  use pm_parameters, only: npartmax
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  
  ! Small routine which produces one ascii file containing all the particle
  ! positions in the entire box, ordered by the particle index idp.
  ! This file can easily be diffed to compare results.
  ! VERY SLOW, only use for debug or test puroses.

  integer :: i, j, info

  ! create the file
  if (myid == 1)then
     open(12, file="particles.txt", action="write", form="formatted")
     write(12,'(A)'), "particle locations and positions"
     close(12)
  end if
     
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif

  do i = 1, npartmax
     do j = 1, npartmax
        if (idp(j) == i)then
           open(12, file="particles.txt", status="old", position="append", action="write", form="formatted")
           write(12, '(A,X,I10,3(X,E16.8E2))'),"xp:", i, xp(j, :)
           write(12, '(A,X,I10,3(X,E16.8E2))'),"vp:", i, vp(j, :)
           close(12)
        end if
     end do

#ifndef WITHOUTMPI
     call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif
  end do
  
end subroutine part_to_ascii
