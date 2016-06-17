subroutine output_poisson(filename)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename

  integer::ilevel,igrid,ilun
  character(LEN=5)::nchar
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering output_poisson'
  ilun=ncpu+myid+10
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,access="stream"&
       & ,action="write",form='unformatted')
  write(ilun)ndim
  write(ilun)ndim+1
  write(ilun)levelmin
  write(ilun)nlevelmax
  do ilevel=levelmin,nlevelmax
     write(ilun)noct(ilevel)
  enddo
#ifdef GRAV
  do ilevel=levelmin,nlevelmax
     do igrid=head(ilevel),tail(ilevel)
        write(ilun)grid(igrid)%phi
        write(ilun)grid(igrid)%f
     end do
  enddo
#endif
  close(ilun)
     
end subroutine output_poisson





