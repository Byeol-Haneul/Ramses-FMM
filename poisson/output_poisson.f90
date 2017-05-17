!#########################################################
!#########################################################
!#########################################################
!#########################################################
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
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine output_poisson_2(r,g,m,filename)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  character(LEN=80)::filename

  integer::ilevel,igrid,ilun
  character(LEN=5)::nchar
  character(LEN=80)::fileloc

  if(r%verbose)write(*,*)'Entering output_poisson'
  ilun=g%ncpu+g%myid+10
  call title(g%myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,access="stream"&
       & ,action="write",form='unformatted')
  write(ilun)ndim
  write(ilun)ndim+1
  write(ilun)r%levelmin
  write(ilun)r%nlevelmax
  do ilevel=r%levelmin,r%nlevelmax
     write(ilun)m%noct(ilevel)
  enddo
#ifdef GRAV
  do ilevel=r%levelmin,r%nlevelmax
     do igrid=m%head(ilevel),m%tail(ilevel)
        write(ilun)m%grid(igrid)%phi
        write(ilun)m%grid(igrid)%f
     end do
  enddo
#endif
  close(ilun)
     
end subroutine output_poisson_2





