subroutine output_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
  character(LEN=80)::filename

  integer::ilevel,igrid,ilun
  character(LEN=5)::nchar
  character(LEN=80)::fileloc

#ifdef HYDRO

  if(verbose)write(*,*)'Entering output_hydro'
  ilun=ncpu+myid+10     
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,access="stream"&
       & ,action="write",form='unformatted')
  write(ilun)ndim
  write(ilun)nvar
  write(ilun)levelmin
  write(ilun)nlevelmax
  do ilevel=levelmin,nlevelmax
     write(ilun)noct(ilevel)
  enddo
  do ilevel=levelmin,nlevelmax
     do igrid=head(ilevel),tail(ilevel)
        write(ilun)grid(igrid)%uold
     end do
  enddo
  close(ilun)

#endif
     
end subroutine output_hydro
!###################################################
!###################################################
!###################################################
!###################################################
subroutine output_hydro_2(r,g,m,filename)
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

#ifdef HYDRO

  if(r%verbose)write(*,*)'Entering output_hydro'
  ilun=g%ncpu+g%myid+10     
  call title(g%myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,access="stream"&
       & ,action="write",form='unformatted')
  write(ilun)ndim
  write(ilun)nvar
  write(ilun)r%levelmin
  write(ilun)r%nlevelmax
  do ilevel=r%levelmin,r%nlevelmax
     write(ilun)m%noct(ilevel)
  enddo
  do ilevel=r%levelmin,r%nlevelmax
     do igrid=m%head(ilevel),m%tail(ilevel)
        write(ilun)m%grid(igrid)%uold
     end do
  enddo
  close(ilun)

#endif
     
end subroutine output_hydro_2
!###################################################
!###################################################
!###################################################
!###################################################
subroutine file_descriptor_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  character(LEN=80)::filename
  character(LEN=80)::fileloc
  integer::ivar,ilun

  if(verbose)write(*,*)'Entering file_descriptor_hydro'

  ilun=11

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write run parameters
  write(ilun,'("nvar        =",I11)')nvar
  ivar=1
  write(ilun,'("variable #",I2,": density")')ivar
  ivar=2
  write(ilun,'("variable #",I2,": velocity_x")')ivar
  if(ndim>1)then
     ivar=3
     write(ilun,'("variable #",I2,": velocity_y")')ivar
  endif
  if(ndim>2)then
     ivar=4
     write(ilun,'("variable #",I2,": velocity_z")')ivar
  endif
#if NENER>0
  ! Non-thermal pressures
  do ivar=ndim+2,ndim+1+nener
     write(ilun,'("variable #",I2,": non_thermal_pressure_",I1)')ivar,ivar-ndim-1
  end do
#endif
  ivar=ndim+2+nener
  write(ilun,'("variable #",I2,": thermal_pressure")')ivar
#if NVAR>NDIM+2+NENER
  ! Passive scalars
  do ivar=ndim+3+nener,nvar
     write(ilun,'("variable #",I2,": passive_scalar_",I1)')ivar,ivar-ndim-2-nener
  end do
#endif

  close(ilun)

end subroutine file_descriptor_hydro
!###################################################
!###################################################
!###################################################
!###################################################
subroutine file_descriptor_hydro_2(r,g,filename)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nvar,nener
  use amr_commons, only: run_t,global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  character(LEN=80)::filename

  character(LEN=80)::fileloc
  integer::ivar,ilun

  if(r%verbose)write(*,*)'Entering file_descriptor_hydro'

  ilun=11

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write run parameters
  write(ilun,'("nvar        =",I11)')nvar
  ivar=1
  write(ilun,'("variable #",I2,": density")')ivar
  ivar=2
  write(ilun,'("variable #",I2,": velocity_x")')ivar
  if(ndim>1)then
     ivar=3
     write(ilun,'("variable #",I2,": velocity_y")')ivar
  endif
  if(ndim>2)then
     ivar=4
     write(ilun,'("variable #",I2,": velocity_z")')ivar
  endif
#if NENER>0
  ! Non-thermal pressures
  do ivar=ndim+2,ndim+1+nener
     write(ilun,'("variable #",I2,": non_thermal_pressure_",I1)')ivar,ivar-ndim-1
  end do
#endif
  ivar=ndim+2+nener
  write(ilun,'("variable #",I2,": thermal_pressure")')ivar
#if NVAR>NDIM+2+NENER
  ! Passive scalars
  do ivar=ndim+3+nener,nvar
     write(ilun,'("variable #",I2,": passive_scalar_",I1)')ivar,ivar-ndim-2-nener
  end do
#endif

  close(ilun)

end subroutine file_descriptor_hydro_2





