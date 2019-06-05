module output_hydro_module
contains
!###################################################
!###################################################
!###################################################
!###################################################
recursive subroutine r_output_hydro(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use amr_parameters, only: flen
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array
  
  character(LEN=flen)::filename
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_HYDRO,pst%iUpper+1,input_size,output_size,input_array)
     call r_output_hydro(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     filename=transfer(input_array,filename)
     call output_hydro(pst%s%r,pst%s%g,pst%s%m,filename)
  endif

end subroutine r_output_hydro
!###################################################
!###################################################
!###################################################
!###################################################
subroutine output_hydro(r,g,m,filename)
  use amr_parameters, only: ndim,flen
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  character(LEN=flen)::filename

  integer::ilevel,igrid,ilun
  character(LEN=5)::nchar
  character(LEN=flen)::fileloc

#ifdef HYDRO

  ilun=10!+g%ncpu+g%myid
  call title(g%myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,access="stream",action="write",form='unformatted')
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
     
end subroutine output_hydro
!###################################################
!###################################################
!###################################################
!###################################################
subroutine file_descriptor_hydro(r,filename)
  use amr_parameters, only: ndim,flen
  use hydro_parameters, only: nvar,nener
  use amr_commons, only: run_t
  implicit none
  type(run_t)::r
  character(LEN=flen)::filename

  character(LEN=flen)::fileloc
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

end subroutine file_descriptor_hydro
end module output_hydro_module
