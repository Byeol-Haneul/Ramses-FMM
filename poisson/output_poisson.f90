module output_poisson_module
  use amr_parameters, only: flen
  type :: in_output_poisson_t
    character(LEN=flen)::filename
  end type in_output_poisson_t
contains
!#########################################################
!#########################################################
!#########################################################
!#########################################################
recursive subroutine r_output_poisson(pst,input,input_size)
  use mdl_module
  use amr_parameters, only: flen
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(in_output_poisson_t)::input
  
  integer::rID
  
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_POISSON,pst%iUpper+1,input_size,0,input)
     call r_output_poisson(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     if(index(input%filename,'output')==0)then
        call backup_poisson(pst%s%r,pst%s%g,pst%s%m,pst%s%mdl,input%filename)
     else
        call output_poisson(pst%s,input%filename)
     endif
  endif
end subroutine r_output_poisson
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine output_poisson(s,filename)
  use amr_parameters, only: ndim,twotondim,flen
  use ramses_commons, only: ramses_t,open_file,close_file
  use mdl_module
  implicit none
  type(ramses_t)::s
  character(LEN=flen)::filename
  !-----------------------------------
  ! Output grav data in file
  !-----------------------------------
  integer::ilevel,igrid,ilun
  integer(kind=8),dimension(s%r%levelmin:s%r%nlevelmax)::nskip
  real(kind=4),dimension(1:twotondim,1:ndim)::f
  real(kind=4),dimension(1:twotondim)::phi

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

#ifdef GRAV

  call open_file(s,filename,nskip,ilun)

  do ilevel=r%levelmin,r%nlevelmax
     write(ilun,POS=nskip(ilevel))
     do igrid=m%head(ilevel),m%tail(ilevel)
        phi=real(m%grid(igrid)%phi,kind=4)
        f(1:twotondim,1:ndim)=real(m%grid(igrid)%f(1:twotondim,1:ndim),kind=4)
        write(ilun)phi
        write(ilun)f
     end do
  enddo

  call close_file(s,filename,nskip,ilun)

#endif

  end associate

end subroutine output_poisson
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine backup_poisson(r,g,m,mdl,filename)
  use amr_parameters, only: ndim,flen
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  use mdl_module
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(mdl_t)::mdl
  character(LEN=flen)::filename

  integer::ilevel,igrid,ilun,ierr
  character(LEN=5)::nchar
  character(LEN=flen)::fileloc
  logical::file_exist

  ilun=10+mdl_core(mdl)
  call title(g%myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  inquire(file=fileloc, exist=file_exist)
  if (file_exist) then
     open(unit=ilun,file=fileloc,iostat=ierr)
     close(ilun,status="delete")
  end if
  open(unit=ilun,file=fileloc,access="stream",action="write",form='unformatted')
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
     
end subroutine backup_poisson
!###################################################
!###################################################
!###################################################
!###################################################
subroutine file_descriptor_poisson(r,filename,write_bkp_file)
  use amr_parameters, only: ndim,flen
  use hydro_parameters, only: nvar,nener
  use amr_commons, only: run_t
  implicit none
  type(run_t)::r
  character(LEN=flen)::filename
  logical::write_bkp_file

  character(LEN=flen)::fileloc
  integer::ivar,ilun

  if(r%verbose)write(*,*)'Entering file_descriptor_poisson'

  ilun=11

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write variable names in backup file or in output file
  write(ilun,'("nvar        =",I11)')ndim+1
  ivar=1
  write(ilun,'("variable #",I2,": potential")')ivar
  ivar=2
  write(ilun,'("variable #",I2,": accel_x")')ivar
  if(ndim>1)then
     ivar=3
     write(ilun,'("variable #",I2,": accel_y")')ivar
  endif
  if(ndim>2)then
     ivar=4
     write(ilun,'("variable #",I2,": accel_z")')ivar
  endif

  close(ilun)

end subroutine file_descriptor_poisson
!#########################################################
!#########################################################
!#########################################################
!#########################################################
end module output_poisson_module
