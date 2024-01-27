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
        call output_poisson(pst%s%r,pst%s%g,pst%s%m,pst%s%mdl,input%filename)
     endif
  endif
end subroutine r_output_poisson
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine output_poisson(r,g,m,mdl,filename)
  use amr_parameters, only: ndim,twotondim,flen
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
  real(kind=4),dimension(1:twotondim,1:ndim)::f
  real(kind=4),dimension(1:twotondim)::phi

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
        phi=real(m%grid(igrid)%phi,kind=4)
        f(1:twotondim,1:ndim)=real(m%grid(igrid)%f(1:twotondim,1:ndim),kind=4)
        write(ilun)phi
        write(ilun)f
     end do
  enddo
#endif
  close(ilun)

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
!#########################################################
!#########################################################
!#########################################################
!#########################################################
end module output_poisson_module
