!#########################################################
!#########################################################
!#########################################################
!#########################################################
recursive subroutine r_output_poisson(pst,input_size,output_size,input_array)
  use amr_parameters, only: flen
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  
  character(LEN=flen)::filename
  
  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_OUTPUT_POISSON,pst%iUpper+1,input_size,output_size,input_array)
     call r_output_poisson(pst%pLower,input_size,output_size,input_array)
  else
     filename=transfer(input_array,filename)
     call output_poisson(pst%s%r,pst%s%g,pst%s%m,filename)
  endif
end subroutine r_output_poisson
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine output_poisson(r,g,m,filename)
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

  ilun=g%ncpu+g%myid+10
  call title(g%myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
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
     
end subroutine output_poisson





