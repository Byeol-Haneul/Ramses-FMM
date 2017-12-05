!#########################################################
!#########################################################
!#########################################################
!#########################################################
recursive subroutine r_output_poisson(r,g,m,p,mdl,cpu_range,input_size,output_size,input_array)
  use amr_parameters, only: flen
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size
  integer,dimension(1:input_size)::input_array
  
  integer::next_range,next_cpu
  character(LEN=flen)::filename
  
  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_OUTPUT_POISSON,next_cpu,next_range,input_size,output_size,input_array)
     call r_output_poisson(r,g,m,p,mdl,next_range,input_size,output_size,input_array)
  else
     filename=transfer(input_array,filename)
     call output_poisson(r,g,m,filename)
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





