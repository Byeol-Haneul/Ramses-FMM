!#######################################################
!#######################################################
!#######################################################
!#######################################################
recursive subroutine r_output_part(r,g,m,p,mdl,cpu_range,input_size,output_size,input_array)
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
     call mdl_send_request(mdl,MDL_OUTPUT_PART,next_cpu,next_range,input_size,output_size,input_array)
     call r_output_part(r,g,m,p,mdl,next_range,input_size,output_size,input_array)
  else
     filename=transfer(input_array,filename)
     call output_part(r,g,p,filename)
  endif

end subroutine r_output_part
!#######################################################
!#######################################################
!#######################################################
!#######################################################
subroutine output_part(r,g,p,filename)
  use amr_parameters, only: ndim,dp,i8b,flen
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  character(LEN=flen)::filename

  integer::i,idim,ilun
  character(LEN=flen)::fileloc
  character(LEN=5)::nchar
  real(dp),allocatable,dimension(:)::xdp
  integer(i8b),allocatable,dimension(:)::ii8
  integer,allocatable,dimension(:)::ll
  
  ilun=2*g%ncpu+g%myid+10

  call title(g%myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=TRIM(fileloc),access="stream",action="write",form='unformatted')
  rewind(ilun)
  ! Write header
  write(ilun)ndim
  write(ilun)p%npart
  ! Write position
  allocate(xdp(1:p%npart))
  do idim=1,ndim
     do i=1,p%npart
        xdp(i)=p%xp(i,idim)
     end do
     write(ilun)xdp
  end do
  ! Write velocity
  do idim=1,ndim
     do i=1,p%npart
        xdp(i)=p%vp(i,idim)
     end do
     write(ilun)xdp
  end do
  ! Write mass
  do i=1,p%npart
     xdp(i)=p%mp(i)
  end do
  write(ilun)xdp
  deallocate(xdp)
  ! Write identity
  allocate(ii8(1:p%npart))
  do i=1,p%npart
     ii8(i)=p%idp(i)
  end do
  write(ilun)ii8
  deallocate(ii8)
  ! Write level
  allocate(ll(1:p%npart))
  do i=1,p%npart
     ll(i)=p%levelp(i)
  end do
  write(ilun)ll
  deallocate(ll)

#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Write potential (optional)
  allocate(xdp(1:p%npart))
  do i=1,p%npart
     xdp(i)=p%phip(i)
  end do
  write(ilun)xdp
  deallocate(xdp)
#endif

  close(ilun)

end subroutine output_part

