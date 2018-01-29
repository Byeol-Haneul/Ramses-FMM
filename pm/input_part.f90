!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part(s)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from many different initial conditions file formats.
  !--------------------------------------------------------------------
  integer,dimension(1:2)::output_array

  ! Input particle properties from files
  select case(s%r%filetype)
  case ('grafic')
     call m_input_part_grafic(s)
  case('ascii')
     call m_input_part_ascii(s)
  case('gadget')
     write(*,*)'Gadget format not supported yet'
     call mdl_abort
  case('restart')
     call m_input_part_restart(s)
  case DEFAULT
     write(*,*) 'Unsupported format file ' // s%r%filetype
     call mdl_abort
  end select
  
  ! Compute minimum particle mass
  call r_mass_min_part(s,s%g%ncpu,1,2,s%r%levelmin,output_array)
  
  ! Broadcast minimum particle mass
  call r_broadcast_mp_min(s,s%g%ncpu,2,0,output_array)
  
  ! Computing maximum particle count (only in master)
  call r_npart_max(s,s%g%ncpu,1,1,s%r%levelmin,s%p%npart_max)
  
end subroutine m_input_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_mass_min_part(s,cpu_range,input_size,output_size,ilevel,output_array)
  use ramses_commons, only: ramses_t
  use mdl_parameters
  implicit none
  type(ramses_t)::s
  integer::cpu_range,input_size,output_size
  integer::ilevel
  integer,dimension(1:output_size)::output_array

  integer::next_range,next_cpu
  integer,dimension(1:output_size)::next_output_array
  real(kind=8)::mp_min,next_mp_min

  next_range=cpu_range/2
  next_cpu=s%g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(s%mdl,MDL_MASS_MIN_PART,next_cpu,next_range,input_size,output_size,ilevel)
     call r_mass_min_part(s,next_range,input_size,output_size,ilevel,output_array)
     call mdl_get_reply(s%mdl,next_cpu,output_size,next_output_array)
     mp_min=transfer(output_array(1:2),mp_min)
     next_mp_min=transfer(next_output_array(1:2),next_mp_min)
     mp_min=MIN(mp_min,next_mp_min)
     output_array(1:2)=transfer(mp_min,output_array)
  else
     mp_min=MINVAL(s%p%mp(1:s%p%npart))
     output_array(1:2)=transfer(mp_min,output_array)
  endif

end subroutine r_mass_min_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_broadcast_mp_min(s,cpu_range,input_size,output_size,input_array)
  use ramses_commons, only: ramses_t
  use mdl_parameters
  implicit none
  type(ramses_t)::s
  integer::cpu_range,input_size,output_size
  integer,dimension(1:input_size)::input_array

  integer::next_range,next_cpu
  real(kind=8)::mp_min

  next_range=cpu_range/2
  next_cpu=s%g%myid+next_range
  
  if(next_range>0)then
     call mdl_send_request(s%mdl,MDL_BROADCAST_MP_MIN,next_cpu,next_range,input_size,output_size,input_array)
     call r_broadcast_mp_min(s,next_range,input_size,output_size,input_array)
  else
     s%g%mp_min=transfer(input_array(1:2),mp_min)
  endif

end subroutine r_broadcast_mp_min
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_npart_max(s,cpu_range,input_size,output_size,ilevel,npart_max)
  use ramses_commons, only: ramses_t
  use mdl_parameters
  implicit none
  type(ramses_t)::s
  integer::cpu_range,input_size,output_size
  integer::ilevel,npart_max

  integer::next_range,next_cpu
  integer::next_npart_max

  next_range=cpu_range/2
  next_cpu=s%g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(s%mdl,MDL_NPART_MAX,next_cpu,next_range,input_size,output_size,ilevel)
     call r_npart_max(s,next_range,input_size,output_size,ilevel,npart_max)
     call mdl_get_reply(s%mdl,next_cpu,output_size,next_npart_max)
     npart_max=MAX(npart_max,next_npart_max)
  else
     npart_max=s%p%npart
  endif

end subroutine r_npart_max
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
