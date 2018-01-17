!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part(r,g,m,p,mdl)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from many different initial conditions file formats.
  !--------------------------------------------------------------------
  integer,dimension(1:2)::output_array

  ! Input particle properties from files
  select case(r%filetype)
  case ('grafic')
     call m_input_part_grafic(r,g,m,p,mdl)
  case('ascii')
     call m_input_part_ascii(r,g,m,p,mdl)
  case('gadget')
     write(*,*)'Gadget format not supported yet'
     call mdl_abort
  case('restart')
     call m_input_part_restart(r,g,m,p,mdl)
  case DEFAULT
     write(*,*) 'Unsupported format file ' // r%filetype
     call mdl_abort
  end select
  
  ! Compute minimum particle mass
  call r_mass_min_part(r,g,m,p,mdl,g%ncpu,1,2,r%levelmin,output_array)
  
  ! Broadcast minimum particle mass
  call r_broadcast_mp_min(r,g,m,p,mdl,g%ncpu,2,0,output_array)
  
  ! Computing maximum particle count (only in master)
  call r_npart_max(r,g,m,p,mdl,g%ncpu,1,1,r%levelmin,p%npart_max)
  
end subroutine m_input_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_mass_min_part(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel,output_array)
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
  integer::ilevel
  integer,dimension(1:output_size)::output_array

  integer::next_range,next_cpu
  integer,dimension(1:output_size)::next_output_array
  real(kind=8)::mp_min,next_mp_min

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_MASS_MIN_PART,next_cpu,next_range,input_size,output_size,ilevel)
     call r_mass_min_part(r,g,m,p,mdl,next_range,input_size,output_size,ilevel,output_array)
     call mdl_get_reply(mdl,next_cpu,output_size,next_output_array)
     mp_min=transfer(output_array(1:2),mp_min)
     next_mp_min=transfer(next_output_array(1:2),next_mp_min)
     mp_min=MIN(mp_min,next_mp_min)
     output_array(1:2)=transfer(mp_min,output_array)
  else
     mp_min=MINVAL(p%mp(1:p%npart))
     output_array(1:2)=transfer(mp_min,output_array)
  endif

end subroutine r_mass_min_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_broadcast_mp_min(r,g,m,p,mdl,cpu_range,input_size,output_size,input_array)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  use hilbert
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size
  integer,dimension(1:input_size)::input_array

  integer::next_range,next_cpu
  real(kind=8)::mp_min

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_BROADCAST_MP_MIN,next_cpu,next_range,input_size,output_size,input_array)
     call r_broadcast_mp_min(r,g,m,p,mdl,next_range,input_size,output_size,input_array)
  else
     g%mp_min=transfer(input_array(1:2),mp_min)
  endif

end subroutine r_broadcast_mp_min
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_npart_max(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel,npart_max)
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
  integer::ilevel,npart_max

  integer::next_range,next_cpu
  integer::next_npart_max

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_NPART_MAX,next_cpu,next_range,input_size,output_size,ilevel)
     call r_npart_max(r,g,m,p,mdl,next_range,input_size,output_size,ilevel,npart_max)
     call mdl_get_reply(mdl,next_cpu,output_size,next_npart_max)
     npart_max=MAX(npart_max,next_npart_max)
  else
     npart_max=p%npart
  endif

end subroutine r_npart_max
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
