!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from many different initial conditions file formats.
  !--------------------------------------------------------------------
  integer,dimension(1:2)::output_array,dummy

  ! Input particle properties from files
  select case(pst%s%r%filetype)
  case ('grafic')
     call m_input_part_grafic(pst)
  case('ascii')
     call m_input_part_ascii(pst)
  case('gadget')
     write(*,*)'Gadget format not supported yet'
     call mdl_abort(pst%s%mdl)
  case('restart')
     call m_input_part_restart(pst)
  case DEFAULT
     write(*,*) 'Unsupported format file ' // pst%s%r%filetype
     call mdl_abort(pst%s%mdl)
  end select
  
  ! Compute minimum particle mass
  call r_mass_min_part(pst,pst%s%r%levelmin,1,output_array,2)
  
  ! Broadcast minimum particle mass
  call r_broadcast_mp_min(pst,output_array,2,dummy,0)
  
  ! Computing maximum particle count (only in master)
  call r_npart_max(pst,pst%s%r%levelmin,1,pst%s%p%npart_max,1)
  
end subroutine m_input_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_mass_min_part(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer,dimension(1:output_size)::next_output_array
  integer::ilevel
  real(kind=8)::mp_min,next_mp_min

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_MASS_MIN_PART,pst%iUpper+1,input_size,output_size,input_array)
     call r_mass_min_part(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size,next_output_array)
     mp_min=transfer(output_array(1:2),mp_min)
     next_mp_min=transfer(next_output_array(1:2),next_mp_min)
     mp_min=MIN(mp_min,next_mp_min)
     output_array(1:2)=transfer(mp_min,output_array)
  else
     ilevel=input_array(1)
     mp_min=MINVAL(pst%s%p%mp(1:pst%s%p%npart))
     output_array(1:2)=transfer(mp_min,output_array)
  endif

end subroutine r_mass_min_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_broadcast_mp_min(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  real(kind=8)::mp_min
  
  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_BROADCAST_MP_MIN,pst%iUpper+1,input_size,output_size,input_array)
     call r_broadcast_mp_min(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     pst%s%g%mp_min=transfer(input_array(1:2),mp_min)
  endif

end subroutine r_broadcast_mp_min
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_npart_max(pst,ilevel,input_size,npart_max,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::ilevel
  integer,dimension(1:output_size)::npart_max

  integer,dimension(1:output_size)::next_npart_max

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_NPART_MAX,pst%iUpper+1,input_size,output_size,ilevel)
     call r_npart_max(pst%pLower,ilevel,input_size,npart_max,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size,next_npart_max)
     npart_max(1)=MAX(npart_max(1),next_npart_max(1))
  else
     npart_max(1)=pst%s%p%npart
  endif

end subroutine r_npart_max
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
