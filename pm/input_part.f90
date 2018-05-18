!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part(pst)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from many different initial conditions file formats.
  !--------------------------------------------------------------------
  integer,dimension(1:2)::output_array

  ! Input particle properties from files
  select case(pst%s%r%filetype)
  case ('grafic')
     call m_input_part_grafic(pst)
  case('ascii')
     call m_input_part_ascii(pst)
  case('gadget')
     write(*,*)'Gadget format not supported yet'
     call mdl_abort
  case('restart')
     call m_input_part_restart(pst)
  case DEFAULT
     write(*,*) 'Unsupported format file ' // pst%s%r%filetype
     call mdl_abort
  end select
  
  ! Compute minimum particle mass
  call r_mass_min_part(pst,1,2,pst%s%r%levelmin,output_array)
  
  ! Broadcast minimum particle mass
  call r_broadcast_mp_min(pst,2,0,output_array)
  
  ! Computing maximum particle count (only in master)
  call r_npart_max(pst,1,1,pst%s%r%levelmin,pst%s%p%npart_max)
  
end subroutine m_input_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_mass_min_part(pst,input_size,output_size,ilevel,output_array)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer::ilevel
  integer,dimension(1:output_size)::output_array

  integer,dimension(1:output_size)::next_output_array
  real(kind=8)::mp_min,next_mp_min

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_MASS_MIN_PART,pst%iUpper+1,input_size,output_size,ilevel)
     call r_mass_min_part(pst%pLower,input_size,output_size,ilevel,output_array)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size,next_output_array)
     mp_min=transfer(output_array(1:2),mp_min)
     next_mp_min=transfer(next_output_array(1:2),next_mp_min)
     mp_min=MIN(mp_min,next_mp_min)
     output_array(1:2)=transfer(mp_min,output_array)
  else
     mp_min=MINVAL(pst%s%p%mp(1:pst%s%p%npart))
     output_array(1:2)=transfer(mp_min,output_array)
  endif

end subroutine r_mass_min_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_broadcast_mp_min(pst,input_size,output_size,input_array)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array

  real(kind=8)::mp_min
  
  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_BROADCAST_MP_MIN,pst%iUpper+1,input_size,output_size,input_array)
     call r_broadcast_mp_min(pst%pLower,input_size,output_size,input_array)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     pst%s%g%mp_min=transfer(input_array(1:2),mp_min)
  endif

end subroutine r_broadcast_mp_min
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_npart_max(pst,input_size,output_size,ilevel,npart_max)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer::ilevel,npart_max

  integer::next_npart_max

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_NPART_MAX,pst%iUpper+1,input_size,output_size,ilevel)
     call r_npart_max(pst%pLower,input_size,output_size,ilevel,npart_max)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size,next_npart_max)
     npart_max=MAX(npart_max,next_npart_max)
  else
     npart_max=pst%s%p%npart
  endif

end subroutine r_npart_max
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
