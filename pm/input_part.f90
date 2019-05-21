module input_part_module
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use input_part_grafic_module, only: m_input_part_grafic
  use input_part_ascii_module, only: m_input_part_ascii
  use input_part_restart_module, only: m_input_part_restart
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from many different initial conditions file formats.
  !--------------------------------------------------------------------
  integer,dimension(1:2)::output_array,dummy
  real(kind=8)::mp_min

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
  call r_mass_min_part(pst,pst%s%r%levelmin,1,mp_min,2)
  
  ! Broadcast minimum particle mass
  call r_broadcast_mp_min(pst,mp_min,2,dummy,0)
  
  ! Computing maximum particle count (only in master)
  call r_npart_max(pst,pst%s%r%levelmin,1,pst%s%p%npart_max,1)
  
end subroutine m_input_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_mass_min_part(pst,ilevel,input_size,mp_min,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer::ilevel
  real(kind=8)::mp_min,next_mp_min
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_MASS_MIN_PART,pst%iUpper+1,input_size,output_size,ilevel)
     call r_mass_min_part(pst%pLower,ilevel,input_size,mp_min,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_mp_min)
     mp_min=MIN(mp_min,next_mp_min)
  else
     mp_min=MINVAL(pst%s%p%mp(1:pst%s%p%npart))
  endif

end subroutine r_mass_min_part
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
recursive subroutine r_broadcast_mp_min(pst,mp_min,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:output_size)::output_array
  real(kind=8)::mp_min

  integer::rID
  
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BROADCAST_MP_MIN,pst%iUpper+1,input_size,output_size,mp_min)
     call r_broadcast_mp_min(pst%pLower,mp_min,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     pst%s%g%mp_min=mp_min
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
  integer::ilevel
  integer::npart_max

  integer::next_npart_max
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_NPART_MAX,pst%iUpper+1,input_size,output_size,ilevel)
     call r_npart_max(pst%pLower,ilevel,input_size,npart_max,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_npart_max)
     npart_max=MAX(npart_max,next_npart_max)
  else
     npart_max=pst%s%p%npart
  endif

end subroutine r_npart_max
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module input_part_module