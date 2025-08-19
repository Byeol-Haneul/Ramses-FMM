module move_fine_module
  use rho_fine_module, only: cic_weight, cic_index, tsc_weight, tsc_index, pcs_weight, pcs_index
contains
!################################################################
!################################################################
!################################################################
!################################################################
subroutine m_kick_drift_part(pst,ilevel,action_part)
  use amr_parameters, only: ndim, twotondim
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  integer::ilevel
  integer::action_part
  !--------------------------------------------------------------
  ! Move particles according to kick-drift leap frog scheme.
  !--------------------------------------------------------------
  integer,dimension(1:2)::input_array
  integer::dummy(2)

  if(pst%s%m%noct_tot(ilevel)==0)return
  if(pst%s%r%verbose)write(*,'("   Entering kick_drift_part for level",i2," and action=",i2)')ilevel,action_part

  input_array(1)=ilevel
  input_array(2)=action_part
  call r_kick_drift_part(pst,input_array,2,dummy,0)

end subroutine m_kick_drift_part
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_kick_drift_part(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel
  integer::action_part
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_KICK_DRIFT_PART,pst%iUpper+1,input_size,output_size,input_array)
     call r_kick_drift_part(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     ilevel=input_array(1)
     action_part=input_array(2)
     ! Force interpolation for various components (DM particles, star, sink, tree, trac)
     ! based on their respective deposition schemes (CIC 1, TSC 2 or PCS 3)
     if(pst%s%r%part)then
        if(pst%s%r%part_force_interpolation_scheme==1)then
           call cic_kick_drift_part(pst%s,pst%s%p   ,ilevel,action_part)
        elseif(pst%s%r%part_force_interpolation_scheme==2)then
           call tsc_kick_drift_part(pst%s,pst%s%p   ,ilevel,action_part)
        elseif(pst%s%r%part_force_interpolation_scheme==3)then
           call pcs_kick_drift_part(pst%s,pst%s%p   ,ilevel,action_part)
        endif
     endif
     if(pst%s%r%star)then
        if(pst%s%r%star_force_interpolation_scheme==1)then
           call cic_kick_drift_part(pst%s,pst%s%star,ilevel,action_part)
        elseif(pst%s%r%star_force_interpolation_scheme==2)then
           call tsc_kick_drift_part(pst%s,pst%s%star,ilevel,action_part)
        elseif(pst%s%r%star_force_interpolation_scheme==3)then
           call pcs_kick_drift_part(pst%s,pst%s%star,ilevel,action_part)
        endif
     endif
     if(pst%s%r%sink)then
        if(pst%s%r%sink_force_interpolation_scheme==1)then
           call cic_kick_drift_part(pst%s,pst%s%sink,ilevel,action_part)
        elseif(pst%s%r%sink_force_interpolation_scheme==2)then
           call tsc_kick_drift_part(pst%s,pst%s%sink,ilevel,action_part)
        elseif(pst%s%r%sink_force_interpolation_scheme==3)then
           call pcs_kick_drift_part(pst%s,pst%s%sink,ilevel,action_part)
        endif
     endif
     if(pst%s%r%tree)then
        if(pst%s%r%tree_force_interpolation_scheme==1)then
           call cic_kick_drift_part(pst%s,pst%s%tree,ilevel,action_part)
        elseif(pst%s%r%tree_force_interpolation_scheme==2)then
           call tsc_kick_drift_part(pst%s,pst%s%tree,ilevel,action_part)
        elseif(pst%s%r%tree_force_interpolation_scheme==3)then
           call pcs_kick_drift_part(pst%s,pst%s%tree,ilevel,action_part)
        endif
     endif
     if(pst%s%r%trac)then
        if(pst%s%r%trac_force_interpolation_scheme==1)then
           call cic_kick_drift_part(pst%s,pst%s%trac,ilevel,action_part)
        elseif(pst%s%r%trac_force_interpolation_scheme==2)then
           call tsc_kick_drift_part(pst%s,pst%s%trac,ilevel,action_part)
        elseif(pst%s%r%trac_force_interpolation_scheme==3)then
           call pcs_kick_drift_part(pst%s,pst%s%trac,ilevel,action_part)
        endif
     endif
  endif

end subroutine r_kick_drift_part
!################################################################
!################################################################
!################################################################
!################################################################
! Reuse underlying interpolation kernels unchanged
include "../../../pm/move_fine.f90"
!################################################################
!################################################################
!################################################################
!################################################################
end module move_fine_module


