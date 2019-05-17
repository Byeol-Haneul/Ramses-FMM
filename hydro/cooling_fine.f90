!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_cooling_fine(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_COOLING_FINE,pst%iUpper+1,input_size,output_size,input_array)
     call r_cooling_fine(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     ilevel=input_array(1)
     call cooling_fine(pst%s%r,pst%s%g,pst%s%m,ilevel)
  endif

end subroutine r_cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cooling_fine(r,g,m,ilevel)
  use amr_parameters, only: ndim,twotondim,dp
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::igrid,ind,idim
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::d,nH,T2,ekin,etot,eint

#ifdef HYDRO

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim

        if(.NOT. m%grid(igrid)%refined(ind))then

           d=max(m%grid(igrid)%uold(ind,1),r%smallr)
           etot=m%grid(igrid)%uold(ind,ndim+2)
           ekin=0.0
           do idim=1,ndim
              ekin=ekin+0.5*m%grid(igrid)%uold(ind,idim+1)**2/d
           end do
           eint=etot-ekin
           T2=(r%gamma-1.0)*(eint/d)*scale_T2
           nH=d*scale_nH
           
           ! Set isothermal temperature in Kelvins. 
           T2=r%T2_star*(1.0+(nH/r%n_star)**(r%g_star-1.0))
           
           eint=d*(T2/scale_T2/(r%gamma-1.0))
           etot=ekin+eint
           m%grid(igrid)%uold(ind,ndim+2)=etot

        endif

     end do
  end do

#endif

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
