!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_phi(m,igrid_nbor,ind_nbor,ccc,bbb,tfrac,phi_int)
  use amr_parameters, only: ndim,dp,twotondim,threetondim
  use amr_commons, only: mesh_t
  implicit none
  type(mesh_t)::m
  integer,dimension(1:threetondim)::igrid_nbor,ind_nbor
  integer,dimension(1:8,1:8)::ccc
  real(dp),dimension(1:8)::bbb
  real(dp)::tfrac
  real(dp),dimension(1:twotondim)::phi_int
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! - boundary conditions for solving poisson equation at fine level
  ! - computing force (gradient_phi) at fine level for cells close to boundary
  ! Interpolation is performed in space (using CIC) and - if adaptive 
  ! timestepping is on - also in time (using linear extrapolation 
  ! of the change in phi during the last coarse step onto the first fine step)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::ind,ind_average,ind_father
  integer::igrid_nbr,ind_nbr,igrid_cen,ind_cen
  real(dp)::coeff,add
#ifdef GRAV
  ! Store central cell
  igrid_cen=igrid_nbor(threetondim/2+1)
  ind_cen=ind_nbor(threetondim/2+1)

  ! Third order phi interpolation
  do ind=1,twotondim
     phi_int(ind)=0d0
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbb(ind_average)
        igrid_nbr=igrid_nbor(ind_father)
        ind_nbr=ind_nbor(ind_father)
        if (igrid_nbr==0) then 
           write(*,*)'no all neighbors present in interpol_phi...'
           write(*,*)igrid_nbor
           stop
           add=coeff*(m%grid(igrid_cen)%phi(ind_cen)+&
                & (m%grid(igrid_cen)%phi(ind_cen)-m%grid(igrid_cen)%phi_old(ind_cen))*tfrac)
        else
           add=coeff*(m%grid(igrid_nbr)%phi(ind_nbr)+&
                & (m%grid(igrid_nbr)%phi(ind_nbr)-m%grid(igrid_nbr)%phi_old(ind_nbr))*tfrac)
        endif
        phi_int(ind)=phi_int(ind)+add
     end do
  end do
#endif
end subroutine interpol_phi
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_save_phi_old(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel)
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

  integer::next_range,next_cpu

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_SAVE_PHI_OLD,next_cpu,next_range,input_size,output_size,ilevel)
     call r_save_phi_old(r,g,m,p,mdl,next_range,input_size,output_size,ilevel)
  else
     call save_phi_old(m,ilevel)
  endif

end subroutine r_save_phi_old
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine save_phi_old(m,ilevel)
  use amr_parameters, only: ndim,dp,twotondim,threetondim
  use amr_commons, only: mesh_t
  implicit none
  type(mesh_t)::m
  integer ilevel
  ! Save the old potential for time extrapolation in case of subcycling
  integer::ind,igrid

#ifdef GRAV
  ! Loop over level grids
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
        ! Save phi      
        m%grid(igrid)%phi_old(ind)=m%grid(igrid)%phi(ind)
     end do
  end do
#endif

end subroutine save_phi_old
!###########################################################
!###########################################################
!###########################################################
!###########################################################
