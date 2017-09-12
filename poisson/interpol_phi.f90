!###########################################################
!###########################################################
!###########################################################
!###########################################################
#ifdef TOTO
subroutine interpol_phi(igrid_nbor,ind_nbor,ccc,bbb,tfrac,phi_int)
  use amr_commons
  implicit none
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
           add=coeff*(grid(igrid_cen)%phi(ind_cen)+&
                & (grid(igrid_cen)%phi(ind_cen)-grid(igrid_cen)%phi_old(ind_cen))*tfrac)
        else
           add=coeff*(grid(igrid_nbr)%phi(ind_nbr)+&
                & (grid(igrid_nbr)%phi(ind_nbr)-grid(igrid_nbr)%phi_old(ind_nbr))*tfrac)
        endif
        phi_int(ind)=phi_int(ind)+add
     end do
  end do
#endif
 end subroutine interpol_phi
#endif
 !###########################################################
!###########################################################
!###########################################################
!###########################################################
#ifdef TOTO
 subroutine save_phi_old(ilevel)
  use amr_commons
  implicit none
  integer ilevel
  ! Save the old potential for time extrapolation in case of subcycling
  integer::ind,igrid
#ifdef GRAV
  ! Loop over level grids
  do igrid=head(ilevel),tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
        ! Save phi      
        grid(igrid)%phi_old(ind)=grid(igrid)%phi(ind)
     end do
  end do
#endif
end subroutine save_phi_old
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_phi_2(r,g,m,igrid_nbor,ind_nbor,ccc,bbb,tfrac,phi_int)
  use amr_parameters, only: ndim,dp,twotondim,threetondim
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
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
end subroutine interpol_phi_2
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine save_phi_old_2(r,g,m,ilevel)
  use amr_parameters, only: ndim,dp,twotondim,threetondim
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
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
end subroutine save_phi_old_2
