!################################################################
!################################################################
!################################################################
!################################################################
subroutine synchro_hydro_fine(ilevel,dteff)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  real(dp)::dteff
  !--------------------------------------------------------------
  ! Add gravity source terms to uold with time step dteff.
  !--------------------------------------------------------------
  integer::igrid,ind
  integer::idim,neul=ndim+2
  real(dp)::ener

  if(.not. poisson)return
  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel,dteff

  ! Loop over octs
  do igrid=head(ilevel),tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim

        ! Remove kinetic energy from total energy
        ener=fluid(igrid)%uold(ind,neul)
        do idim=1,ndim
           ener=ener-0.5*fluid(igrid)%uold(ind,idim+1)**2/max(fluid(igrid)%uold(ind,1),smallr)
        end do
  
        ! Update momentum
        do idim=1,ndim
           fluid(igrid)%uold(ind,idim+1)=fluid(igrid)%uold(ind,idim+1)+&
                & max(fluid(igrid)%uold(ind,1),smallr)*grav(igrid)%f(ind,idim)*dteff
        end do

        ! Update total energy
        do idim=1,ndim
           ener=ener+0.5*fluid(igrid)%uold(ind,idim+1)**2/max(fluid(igrid)%uold(ind,1),smallr)
        end do
        fluid(igrid)%uold(ind,neul)=ener

     end do
     ! End loop over cells
  end do
  ! End loop over grids

111 format('   Entering synchro_hydro_fine for level',i2,' and time step dt=',1PE12.5)

end subroutine synchro_hydro_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_gravity_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------
  ! This routine adds to unew the gravity source terms to unew
  ! with only half a time step. Only the momentum and the
  ! total energy are modified in array unew.
  !--------------------------------------------------------------
  integer::igrid,ivar,ind
  real(dp)::d,u,v,w,e_kin,e_prim,d_old,fact

  if(.not. poisson)return
  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Add gravity source term at time t with half time step
  do igrid=head(ilevel),tail(ilevel)
     do ind=1,twotondim

        d=max(fluid(igrid)%unew(ind,1),smallr)
        u=0.0; v=0.0; w=0.0
        if(ndim>0)u=fluid(igrid)%unew(ind,2)/d
        if(ndim>1)v=fluid(igrid)%unew(ind,3)/d
        if(ndim>2)w=fluid(igrid)%unew(ind,4)/d
        e_kin=0.5*d*(u**2+v**2+w**2)
        e_prim=fluid(igrid)%unew(ind,ndim+2)-e_kin
        d_old=max(fluid(igrid)%uold(ind,1),smallr)
        fact=d_old/d*0.5*dtnew(ilevel)
        if(ndim>0)then
           u=u+grav(igrid)%f(ind,1)*fact
           fluid(igrid)%unew(ind,2)=d*u
        endif
        if(ndim>1)then
           v=v+grav(igrid)%f(ind,2)*fact
           fluid(igrid)%unew(ind,3)=d*v
        end if
        if(ndim>2)then
           w=w+grav(igrid)%f(ind,3)*fact
           fluid(igrid)%unew(ind,4)=d*w
        endif
        e_kin=0.5*d*(u**2+v**2+w**2)
        fluid(igrid)%unew(ind,ndim+2)=e_prim+e_kin
     end do
  end do

111 format('   Entering add_gravity_source_terms for level ',i2)

end subroutine add_gravity_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################

