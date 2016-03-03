!################################################################
!################################################################
!################################################################
!################################################################
subroutine synchro_hydro_fine(ilevel,dteff)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  real(dp)::dteff
  !-------------------------------------------------------------------
  ! Update velocity  from gravitational acceleration
  !-------------------------------------------------------------------
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
        ener=grid(igrid)%uold(ind,neul)
        do idim=1,ndim
           ener=ener-0.5*grid(igrid)%uold(ind,idim+1)**2/grid(igrid)%uold(ind,1)
        end do
  
        ! Update momentum
        do idim=1,ndim
           grid(igrid)%uold(ind,idim+1)=grid(igrid)%uold(ind,idim+1)+&
                & grid(igrid)%uold(ind,1)*grid(igrid)%f(ind,idim)*dteff
        end do
  
        ! Update total energy
        do idim=1,ndim
           ener=ener+0.5*grid(igrid)%uold(ind,idim+1)**2/grid(igrid)%uold(ind,1)
        end do
        grid(igrid)%uold(ind,neul)=ener

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

