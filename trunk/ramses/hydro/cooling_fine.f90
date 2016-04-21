subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::igrid,ind,idim
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::d,nH,T2,ekin,etot,eint

  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  do igrid=head(ilevel),tail(ilevel)
     do ind=1,twotondim

        if(.NOT. grid(igrid)%refined(ind))then

           d=max(grid(igrid)%uold(ind,1),smallr)
           etot=grid(igrid)%uold(ind,ndim+2)
           ekin=0.0
           do idim=1,ndim
              ekin=ekin+0.5*grid(igrid)%uold(ind,idim+1)**2/d
           end do
           eint=etot-ekin
           T2=(gamma-1.0)*(eint/d)*scale_T2
           nH=d*scale_nH
           
           ! Set isothermal temperature in Kelvins. 
           T2=T2_star*(1.0+(nH/n_star)**(g_star-1.0))
           
           eint=d*(T2/scale_T2/(gamma-1.0))
           etot=ekin+eint
           grid(igrid)%uold(ind,ndim+2)=etot

        endif

     end do
  end do

111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
