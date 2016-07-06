subroutine poisson_flag(ilevel)
  use amr_commons
  use hydro_commons
  use pm_commons, ONLY: mp_min
  use hash
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine flag for refinement cells that satisfies
  ! some user-defined physical criteria at the level ilevel. 
  ! -------------------------------------------------------------------
  real(dp)::dx_loc,vol_loc,d_scale,dp_scale,factG
  real(dp),dimension(1:nvar),save::uu
  integer::igrid,ind,ivar
  logical::ok

  if(ilevel==nlevelmax)return
  if(noct_tot(ilevel)==0)return

  if(        m_refine(ilevel).LE.-1.0 .and.&
       & jeans_refine(ilevel).LE.-1.0 )return

  ! Constants
  factG=1
  if(cosmo)factG=3d0/4d0/twopi*omega_m*aexp
  dx_loc=boxlen/2**ilevel
  vol_loc=dx_loc**3
  d_scale=mass_sph/vol_loc
  dp_scale=mp_min/vol_loc

  ! Loop over active grids
  do igrid=head(ilevel),tail(ilevel)

     ! Loop over cells
     do ind=1,twotondim

        ! Initialize refinement to false
        ok=.false.

        ! Flag cells with density beyond the threshold
#ifdef HYDRO
        if(mass_sph>0)then
           if(m_refine(ilevel)>=0)then
              ok=(ok .or. grid(igrid)%uold(ind,1)>=m_refine(ilevel)*d_scale)
           endif
        endif
#endif
        
        ! Flag cells with density beyond the threshold
#ifdef GRAV
        if(pic.and.mp_min>0)then
           if(m_refine(ilevel)>=0)then
              ok=(ok .or. grid(igrid)%rho(ind)>=m_refine(ilevel)*dp_scale)
           endif
        endif
#endif

        if(jeans_refine(ilevel)>=0)then
           ! Gather hydro variables
#ifdef HYDRO
           do ivar=1,nvar
              uu(ivar)=grid(igrid)%uold(ind,ivar)
           end do
           call jeans_length_refine(uu,factG,dx_loc,jeans_refine(ilevel),ok)
#endif
        endif
        
        ! Count only newly flagged cells
        if(grid(igrid)%flag1(ind)==0.and.ok)nflag=nflag+1
        if(ok)grid(igrid)%flag1(ind)=1

     end do
     ! End loop over cells
  end do
  ! End loop over grids

end subroutine poisson_flag
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine jeans_length_refine(uu,factG,size_cell,n_jeans,ok)
  use amr_parameters
  use hydro_commons
  use const
  implicit none
  real(dp)::uu(1:nvar)
  real(dp)::n_jeans,factG,size_cell
  logical ::ok
  ! 
  real(dp)::lamb_jeans
  real(dp)::dens,tempe,emag,etherm

  ! compute the thermal energy
  dens = max( uu(1) , smallr )
  etherm = uu(ndim+2)
  etherm = etherm - 0.5d0*uu(2)**2/dens
#if NDIM > 1
  etherm = etherm - 0.5d0*uu(3)**2/dens
#endif
#if NDIM > 2
  etherm = etherm - 0.5d0*uu(4)**2/dens
#endif
#if NENER>0
  do irad=1,nener
     etherm=etherm-uu(ndim+2+irad)
  end do
#endif

  ! compute the temperature
  tempe =  etherm / dens * (gamma - one)
  ! prevent numerical crash due to negative temperature
  tempe = max( tempe , smallc**2 )

  ! compute the Jeans length
  lamb_jeans = sqrt( tempe * twopi / two / dens / factG )

  ! the Jeans length must be smaller
  ! than n_jeans times the size of the pixel
  ok = ok .or. ( n_jeans*size_cell >= lamb_jeans )

end subroutine jeans_length_refine
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################

