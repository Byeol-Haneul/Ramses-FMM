!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine poisson_flag(s,ilevel)
  use amr_parameters, only: dp,ndim,twotondim,twopi
  use ramses_commons, only: ramses_t
  use hydro_parameters, only: nvar
  use hash
  implicit none
  type(ramses_t)::s
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine flag for refinement cells that satisfies
  ! some user-defined physical criteria at the level ilevel. 
  ! -------------------------------------------------------------------
  real(dp)::dx_loc,vol_loc,d_scale,factG
  real(dp),dimension(1:nvar)::uu
  integer::igrid,ind,ivar
  logical::ok

  associate(r=>s%r,g=>s%g,m=>s%m)

  if(        r%m_refine(ilevel).LE.-1.0 .and.&
       & r%jeans_refine(ilevel).LE.-1.0 )return

  ! Constants
  factG=1
  if(r%cosmo)factG=3d0/4d0/twopi*g%omega_m*g%aexp
  dx_loc=r%boxlen/2**ilevel
  vol_loc=dx_loc**3
  d_scale=r%mass_sph/vol_loc

  ! Loop over active grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Loop over cells
     do ind=1,twotondim

        ! Initialize refinement to false
        ok=.false.

        ! Flag cells with density beyond the threshold
#ifndef GRAV
#ifdef HYDRO
        if(r%mass_sph>0.and.r%m_refine(ilevel)>=0)then
           ok=(ok .or. m%grid(igrid)%uold(ind,1)>=r%m_refine(ilevel)*d_scale)
        endif
#endif
#endif
        ! Flag cells with density beyond the threshold
#ifdef GRAV
        if(r%m_refine(ilevel)>=0)then
           ok=(ok .or. m%grid(igrid)%nref(ind)>=r%m_refine(ilevel))
        endif
#endif

        if(r%jeans_refine(ilevel)>=0)then
           ! Gather hydro variables
#ifdef HYDRO
           do ivar=1,nvar
              uu(ivar)=m%grid(igrid)%uold(ind,ivar)
           end do
           call jeans_length_refine(r,uu,factG,dx_loc,r%jeans_refine(ilevel),ok)
#endif
        endif
        
        ! Count only newly flagged cells
        if(m%grid(igrid)%flag1(ind)==0.and.ok)g%nflag=g%nflag+1
        if(ok)m%grid(igrid)%flag1(ind)=1

     end do
     ! End loop over cells
  end do
  ! End loop over grids

  end associate

end subroutine poisson_flag
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine jeans_length_refine(r,uu,factG,size_cell,n_jeans,ok)
  use amr_parameters, only: ndim,dp,twopi
  use amr_commons, only: run_t
  use hydro_parameters, only: nvar
  use const
  implicit none
  type(run_t)::r
  real(dp)::uu(1:nvar)
  real(dp)::n_jeans,factG,size_cell
  logical ::ok
  ! 
  real(dp)::lamb_jeans
  real(dp)::dens,tempe,etherm

  ! compute the thermal energy
  dens = max( uu(1) , r%smallr )
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
  tempe =  etherm / dens * (r%gamma - one)
  ! prevent numerical crash due to negative temperature
  tempe = max( tempe , r%smallc**2 )

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

