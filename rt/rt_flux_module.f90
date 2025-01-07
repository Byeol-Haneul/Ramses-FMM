MODULE rt_flux_module

  use rt_parameters, only: nrtvar
  use amr_parameters, only: ndim, dp
  implicit none

  private   ! default
  integer,parameter::ifrt1=0                                                      ! 0
  integer,parameter::jfrt1=1-ndim/2                                               ! 0 or 1
  integer,parameter::kfrt1=1-ndim/3                                               ! 0 or 1

  public rt_unsplit

CONTAINS

!************************************************************************
subroutine cmp_flux_tensors(uin, iP0, F, rt_c &
     &                     ,iu1,iu2,ju1,ju2,ku1,ku2,if2,jf2,kf2)
  !------------------------------------------------------------------------  
  ! Compute central fluxes for a photon group, for each cell in a vector 
  ! of grids. 
  ! The flux tensor is a three by four tensor (2*3 and 1*2 in 1D and 2D, 
  ! respectively) where the first column is photon flux (x,y,z) and 
  ! the other three columns compose the Eddington tensor (see Aubert & 
  ! Teyssier '08), times c^2. 
  ! input/output:
  ! uin       => RT variables of all cells in a vector of grids
  !              (photon energy densities and photon fluxes).
  ! iP0       => Starting index of photon group among the RT variables.
  ! F        <=  Group flux tensors for all the cells.
  !------------------------------------------------------------------------
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nrtvar)::uin
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nDim+1,1:ndim)::F
  real(dp)::rt_c
  integer::iu1,iu2,ju1,ju2,ku1,ku2
  integer::if2,jf2,kf2
  integer::iP0
  !---------------------------------------------------
  real(dp),dimension(1:ndim)::pflux, u
  real(dp)::Np, Np_c_sq, pflux_sq, chi, iterm, oterm
  integer::i, j, k, p, q, nedge
  !------------------------------------------------------------------------
  ! Loop (N+2)X(N+2)X(N+2) cells in grid, where N=2**(nsuperoct+1) = 2 by 
  ! default. All dimension indices go from 0 to N+1.
  ! We only need to calculate tensors for those cells which have faces to
  ! the NXNXN center cells, so by skipping the 'corners' we are reduced
  ! to fewer cells to calculate (by half for the default N=2).
  do k = kfrt1, kf2
  do j = jfrt1, jf2
  do i = ifrt1, if2

     nedge=0                 ! Check if we're at a corner and if so, cycle
     if(mod(i,if2).eq.0) nedge=nedge+1
     if(ndim.gt.1 .and. mod(j,jf2).eq.0) nedge=nedge+1
     if(ndim.gt.2 .and. mod(k,kf2).eq.0) nedge=nedge+1
     if(nedge.ge.2) cycle

     Np =   uin(i, j, k, iP0)      !          Photon density in cell
     pflux= uin(i, j, k, iP0+1:iP0+ndim)  !       Photon flux vector
     if(Np .lt. 0d0) then
        write(*,*)'negative photon density in cmp_eddington. -EXITING-'
        stop
     endif
     F(i,j,k,1,1:nDim)= pflux            !   First col is photon flux
     ! Rest is Eddington tensor...initalize it to zero
     F(i,j,k,2:ndim+1,1:nDim) = 0d0   
     Np_c_sq = Np * rt_c**2 * Np
     if(Np_c_sq .eq. 0d0) cycle            !Zero density => no pressure
        
     pflux_sq = sum(pflux**2)              !  Sq. photon flux magnitude
     u(:) = 0d0                            !           Flux unit vector
     if(pflux_sq .gt. 0d0) u(:) = pflux/sqrt(pflux_sq)
     pflux_sq = pflux_sq/Np_c_sq           !      Reduced flux, squared
     chi = max(4d0-3d0*pflux_sq, 0d0)      !           Eddington factor
     chi = (3d0+ 4d0*pflux_sq)/(5d0 + 2d0*sqrt(chi))

     iterm = (1d0-chi)/2d0                 !    Identity term in tensor
     oterm = (3d0*chi-1d0)/2d0             !         Outer product term
     do p = 1, ndim
        do q = 1, ndim
           F(i,j,k,p+1,q) = oterm * u(p) * u(q)
        enddo
        F(i,j,k,p+1,p) = F(i,j,k,p+1,p) + iterm
     enddo
     F(i, j, k, 2:ndim+1, 1:ndim) =                                &
          F(i, j, k, 2:ndim+1, 1:ndim) * rt_c**2 * Np
  enddo
  enddo
  enddo

end subroutine cmp_flux_tensors

!************************************************************************
FUNCTION cmp_face(fdn, fup, udn, uup, rt_c)
  
! Compute intercell fluxes for all (four) RT variables, using the
! Harten-Lax-van Leer method (see eq. 30 in Aubert&Teyssier(2008).
! fdn    => flux function in the cell downwards from the border
! fup    => flux function in the cell upwards from the border
! udn    => state (photon density and flux downwards from the border
! uup    => state (photon density and flux upwards from the border
! returns      flux vector for the given state variables, i.e. line nr dim
!              in the 3*4 flux function tensor
!------------------------------------------------------------------------
  real(dp),dimension(nDim+1)::fdn, fup, udn, uup, cmp_face
  real(dp)::rt_c
!------------------------------------------------------------------------
  cmp_face = ( fdn + fup - rt_c*( uup-udn )) / 2d0
  return
END FUNCTION cmp_face


!************************************************************************
SUBROUTINE rt_unsplit(uin,flux,rt_c,dx,dy,dz,dt &
                        ,iu1,iu2,ju1,ju2,ku1,ku2,if1,if2,jf1,jf2,kf1,kf2)

!  Compute intercell fluxes for one photon group in all dimensions,
!  using the Eddington tensor with the M1 closure relation.
!  The intercell fluxes are the right-hand sides of the equations:
!      dN/dt = - nabla(F),
!      dF/dt = - nabla(c^2*P),
!  where N is photon density, F is photon flux, c the light speed and P
!  the Eddington pressure tensor. A flux at index i,j,k represents
!  flux across the lower faces of that cell, i.e. at i-1/2 etc.
!
!  inputs/outputs
!  uin         => input states
!  flux       <=  return fluxes in the 3 coord directions.
!  dx,dy,dz    => (dx,dy,dz)
!  dt          => time step
!
!  other vars
!  iu1,iu2     |First and last index of input array,
!  ju1,ju2     |cell centered,
!  ku1,ku2     |including buffer cells (6x6x6).
!  if1,if2     |First and last index of output array,
!  jf1,jf2     |edge centered, for active
!  kf1,kf2     |cells only (3x3x3).
!------------------------------------------------------------------------
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nrtvar)::       uin
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2,1:nrtvar,1:ndim)::flux
  real(dp)::dx, dy, dz, dt, rt_c
  integer::iu1,iu2,ju1,ju2,ku1,ku2
  integer::if1,if2,jf1,jf2,kf1,kf2
  ! Central fluxes:
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2, ndim+1, ndim)::   cFlx
  ! Upwards and downwards fluxes and states of the group
  real(dp),dimension(nDim+1),save:: fdn, fup, udn, uup
  real(dp)::dtdx
  integer ::i, j, k, iP0, iP1
!------------------------------------------------------------------------
  iP0=1                                    ! For now just using one group
  iP1=iP0+nDim                             ! end index of photon group

  ! compute flux tensors for all the cells with correction
  call cmp_flux_tensors(uin, iP0, cFlx, rt_c &
                        ,iu1,iu2,ju1,ju2,ku1,ku2,if2,jf2,kf2)

  ! Solve for 1D flux in X direction
  !----------------------------------------------------------------------
  dtdx=dt/dx
  do i=if1,if2                                 !
  do j=jf1,jf2                                 !        each cell in grid
  do k=kf1,kf2                                 !
     fdn = cFlx(i-1, j, k, :, 1    )    !
     fup = cFlx(i,   j, k, :, 1    )   !  upwards and downwards
     udn = uin( i-1, j, k, iP0:iP1 )   !  conditions
     uup = uin( i,   j, k, iP0:iP1 )    !
     flux(i, j, k, iP0:iP1, 1)=&
          cmp_face( fdn, fup, udn, uup, rt_c )*dtdx
  end do
  end do
  end do

  ! Solve for 1D flux in Y direction
  !----------------------------------------------------------------------
#if NDIM>1
  dtdx=dt/dy
  do i=if1,if2
  do j=jf1,jf2
  do k=kf1,kf2
     fdn = cFlx(i, j-1, k, :, 2    )
     fup = cFlx(i, j,   k, :, 2    )
     udn = uin( i, j-1, k, iP0:iP1 )
     uup = uin( i, j,   k, iP0:iP1 )
     flux(i,j,k,iP0:iP1,2)=&
          cmp_face( fdn, fup, udn, uup, rt_c )*dtdx
  end do
  end do
  end do
#endif

  ! Solve for 1D flux in Z direction
  !----------------------------------------------------------------------
#if NDIM>2
  dtdx=dt/dz
  do i=if1,if2
  do j=jf1,jf2
  do k=kf1,kf2
     fdn = cFlx(i, j, k-1, :, 3    )
     fup = cFlx(i, j, k,   :, 3    )
     udn = uin( i, j, k-1, iP0:iP1 )
     uup = uin( i, j, k,   iP0:iP1 )
     flux(i,j,k,iP0:iP1,3)=&
          cmp_face( fdn, fup, udn, uup, rt_c )*dtdx
  end do
  end do
  end do
#endif

end subroutine rt_unsplit

END MODULE rt_flux_module
