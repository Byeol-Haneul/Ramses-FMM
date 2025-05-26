MODULE rt_flux_module

  use rt_parameters, only: nrtvar, nrtgrp
  use amr_parameters, only: ndim, dp
  implicit none

  private   ! default
  integer,parameter::ifrt1=0                                                      ! 0
  integer,parameter::jfrt1=1-ndim/2                                               ! 0 or 1
  integer,parameter::kfrt1=1-ndim/3                                               ! 0 or 1


  real(dp), dimension(100), parameter :: Fx_arr = (/ &
    5.00000000E-09_dp, 6.43452436E-09_dp, 8.28062074E-09_dp, 1.06563712E-08_dp, &
    1.37137360E-08_dp, 1.76482736E-08_dp, 2.27116493E-08_dp, 2.92277321E-08_dp, &
    3.76133108E-08_dp, 4.84047530E-08_dp, 6.22923124E-08_dp, 8.01642803E-08_dp, &
    1.03163803E-07_dp, 1.32762000E-07_dp, 1.70852065E-07_dp, 2.19870355E-07_dp, &
    2.82952231E-07_dp, 3.64132604E-07_dp, 4.68604022E-07_dp, 6.03048798E-07_dp, &
    7.76066436E-07_dp, 9.98723678E-07_dp, 1.28526237E-06_dp, 1.65401040E-06_dp, &
    2.12855404E-06_dp, 2.73924656E-06_dp, 3.52514975E-06_dp, 4.53653238E-06_dp, &
    5.83808562E-06_dp, 7.51306083E-06_dp, 9.66859458E-06_dp, 1.24425615E-05_dp, &
    1.60123930E-05_dp, 2.06064265E-05_dp, 2.65185106E-05_dp, 3.41268005E-05_dp, &
    4.39179458E-05_dp, 5.65182184E-05_dp, 7.27335705E-05_dp, 9.36011860E-05_dp, &
    1.20455822E-04_dp, 1.55015183E-04_dp, 1.99489793E-04_dp, 2.56724383E-04_dp, &
    3.30379852E-04_dp, 4.25167425E-04_dp, 5.47149998E-04_dp, 7.04129929E-04_dp, &
    9.06148089E-04_dp, 1.16612608E-03_dp, 1.50069266E-03_dp, 1.93124726E-03_dp, &
    2.48532847E-03_dp, 3.19837484E-03_dp, 4.11599034E-03_dp, 5.29685859E-03_dp, &
    6.81649038E-03_dp, 8.77204096E-03_dp, 1.12884973E-02_dp, 1.45266148E-02_dp, &
    1.86930773E-02_dp, 2.40534548E-02_dp, 3.09486337E-02_dp, 3.98154359E-02_dp, &
    5.12120423E-02_dp, 6.58483623E-02_dp, 8.46201778E-02_dp, 1.08642800E-01_dp, &
    1.39273292E-01_dp, 1.78096919E-01_dp, 2.26829255E-01_dp, 2.87048524E-01_dp, &
    3.59638450E-01_dp, 4.43856220E-01_dp, 5.36197197E-01_dp, 6.29835311E-01_dp, &
    7.15943180E-01_dp, 7.87239502E-01_dp, 8.41337068E-01_dp, 8.80622192E-01_dp, &
    9.09199036E-01_dp, 9.30424054E-01_dp, 9.46457421E-01_dp, 9.58685105E-01_dp, &
    9.68062145E-01_dp, 9.75278922E-01_dp, 9.80846883E-01_dp, 9.85150356E-01_dp, &
    9.88480827E-01_dp, 9.91060778E-01_dp, 9.93060789E-01_dp, 9.94612076E-01_dp, &
    9.95815819E-01_dp, 9.96750181E-01_dp, 9.97475624E-01_dp, 9.98038968E-01_dp, &
    9.98476499E-01_dp, 9.98816353E-01_dp, 9.99080359E-01_dp, 9.99285459E-01_dp /)

  real(dp), dimension(100), parameter :: b_arr = (/ &
    1.00000000E-08_dp, 1.28690487E-08_dp, 1.65612415E-08_dp, 2.13127423E-08_dp, &
    2.74274719E-08_dp, 3.52965472E-08_dp, 4.54232986E-08_dp, 5.84554642E-08_dp, &
    7.52266217E-08_dp, 9.68095059E-08_dp, 1.24584625E-07_dp, 1.60328561E-07_dp, &
    2.06327606E-07_dp, 2.65524001E-07_dp, 3.41704130E-07_dp, 4.39740709E-07_dp, &
    5.65904461E-07_dp, 7.28265208E-07_dp, 9.37208044E-07_dp, 1.20609760E-06_dp, &
    1.55213287E-06_dp, 1.99744736E-06_dp, 2.57052473E-06_dp, 3.30802080E-06_dp, &
    4.25710808E-06_dp, 5.47849313E-06_dp, 7.05029949E-06_dp, 9.07306476E-06_dp, &
    1.16761712E-05_dp, 1.50261217E-05_dp, 1.93371892E-05_dp, 2.48851229E-05_dp, &
    3.20247859E-05_dp, 4.12128530E-05_dp, 5.30370213E-05_dp, 6.82536011E-05_dp, &
    8.78358917E-05_dp, 1.13036437E-04_dp, 1.45467141E-04_dp, 1.87202373E-04_dp, &
    2.40911645E-04_dp, 3.10030370E-04_dp, 3.98979594E-04_dp, 5.13448783E-04_dp, &
    6.60759740E-04_dp, 8.50334928E-04_dp, 1.09430016E-03_dp, 1.40826021E-03_dp, &
    1.81229692E-03_dp, 2.33225374E-03_dp, 3.00138870E-03_dp, 3.86250173E-03_dp, &
    4.97067230E-03_dp, 6.39678239E-03_dp, 8.23205042E-03_dp, 1.05938658E-02_dp, &
    1.36332975E-02_dp, 1.75447570E-02_dp, 2.25784332E-02_dp, 2.90562957E-02_dp, &
    3.73926884E-02_dp, 4.81208329E-02_dp, 6.19269343E-02_dp, 7.96940734E-02_dp, &
    1.02558691E-01_dp, 1.31983279E-01_dp, 1.69849925E-01_dp, 2.18580696E-01_dp, &
    2.81292563E-01_dp, 3.61996769E-01_dp, 4.65855406E-01_dp, 5.99511591E-01_dp, &
    7.71514387E-01_dp, 9.92865623E-01_dp, 1.27772361E+00_dp, 1.64430873E+00_dp, &
    2.11606892E+00_dp, 2.72317940E+00_dp, 3.50447284E+00_dp, 4.50992317E+00_dp, &
    5.80384209E+00_dp, 7.46899266E+00_dp, 9.61188304E+00_dp, 1.23695791E+01_dp, &
    1.59184716E+01_dp, 2.04855587E+01_dp, 2.63629652E+01_dp, 3.39266284E+01_dp, &
    4.36603433E+01_dp, 5.61867085E+01_dp, 7.23069489E+01_dp, 9.30521648E+01_dp, &
    1.19749284E+02_dp, 1.54105937E+02_dp, 1.98319681E+02_dp, 2.55218564E+02_dp, &
    3.28442013E+02_dp, 4.22673627E+02_dp, 5.43940749E+02_dp, 7.00000000E+02_dp /)


  public rt_unsplit

CONTAINS

function interp_b(Fx) result(b)
  use amr_parameters, only: dp
  implicit none
  real(dp), intent(in) :: Fx
  real(dp) :: b
  integer :: i, n

  n = size(Fx_arr)

  ! Handle extrapolation on the left
  if (Fx <= Fx_arr(1)) then
    b = b_arr(1)
    return
  end if

  ! Handle extrapolation on the right
  if (Fx >= Fx_arr(n)) then
    b = b_arr(n)
    return
  end if

  ! Find i such that Fx_arr(i) <= Fx < Fx_arr(i+1)
  do i = 1, n - 1
    if (Fx < Fx_arr(i + 1)) then
      b = b_arr(i) + (Fx - Fx_arr(i)) * (b_arr(i+1) - b_arr(i)) / (Fx_arr(i+1) - Fx_arr(i))
      return
    end if
  end do

end function interp_b

function bessel_i0(x) result(iv0)
  implicit none
  real(kind=8), intent(in) :: x
  real(kind=8) :: iv0, ax, y
  real(kind=8), parameter :: p1 = 1.0d0,         p2 = 3.5156229d0
  real(kind=8), parameter :: p3 = 3.0899424d0,   p4 = 1.2067492d0
  real(kind=8), parameter :: p5 = 0.2659732d0,   p6 = 0.0360768d0
  real(kind=8), parameter :: p7 = 0.0045813d0
  real(kind=8), parameter :: q1 = 0.39894228d0,  q2 = 0.01328592d0
  real(kind=8), parameter :: q3 = 0.00225319d0,  q4 = -0.00157565d0
  real(kind=8), parameter :: q5 = 0.00916281d0,  q6 = -0.02057706d0
  real(kind=8), parameter :: q7 = 0.02635537d0,  q8 = -0.01647633d0
  real(kind=8), parameter :: q9 = 0.00392377d0

  ax = abs(x)
  if (ax < 3.75d0) then
    y = (x / 3.75d0)**2
    iv0 = p1 + y*(p2 + y*(p3 + y*(p4 + y*(p5 + y*(p6 + y*p7)))))
  else
    y = 3.75d0 / ax
    iv0 = (exp(ax) / sqrt(ax)) * (q1 + y*(q2 + y*(q3 + y*(q4 + y*(q5 + y*(q6 + y*(q7 + y*(q8 + y*q9))))))))
  end if
end function bessel_i0

function bessel_i1(x) result(iv1)
  implicit none
  real(kind=8), intent(in) :: x
  real(kind=8) :: iv1, ax, y

  real(kind=8), parameter :: p1 = 0.5d0,       p2 = 0.87890594d0
  real(kind=8), parameter :: p3 = 0.51498869d0, p4 = 0.15084934d0
  real(kind=8), parameter :: p5 = 0.02658733d0, p6 = 0.00301532d0
  real(kind=8), parameter :: p7 = 0.00032411d0
  real(kind=8), parameter :: q1 = 0.39894228d0, q2 = -0.03988024d0
  real(kind=8), parameter :: q3 = -0.00362018d0, q4 = 0.00163801d0
  real(kind=8), parameter :: q5 = -0.01031555d0, q6 = 0.02282967d0
  real(kind=8), parameter :: q7 = -0.02895312d0, q8 = 0.01787654d0
  real(kind=8), parameter :: q9 = -0.00420059d0

  ax = abs(x)
  if (ax < 3.75d0) then
    y = (x / 3.75d0)**2
    iv1 = x * (p1 + y*(p2 + y*(p3 + y*(p4 + y*(p5 + y*(p6 + y*p7))))))
  else
    y = 3.75d0 / ax
    iv1 = (exp(ax) / sqrt(ax)) * (q1 + y*(q2 + y*(q3 + y*(q4 + y*(q5 + y*(q6 + y*(q7 + y*(q8 + y*q9))))))))
    if (x < 0.0d0) iv1 = -iv1
  end if
end function bessel_i1

function bessel_i2(x) result(iv2)
  implicit none
  real(kind=8), intent(in) :: x
  real(kind=8) :: iv2, x2, iv0, iv1

  if (abs(x) < 1.0d-4) then
    x2 = x*x
    ! Taylor expansion around x=0: I_2(x) ≈ (1/8)x^2 + (1/96)x^4
    iv2 = 0.125d0 * x2 + (1.0d0 / 96.0d0) * x2*x2
  else
    iv0 = bessel_i0(x)
    iv1 = bessel_i1(x)
    iv2 = (2.0d0 / x) * iv1 - iv0
  end if
end function bessel_i2

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
  !---------------------------------------------------
  ! Additional variables added by Harley for 2D computations
  real(dp),dimension(1:ndim,1:ndim)::rotation_matrix
  real(dp),dimension(1:ndim,1:ndim)::pressure_tensor_2D, pressure_tensor_2D_rot
  real(dp),dimension(1:ndim)::F_norm, F_rot
  real(dp)::F_norm_norm, lagrange_a, lagrange_b
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

     if (ndim.eq.2) then
        ! Use Harley pressure tensor for 2D

        F_norm = pflux / (Np * rt_c)
        F_norm_norm = SQRT(DOT_PRODUCT(F_norm, F_norm))

        if (F_norm_norm.lt.1.d-15) then
          ! For small F --> revert to analytic solution
          pressure_tensor_2D(1,1) = 0.5d0
          pressure_tensor_2D(1,2) = 0.0d0
          pressure_tensor_2D(2,1) = 0.0d0
          pressure_tensor_2D(2,2) = 0.5d0
        else
          ! Step 1: Create the rotation matrix
          ! Use unit direction vector to construct rotation matrix
          rotation_matrix(1,1) = F_norm(1)/F_norm_norm           ! cos(θ)
          rotation_matrix(1,2) = F_norm(2)/F_norm_norm           ! sin(θ)
          rotation_matrix(2,1) = -1.d0 * F_norm(2)/F_norm_norm   ! -sin(θ)
          rotation_matrix(2,2) = F_norm(1)/F_norm_norm           ! cos(θ)

          ! Step 2: Apply the rotation matrix to F/|cN|
          F_rot = matmul(rotation_matrix,F_norm)
          
          if (F_rot(1).gt.9.99285459E-01_dp) then
            ! For large F close to 1, use the analytic solution
            pressure_tensor_2D(1,1) = 1.0d0
            pressure_tensor_2D(1,2) = 0.0d0
            pressure_tensor_2D(2,1) = 0.0d0
            pressure_tensor_2D(2,2) = 0.0d0
          else
            ! Step 3: Interpolate the value of lagrange_b 
            ! and then compute lagrange_a
            lagrange_b = interp_b(F_rot(1))
            lagrange_a = LOG(1.d0 / (2.d0 * ACOS(-1.0d0) * bessel_i0(lagrange_b)))

            ! Step 4: Use lagrange a and b to compute the pressure tensor
            pressure_tensor_2D_rot(1,1) = ACOS(-1.0d0) * EXP(lagrange_a) * (bessel_i0(lagrange_b) + bessel_i2(lagrange_b))
            pressure_tensor_2D_rot(2,2) = ACOS(-1.0d0) * EXP(lagrange_a) * (bessel_i0(lagrange_b) - bessel_i2(lagrange_b))
            pressure_tensor_2D_rot(1,2) = 0.d0
            pressure_tensor_2D_rot(2,1) = 0.d0
          endif

          ! Step 5: Rotate the pressure tensor back to the correct frame
          pressure_tensor_2D = matmul(transpose(rotation_matrix), matmul(pressure_tensor_2D_rot, rotation_matrix))
        endif

        ! Step 6: Store the pressure tensor back in F
        ! TODO(code): ASK JOKI WHY ORIENTATION SEEMS OFF AFTER ROTATION
        F(i,j,k,2,1) = pressure_tensor_2D(2,2)
        F(i,j,k,3,2) = pressure_tensor_2D(1,1)
        F(i,j,k,2,2) = -1.d0 * pressure_tensor_2D(1,2)
        F(i,j,k,3,1) = -1.d0 * pressure_tensor_2D(1,2)
     else
        ! Use Levermore 1984 Eddington factor for 3D
        chi = max(4d0-3d0*pflux_sq, 0d0)   !           Eddington factor
        chi = (3d0+ 4d0*pflux_sq)/(5d0 + 2d0*sqrt(chi))
        iterm = (1d0-chi)/2d0                 !    Identity term in tensor
        oterm = (3d0*chi-1d0)/2d0             !         Outer product term
        do p = 1, ndim
            do q = 1, ndim
              F(i,j,k,p+1,q) = oterm * u(p) * u(q)
            enddo
            F(i,j,k,p+1,p) = F(i,j,k,p+1,p) + iterm
        enddo
     endif

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
SUBROUTINE rt_unsplit(uin,flux,cFlx,rt_c,dx,dy,dz,dt &
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
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:ndim+1,1:ndim)::cFlx
  real(dp)::dx, dy, dz, dt, rt_c
  integer::iu1,iu2,ju1,ju2,ku1,ku2
  integer::if1,if2,jf1,jf2,kf1,kf2
  ! Central fluxes:
  ! Upwards and downwards fluxes and states of the group
  real(dp),dimension(nDim+1),save:: fdn, fup, udn, uup
  real(dp)::dtdx
  integer ::i, j, k, iP0, iP1, igrp
!------------------------------------------------------------------------
  do igrp = 1, nrtgrp

  ! Select group
  iP0 = 1 + (igrp-1)*(ndim+1)
  iP1 = iP0+nDim

  ! Compute flux tensors for all the cells with correction
  call cmp_flux_tensors(uin, iP0, cFlx, rt_c, &
    & iu1, iu2, ju1, ju2, ku1, ku2, if2, jf2, kf2)

  ! Solve for 1D flux in X direction
  !----------------------------------------------------------------------
  dtdx = dt/dx
  do i = if1, if2                                 !
  do j = jf1, jf2                                 !        each cell in grid
  do k = kf1, kf2                                 !
     fdn = cFlx(i-1, j, k, :, 1    )    !
     fup = cFlx(i,   j, k, :, 1    )   !  upwards and downwards
     udn = uin( i-1, j, k, iP0:iP1 )   !  conditions
     uup = uin( i,   j, k, iP0:iP1 )    !
     flux(i, j, k, iP0:iP1, 1) = cmp_face( fdn, fup, udn, uup, rt_c )*dtdx
  end do
  end do
  end do

  ! Solve for 1D flux in Y direction
  !----------------------------------------------------------------------
#if NDIM>1
  dtdx = dt/dy
  do i = if1, if2
  do j = jf1, jf2
  do k = kf1, kf2
     fdn = cFlx(i, j-1, k, :, 2    )
     fup = cFlx(i, j,   k, :, 2    )
     udn = uin( i, j-1, k, iP0:iP1 )
     uup = uin( i, j,   k, iP0:iP1 )
     flux(i,j,k,iP0:iP1,2) = cmp_face( fdn, fup, udn, uup, rt_c )*dtdx
  end do
  end do
  end do
#endif

  ! Solve for 1D flux in Z direction
  !----------------------------------------------------------------------
#if NDIM>2
  dtdx = dt/dz
  do i = if1, if2
  do j = jf1, jf2
  do k = kf1, kf2
     fdn = cFlx(i, j, k-1, :, 3    )
     fup = cFlx(i, j, k,   :, 3    )
     udn = uin( i, j, k-1, iP0:iP1 )
     uup = uin( i, j, k,   iP0:iP1 )
     flux(i,j,k,iP0:iP1,3) = cmp_face( fdn, fup, udn, uup, rt_c )*dtdx
  end do
  end do
  end do
#endif

  end do

end subroutine rt_unsplit

END MODULE rt_flux_module
