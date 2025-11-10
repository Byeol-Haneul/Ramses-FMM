module fmm_taylor
contains
#ifdef FMM
subroutine calc_taylor_from_multipole(R, D0, D1, D2, D3, multipoles, taylor_coeff)
  use amr_parameters, only: ndim
  implicit none

  real(kind=8), intent(in)  :: R(ndim), D0, D1, D2, D3
  real(kind=8), intent(in)  :: multipoles(:)
  real(kind=8), intent(out) :: taylor_coeff(:)

  real(kind=8) :: trM, S, RdotM1
  integer :: i, j, k, idxC3
#if NDIM==1
  real(kind=8) :: MR1
#endif
#if NDIM==2
  real(kind=8) :: MR1, MR2
#endif
#if NDIM==3
  real(kind=8) :: MR1, MR2, MR3
#endif

!==========================================================
! NDIM = 1
!==========================================================
#if NDIM==1
  trM = multipoles(3)
  MR1 = multipoles(3)*R(1)
  S = R(1)*MR1
  RdotM1 = R(1)*multipoles(2)

  taylor_coeff(1) = multipoles(1)*D0 - RdotM1*D1 + 0.5d0*(trM*D1 + S*D2)

  taylor_coeff(2) = multipoles(1)*R(1)*D1 - multipoles(2)*D1 &
                    - R(1)*RdotM1*D2 + 0.5d0*R(1)*trM*D2 + 0.5d0*R(1)*S*D3 + MR1*D2

  taylor_coeff(3) = multipoles(1)*( D1 + R(1)**2*D2 ) &
                    - RdotM1*D2 - 2*R(1)*multipoles(2)*D2 - R(1)**2*RdotM1*D3

  taylor_coeff(4) = multipoles(1)*( 3*R(1)*D2 + R(1)**3*D3 )
#endif

!==========================================================
! NDIM = 2
!==========================================================
#if NDIM==2
  trM = multipoles(4) + multipoles(6)

  MR1 = multipoles(4)*R(1) + multipoles(5)*R(2)
  MR2 = multipoles(5)*R(1) + multipoles(6)*R(2)

  S = R(1)*MR1 + R(2)*MR2
  RdotM1 = R(1)*multipoles(2) + R(2)*multipoles(3)

  ! C0
  taylor_coeff(1) = multipoles(1)*D0 - RdotM1*D1 + 0.5d0*(trM*D1 + S*D2)

  ! C1
  taylor_coeff(2) = multipoles(1)*R(1)*D1 - multipoles(2)*D1 + MR1*D2 &
                    - R(1)*RdotM1*D2 + 0.5d0*R(1)*(trM*D2 + S*D3)
  taylor_coeff(3) = multipoles(1)*R(2)*D1 - multipoles(3)*D1 + MR2*D2 &
                    - R(2)*RdotM1*D2 + 0.5d0*R(2)*(trM*D2 + S*D3)

  ! C2
  taylor_coeff(4) = multipoles(1)*( D1 + R(1)**2*D2 ) &
                    - RdotM1*D2 - 2*R(1)*multipoles(2)*D2 - R(1)**2*RdotM1*D3

  taylor_coeff(5) = multipoles(1)*( R(1)*R(2)*D2 ) &
                    - R(1)*multipoles(3)*D2 - R(2)*multipoles(2)*D2 - R(1)*R(2)*RdotM1*D3

  taylor_coeff(6) = multipoles(1)*( D1 + R(2)**2*D2 ) &
                    - RdotM1*D2 - 2*R(2)*multipoles(3)*D2 - R(2)**2*RdotM1*D3

  ! C3
  taylor_coeff(7) = multipoles(1)*( 3*R(1)*D2 + R(1)**3*D3 )
  taylor_coeff(8) = multipoles(1)*( R(2)*D2 + R(1)*R(1)*R(2)*D3 )
  taylor_coeff(9) = multipoles(1)*( R(1)*D2 + R(1)*R(2)*R(2)*D3 )
  taylor_coeff(10)= multipoles(1)*( 3*R(2)*D2 + R(2)**3*D3 )
#endif

!==========================================================
! NDIM = 3
!==========================================================
#if NDIM==3
  trM = multipoles(5) + multipoles(8) + multipoles(10)

  MR1 = multipoles(5)*R(1) + multipoles(6)*R(2) + multipoles(7)*R(3)
  MR2 = multipoles(6)*R(1) + multipoles(8)*R(2) + multipoles(9)*R(3)
  MR3 = multipoles(7)*R(1) + multipoles(9)*R(2) + multipoles(10)*R(3)

  S = R(1)*MR1 + R(2)*MR2 + R(3)*MR3
  RdotM1 = R(1)*multipoles(2) + R(2)*multipoles(3) + R(3)*multipoles(4)

  !------------------------------------------------------------------
  ! C0
  !------------------------------------------------------------------
  taylor_coeff(1) = multipoles(1)*D0 - RdotM1*D1 + 0.5d0*(trM*D1 + S*D2)

  !------------------------------------------------------------------
  ! C1
  !------------------------------------------------------------------
  taylor_coeff(2) = multipoles(1)*R(1)*D1 - multipoles(2)*D1 - R(1)*RdotM1*D2 + 0.5d0*R(1)*trM*D2 + 0.5d0*R(1)*S*D3 + MR1*D2
  taylor_coeff(3) = multipoles(1)*R(2)*D1 - multipoles(3)*D1 - R(2)*RdotM1*D2 + 0.5d0*R(2)*trM*D2 + 0.5d0*R(2)*S*D3 + MR2*D2
  taylor_coeff(4) = multipoles(1)*R(3)*D1 - multipoles(4)*D1 - R(3)*RdotM1*D2 + 0.5d0*R(3)*trM*D2 + 0.5d0*R(3)*S*D3 + MR3*D2

  !------------------------------------------------------------------
  ! C2 (correct order: (1,1),(1,2),(2,2),(1,3),(2,3),(3,3))
  !------------------------------------------------------------------
  taylor_coeff(5) = multipoles(1)*(D1 + R(1)**2*D2) - RdotM1*D2 - 2*R(1)*multipoles(2)*D2           - R(1)**2*RdotM1*D3
  taylor_coeff(6) = multipoles(1)*(R(1)*R(2)*D2)    - R(1)*multipoles(3)*D2 - R(2)*multipoles(2)*D2 - R(1)*R(2)*RdotM1*D3
  taylor_coeff(7) = multipoles(1)*(D1 + R(2)**2*D2) - RdotM1*D2 - 2*R(2)*multipoles(3)*D2           - R(2)**2*RdotM1*D3
  taylor_coeff(8) = multipoles(1)*(R(1)*R(3)*D2)    - R(1)*multipoles(4)*D2 - R(3)*multipoles(2)*D2 - R(1)*R(3)*RdotM1*D3
  taylor_coeff(9) = multipoles(1)*(R(2)*R(3)*D2)    - R(2)*multipoles(4)*D2 - R(3)*multipoles(3)*D2 - R(2)*R(3)*RdotM1*D3
  taylor_coeff(10)= multipoles(1)*(D1 + R(3)**2*D2) - RdotM1*D2 - 2*R(3)*multipoles(4)*D2           - R(3)**2*RdotM1*D3

  !------------------------------------------------------------------
  ! C3 (correct fully hand-unrolled)
  ! Loop order: (1,1,1),(1,1,2),(1,2,2),(2,2,2),(1,1,3),(1,2,3),(2,2,3),(1,3,3),(2,3,3),(3,3,3)
  !------------------------------------------------------------------
  taylor_coeff(11) = multipoles(1)*( 3*R(1)*D2 + R(1)*R(1)*R(1)*D3 )   ! (1,1,1)
  taylor_coeff(12) = multipoles(1)*( R(2)*D2   + R(1)*R(1)*R(2)*D3 )   ! (1,1,2)
  taylor_coeff(13) = multipoles(1)*( R(1)*D2   + R(1)*R(2)*R(2)*D3 )   ! (1,2,2)
  taylor_coeff(14) = multipoles(1)*( 3*R(2)*D2 + R(2)*R(2)*R(2)*D3 )   ! (2,2,2)
  taylor_coeff(15) = multipoles(1)*( R(3)*D2   + R(1)*R(1)*R(3)*D3 )   ! (1,1,3)
  taylor_coeff(16) = multipoles(1)*(             R(1)*R(2)*R(3)*D3 )   ! (1,2,3)
  taylor_coeff(17) = multipoles(1)*( R(3)*D2   + R(2)*R(2)*R(3)*D3 )   ! (2,2,3)
  taylor_coeff(18) = multipoles(1)*( R(1)*D2   + R(1)*R(3)*R(3)*D3 )   ! (1,3,3)
  taylor_coeff(19) = multipoles(1)*( R(2)*D2   + R(2)*R(3)*R(3)*D3 )   ! (2,3,3)
  taylor_coeff(20) = multipoles(1)*( 3*R(3)*D2 + R(3)*R(3)*R(3)*D3 )   ! (3,3,3)
#endif

end subroutine calc_taylor_from_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine calc_phi_from_multipole(R, D0, D1, D2, multipoles, phi_out)
  use amr_parameters, only: ndim
  implicit none

  real(kind=8), intent(in)  :: R(ndim)
  real(kind=8), intent(in)  :: D0,D1,D2
  real(kind=8), intent(in)  :: multipoles(:)
  real(kind=8), intent(out) :: phi_out

#if NDIM==1
  real(kind=8) :: R1, S
  R1 = R(1)
  S = multipoles(3) * R1 * R1
  phi_out = - (multipoles(1)*D0 - multipoles(2)*R1*D1 + 0.5D0*(multipoles(3)*D1 + S*D2))
#endif

#if NDIM==2
  real(kind=8) :: R1,R2, S
  R1 = R(1); R2 = R(2)
  S = multipoles(4)*R1*R1 + 2.0D0*multipoles(5)*R1*R2 + multipoles(6)*R2*R2
  phi_out = - (multipoles(1)*D0 - (multipoles(2)*R1 + multipoles(3)*R2)*D1 + 0.5D0*((multipoles(4)+multipoles(6))*D1 + S*D2))
#endif

#if NDIM==3
  real(kind=8) :: R1,R2,R3, S
  R1 = R(1); R2 = R(2); R3 = R(3)

  ! S = R^T M2 R
  S = multipoles(5)*R1*R1 + 2.0D0*multipoles(6)*R1*R2 + 2.0D0*multipoles(7)*R1*R3 + &
      multipoles(8)*R2*R2 + 2.0D0*multipoles(9)*R2*R3 + multipoles(10)*R3*R3

  phi_out = - (multipoles(1)*D0 - (multipoles(2)*R1 + multipoles(3)*R2 + multipoles(4)*R3)*D1 + &
               0.5D0*( (multipoles(5)+multipoles(8)+multipoles(10))*D1 + S*D2 ))
#endif

end subroutine calc_phi_from_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine shift_taylor(taylor_in, a, taylor_out)
  use amr_parameters, only: ndim, taylor_size
  implicit none

  real(kind=8), intent(in)  :: taylor_in(1:taylor_size)
  real(kind=8), intent(in)  :: a(ndim)
  real(kind=8), intent(out) :: taylor_out(1:taylor_size)

  integer :: i, j, k, idx, nq, no
  integer :: idxC2, idxC3

  ! Dense tensors
  real(kind=8) :: C0, C0p
  real(kind=8) :: C1(ndim), C1p(ndim)
  real(kind=8) :: C2(ndim,ndim), C2p(ndim,ndim)
  real(kind=8) :: C3(ndim,ndim,ndim), C3p(ndim,ndim,ndim)

  ! Number of independent components
  nq = ndim*(ndim+1)/2
  no = ndim*(ndim+1)*(ndim+2)/6

  ! Offsets in flattened array
  idxC2 = 1 + ndim
  idxC3 = idxC2 + nq

  ! ----------------------------
  ! Unpack flattened -> dense
  ! ----------------------------
  C0 = taylor_in(1)

  do i = 1, ndim
     C1(i) = taylor_in(1+i)
  end do

  ! C2 symmetric unpack
  idx = 0
  do j = 1, ndim
     do i = 1, j
        idx = idx + 1
        C2(i,j) = taylor_in(idxC2+idx)
        C2(j,i) = C2(i,j)
     end do
  end do

  ! C3 symmetric unpack
  idx = 0
  do k = 1, ndim
     do j = 1, k
        do i = 1, j
           idx = idx + 1
           C3(i,j,k) = taylor_in(idxC3+idx)
           C3(i,k,j) = C3(i,j,k)
           C3(j,i,k) = C3(i,j,k)
           C3(j,k,i) = C3(i,j,k)
           C3(k,i,j) = C3(i,j,k)
           C3(k,j,i) = C3(i,j,k)
        end do
     end do
  end do

  ! ----------------------------
  ! Apply shift formulas
  ! ----------------------------
  ! C3 stays the same
  C3p = C3

  ! C2' = C2 + a_k C3_{ijk}
  do i = 1, ndim
     do j = 1, ndim
        C2p(i,j) = C2(i,j)
        do k = 1, ndim
           C2p(i,j) = C2p(i,j) + a(k)*C3(i,j,k)
        end do
     end do
  end do

  ! C1' = C1 + a_j C2_{ij} + 0.5 a_j a_k C3_{ijk}
  do i = 1, ndim
     C1p(i) = C1(i)
     do j = 1, ndim
        C1p(i) = C1p(i) + a(j)*C2(i,j)
        do k = 1, ndim
           C1p(i) = C1p(i) + 0.5d0*a(j)*a(k)*C3(i,j,k)
        end do
     end do
  end do

  ! C0' = C0 + a_i C1^i + 0.5 a_i a_j C2^{ij} + 1/6 a_i a_j a_k C3^{ijk}
  C0p = C0
  do i = 1, ndim
     C0p = C0p +  a(i)*C1(i)
     do j = 1, ndim
        C0p = C0p + 0.5 * a(i)*a(j)*C2(i,j)
        do k = 1, ndim
           C0p = C0p + (1.0/6.0) * a(i)*a(j)*a(k)*C3(i,j,k)
        end do
     end do
  end do

  ! ----------------------------
  ! Pack dense -> flattened
  ! ----------------------------
  taylor_out(1) = C0p

  do i = 1, ndim
     taylor_out(1+i) = C1p(i)
  end do

  idx = 0
  do j = 1, ndim
     do i = 1, j
        idx = idx + 1
        taylor_out(idxC2+idx) = C2p(i,j)
     end do
  end do

  idx = 0
  do k = 1, ndim
     do j = 1, k
        do i = 1, j
           idx = idx + 1
           taylor_out(idxC3+idx) = C3p(i,j,k)
        end do
     end do
  end do
end subroutine shift_taylor
!################################################################
!################################################################
!################################################################
!################################################################
subroutine calc_phi(taylor_in, a, phi_out)
  use amr_parameters, only: ndim, taylor_size
  implicit none

  real(kind=8), intent(in)  :: taylor_in(1:taylor_size)
  real(kind=8), intent(in)  :: a(ndim)
  real(kind=8), intent(out) :: phi_out
#if NDIM==1
  associate( t => taylor_in )
     phi_out = - t(1) - a(1)*t(2) - 0.5d0*a(1)**2*t(3) - (1.0d0/6.0d0)*a(1)**3*t(4)
  end associate
#elif NDIM==2
  associate( t => taylor_in )
     ! t(1) = C0, t(2:3) = C1, t(4:6) = C2, t(7:10) = C3
     phi_out = - t(1) &
               - a(1)*t(2) - a(2)*t(3) &
               - 0.5d0*( a(1)**2*t(4) + 2.0d0*a(1)*a(2)*t(5) + a(2)**2*t(6) ) &
               - (1.0d0/6.0d0)*( a(1)**3*t(7) + 3.0d0*a(1)**2*a(2)*t(8) + 3.0d0*a(1)*a(2)**2*t(9) + a(2)**3*t(10) )
  end associate
#elif NDIM==3
  associate( t => taylor_in )
     ! t(1) = C0, t(2:4) = C1, t(5:10) = C2, t(11:20) = C3
     phi_out = - t(1) &
               - a(1)*t(2) - a(2)*t(3) - a(3)*t(4) &
               - 0.5d0*( a(1)**2*t(5) + 2.0d0*a(1)*a(2)*t(6) + a(2)**2*t(7) &
               + 2.0d0*a(1)*a(3)*t(8) + 2.0d0*a(2)*a(3)*t(9) + a(3)**2*t(10) ) & 
               - (1.0d0/6.0d0)*( a(1)**3*t(11) + 3.0d0*a(1)**2*a(2)*t(12) + 3.0d0*a(1)**2*a(3)*t(15) &
                                + 3.0d0*a(1)*a(2)**2*t(13) + 6.0d0*a(1)*a(2)*a(3)*t(16) + 3.0d0*a(1)*a(3)**2*t(18) &
                                + a(2)**3*t(14) + 3.0d0*a(2)**2*a(3)*t(17) + 3.0d0*a(2)*a(3)**2*t(19) &
                                + a(3)**3*t(20) )
  end associate
#endif
end subroutine calc_phi
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_displacement(p, q, boxlen, r)
  use amr_parameters, only: ndim
  implicit none
  real(kind=8), intent(in)  :: p(ndim), q(ndim)
  real(kind=8), intent(in)  :: boxlen
  real(kind=8), intent(out) :: r(ndim)

  ! Compute raw difference
  r = p - q

  ! Apply periodic boundary conditions (vectorized)
  !where (r >  boxlen / 2.d0) r = r - boxlen
  !where (r < -boxlen / 2.d0) r = r + boxlen
end subroutine get_displacement
#endif
end module fmm_taylor
