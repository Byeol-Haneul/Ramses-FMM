module fmm_taylor
contains
#ifdef FMM
subroutine calc_taylor_from_multipole(R, D0, D1, D2, D3, multipoles, taylor_coeff)
  use amr_parameters, only: ndim, multipole_size, taylor_size
  implicit none

  integer :: nq, no, i, j, k, idxC2, idxC3
  real(kind=8), intent(in)  :: R(ndim)
  real(kind=8), intent(in)  :: D0,D1,D2,D3
  real(kind=8), intent(in)  :: multipoles(1:multipole_size)
  real(kind=8), intent(out) :: taylor_coeff(1:taylor_size)

  ! Multipoles
  real(kind=8) :: M0
  real(kind=8) :: M1(ndim)
  real(kind=8) :: M2(ndim,ndim)

  ! Displacement
  real(kind=8) :: dist

  ! Precompute
  real(kind=8) :: trM, S, MR(ndim), RdotM1

  !------------------------------------------------------------------
  ! Reconstruct multipoles from 1D array
  !------------------------------------------------------------------
  M0 = multipoles(1)
  M1 = multipoles(2:1+ndim)
  M2 = 0.0D0
#if NDIM==1
  M2(1,1) = multipoles(3)
#endif
#if NDIM==2
  M2(1,1) = multipoles(4)
  M2(1,2) = multipoles(5); M2(2,1) = M2(1,2)
  M2(2,2) = multipoles(6)
#endif
#if NDIM==3
  M2(1,1) = multipoles(5)
  M2(1,2) = multipoles(6); M2(2,1) = M2(1,2)
  M2(1,3) = multipoles(7); M2(3,1) = M2(1,3)
  M2(2,2) = multipoles(8)
  M2(2,3) = multipoles(9); M2(3,2) = M2(2,3)
  M2(3,3) = multipoles(10)
#endif

  !------------------------------------------------------------------
  ! Precompute quantities
  !------------------------------------------------------------------
  trM = 0.0D0
  do i = 1, ndim
     trM = trM + M2(i,i)
  end do

  MR = matmul(M2, R)
  S = dot_product(R, MR)
  RdotM1 = sum(R(:) * M1(:))

  !------------------------------------------------------------------
  ! Compute Taylor coefficients
  !------------------------------------------------------------------
  !------------------------------------------------------------------
  !                               C0
  !------------------------------------------------------------------  
  taylor_coeff(1) = M0*D0 - RdotM1*D1 + 0.5D0*(trM*D1 + S*D2)

  !------------------------------------------------------------------
  !                               C1
  !------------------------------------------------------------------  
  do i = 1, ndim
     taylor_coeff(1+i) = M0*R(i)*D1 - M1(i)*D1 - R(i)*RdotM1*D2 + 0.5D0*R(i)*trM*D2 + 0.5D0*R(i)*S*D3 + MR(i)*D2 
  end do

  !------------------------------------------------------------------
  !                               C2
  !------------------------------------------------------------------
  nq = ndim*(ndim+1)/2
  idxC2 = 1 + ndim
  do j = 1, ndim
     do i = 1, j
        idxC2 = idxC2 + 1
        taylor_coeff(idxC2) = M0*( merge(D1, 0.0D0, i==j) + R(i)*R(j)*D2 ) &
                              - merge(RdotM1, 0.0D0, i==j)*D2 - R(i)*M1(j)*D2 - R(j)*M1(i)*D2 - R(i)*R(j)*RdotM1*D3
     end do
  end do

  !------------------------------------------------------------------
  !                               C3
  !------------------------------------------------------------------
  no = ndim*(ndim+1)*(ndim+2)/6
  idxC3 = 1 + ndim + nq
  do k = 1, ndim
     do j = 1, k
        do i = 1, j
           idxC3 = idxC3 + 1
           taylor_coeff(idxC3) = M0*( (merge(R(k),0.0D0,i==j) + merge(R(i),0.0D0,j==k) + merge(R(j),0.0D0,k==i))*D2 + R(i)*R(j)*R(k)*D3 )
        end do
     end do
  end do

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

  integer :: i, j, k, idx, nq, no
  integer :: idxC2, idxC3

  ! Dense tensors
  real(kind=8) :: C0
  real(kind=8) :: C1(ndim)
  real(kind=8) :: C2(ndim,ndim)
  real(kind=8) :: C3(ndim,ndim,ndim)

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

  ! C0' = C0 + a_i C1^i + 0.5 a_i a_j C2^{ij} + 1/6 a_i a_j a_k C3^{ijk}
  phi_out = - C0
  do i = 1, ndim
     phi_out = phi_out - a(i)*C1(i)
     do j = 1, ndim
        phi_out = phi_out - 0.5d0*a(i)*a(j)*C2(i,j)
        do k = 1, ndim
           phi_out = phi_out - (1.0d0/6.0d0)*a(i)*a(j)*a(k)*C3(i,j,k)
        end do
     end do
  end do
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
