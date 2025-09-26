module fmm_taylor
contains
#ifdef FMM
subroutine calc_taylor(xx_j, xx_i, multipoles, taylor_coeff)
  use amr_parameters, only: ndim, multipole_size, taylor_size
  implicit none

  integer :: nq, no, i, j, k, idxC2, idxC3
  real(kind=8), intent(in)  :: xx_j(ndim), xx_i(ndim)
  real(kind=8), intent(in)  :: multipoles(multipole_size)
  real(kind=8), intent(out) :: taylor_coeff(taylor_size)

  ! Multipoles
  real(kind=8) :: M0
  real(kind=8) :: M1(ndim)
  real(kind=8) :: M2(ndim,ndim)

  ! Displacement
  real(kind=8) :: R(ndim), dist

  ! Precompute
  real(kind=8) :: trM, S, MR(ndim), RdotM1
  real(kind=8) :: D0,D1,D2,D3

  !------------------------------------------------------------------
  ! Compute displacement and distance
  !------------------------------------------------------------------
  do i = 1, ndim
     R(i) = xx_i(i) - xx_j(i)
  end do
  dist = sqrt(sum(R(:)**2))
  if (dist == 0.0D0) dist = 1.0D-12

  !------------------------------------------------------------------
  ! Derivatives for g(r) = 1/r
  !------------------------------------------------------------------
  D0 = 1.0D0 / dist
  D1 = -1.0D0 / dist**3
  D2 = 2.0D0 / dist**4
  D3 = -6.0D0 / dist**5

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
  taylor_coeff(1) = M0*D0 - RdotM1*D1 + 0.5D0*trM*D1 + 0.5D0*S*D2

  !------------------------------------------------------------------
  !                               C1
  !------------------------------------------------------------------  
  do i = 1, ndim
     taylor_coeff(1+i) = M0*R(i)*D1 - M1(i)*D1 - R(i)*RdotM1*D2 + 0.5D0*R(i)*trM*D2 + MR(i)*D2 + 0.5D0*R(i)*S*D3
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

end subroutine calc_taylor
!################################################################
!################################################################
!################################################################
!################################################################
subroutine shift_taylor(taylor_in, a, taylor_out)
  use amr_parameters, only: ndim, taylor_size
  implicit none

  real(kind=8), intent(in)  :: taylor_in(taylor_size)
  real(kind=8), intent(in)  :: a(ndim)
  real(kind=8), intent(out) :: taylor_out(taylor_size)

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

  ! C2' = C2 - a_k C3_{ijk}
  do i = 1, ndim
     do j = 1, ndim
        C2p(i,j) = C2(i,j)
        do k = 1, ndim
           C2p(i,j) = C2p(i,j) - a(k)*C3(i,j,k)
        end do
     end do
  end do

  ! C1' = C1 - a_j C2_{ij} + 0.5 a_j a_k C3_{ijk}
  do i = 1, ndim
     C1p(i) = C1(i)
     do j = 1, ndim
        C1p(i) = C1p(i) - a(j)*C2(i,j)
        do k = 1, ndim
           C1p(i) = C1p(i) + 0.5d0*a(j)*a(k)*C3(i,j,k)
        end do
     end do
  end do

  ! C0' = C0 - a_i C1^i + 0.5 a_i a_j C2^{ij} - 1/6 a_i a_j a_k C3^{ijk}
  C0p = C0
  do i = 1, ndim
     C0p = C0p - a(i)*C1(i)
     do j = 1, ndim
        C0p = C0p + 0.5d0*a(i)*a(j)*C2(i,j)
        do k = 1, ndim
           C0p = C0p - (1.0d0/6.0d0)*a(i)*a(j)*a(k)*C3(i,j,k)
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
#endif
end module fmm_taylor
