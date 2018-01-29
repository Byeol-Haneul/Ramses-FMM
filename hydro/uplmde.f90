!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdivu(q,div,dx,dy,dz,&
     & iu1,iu2,ju1,ju2,ku1,ku2,&
     & if1,if2,jf1,jf2,kf1,kf2)
  use amr_parameters, only: dp, ndim
  use hydro_parameters, only: nvar
  use const
  implicit none

  real(dp)::dx, dy, dz
  integer::iu1,iu2,ju1,ju2,ku1,ku2
  integer::if1,if2,jf1,jf2,kf1,kf2
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2)::div

  integer::i, j, k
  real(dp)::factorx, factory, factorz
  real(dp)::ux, vy, wz

  factorx=half**(ndim-1)/dx
  factory=half**(ndim-1)/dy
  factorz=half**(ndim-1)/dz

  do k = kf1, kf2
     do j = jf1, jf2
        do i = if1, if2
           ux = zero; vy=zero; wz=zero
#if NDIM>0
           ux=ux+factorx*(q(i,j,k,2) - q(i-1,j,k,2))
#endif
#if NDIM>1           
           ux=ux+factorx*(q(i  ,j-1,k,2) - q(i-1,j-1,k,2))
           vy=vy+factory*(q(i  ,j  ,k,3) - q(i  ,j-1,k,3)+&
                &         q(i-1,j  ,k,3) - q(i-1,j-1,k,3))
#endif
#if NDIM>2
           ux=ux+factorx*(q(i  ,j  ,k-1,2) - q(i-1,j  ,k-1,2)+&
                &         q(i  ,j-1,k-1,2) - q(i-1,j-1,k-1,2))
           vy=vy+factory*(q(i  ,j  ,k-1,3) - q(i  ,j-1,k-1,3)+&
                &         q(i-1,j  ,k-1,3) - q(i-1,j-1,k-1,3))
           wz=wz+factorz*(q(i  ,j  ,k  ,4) - q(i  ,j  ,k-1,4)+&
                &         q(i  ,j-1,k  ,4) - q(i  ,j-1,k-1,4)+&
                &         q(i-1,j  ,k  ,4) - q(i-1,j  ,k-1,4)+&
                &         q(i-1,j-1,k  ,4) - q(i-1,j-1,k-1,4))
#endif
           div(i,j,k) = ux + vy + wz

        end do
     end do
  end do

end subroutine cmpdivu
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine consup(uin,flux,div,dt,&
     & iu1,iu2,ju1,ju2,ku1,ku2,&
     & if1,if2,jf1,jf2,kf1,kf2,difmag)
  use amr_parameters, only: dp, ndim
  use hydro_parameters, only: nvar
  use const
  implicit none

  real(dp)::dt,difmag
  integer::iu1,iu2,ju1,ju2,ku1,ku2
  integer::if1,if2,jf1,jf2,kf1,kf2
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin 
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2)::div 

  integer:: i, j, k, n
  real(dp)::factor
  real(dp)::div1

  factor=half**(ndim-1)

  ! Add diffusive flux where flow is compressing
  do n = 1, nvar
     do k = kf1, MAX(kf1,ku2-2)
        do j = jf1, MAX(jf1, ju2-2) 
           do i = if1, if2
#if NDIM>0
              div1 = factor*div(i,j,k)
#endif
#if NDIM>1
              div1=div1+factor*div(i,j+1,k)
#endif
#if NDIM>2
              div1=div1+factor*(div(i,j,k+1)+div(i,j+1,k+1))
#endif
              div1 = difmag*min(zero,div1)
              flux(i,j,k,n,1) = flux(i,j,k,n,1) + &
                   &  dt*div1*(uin(i,j,k,n) - uin(i-1,j,k,n))
           end do
        end do
     end do

#if NDIM>1
     do k = kf1, MAX(kf1,ku2-2)
        do j = jf1, jf2
           do i = iu1+2, iu2-2
              div1 = zero
              div1 = div1+factor*(div(i,j,k ) + div(i+1,j,k))
#if NDIM>2
              div1=div1+factor*(div(i,j,k+1) + div(i+1,j,k+1))
#endif
              div1 = difmag*min(zero,div1)
              flux(i,j,k,n,2) = flux(i,j,k,n,2) + &
                   &  dt*div1*(uin(i,j,k,n) - uin(i,j-1,k,n))
           end do
        end do
     end do
#endif

#if NDIM>2
     do k = kf1, kf2
        do j = ju1+2, ju2-2 
           do i = iu1+2, iu2-2 
              div1 =  factor*(div(i,j  ,k) + div(i+1,j  ,k) &
                   &        + div(i,j+1,k) + div(i+1,j+1,k))
              div1 = difmag*min(zero,div1)
              flux(i,j,k,n,3) = flux(i,j,k,n,3) + &
                   &  dt*div1*(uin(i,j,k,n) - uin(i,j,k-1,n))
           end do
        end do
     end do
#endif

  end do

end subroutine consup
