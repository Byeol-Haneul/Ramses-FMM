!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdivu(q,div,dx,dy,dz,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::div

  integer::i, j, k, l
  real(dp)::factorx, factory, factorz
  real(dp),dimension(1:nvector)::ux, vy, wz

  factorx=half**(ndim-1)/dx
  factory=half**(ndim-1)/dy
  factorz=half**(ndim-1)/dz

  do k = kf1, kf2
     do j = jf1, jf2
        do i = if1, if2

           ux = zero; vy=zero; wz=zero

           if(ndim>0)then
              do l=1, ngrid
                 ux(l)=ux(l)+factorx*(q(l,i  ,j  ,k  ,2) - q(l,i-1,j  ,k  ,2))
              end do
           end if

#if NDIM>1           
           if(ndim>1)then
              do l=1, ngrid
                 ux(l)=ux(l)+factorx*(q(l,i  ,j-1,k  ,2) - q(l,i-1,j-1,k  ,2))
                 vy(l)=vy(l)+factory*(q(l,i  ,j  ,k  ,3) - q(l,i  ,j-1,k  ,3)+&
                      &               q(l,i-1,j  ,k  ,3) - q(l,i-1,j-1,k  ,3))
              end do
           end if
#endif
#if NDIM>2
           if(ndim>2)then
              do l=1, ngrid
                 ux(l)=ux(l)+factorx*(q(l,i  ,j  ,k-1,2) - q(l,i-1,j  ,k-1,2)+&
                      &               q(l,i  ,j-1,k-1,2) - q(l,i-1,j-1,k-1,2))
                 vy(l)=vy(l)+factory*(q(l,i  ,j  ,k-1,3) - q(l,i  ,j-1,k-1,3)+&
                      &               q(l,i-1,j  ,k-1,3) - q(l,i-1,j-1,k-1,3))
                 wz(l)=wz(l)+factorz*(q(l,i  ,j  ,k  ,4) - q(l,i  ,j  ,k-1,4)+&
                      &               q(l,i  ,j-1,k  ,4) - q(l,i  ,j-1,k-1,4)+&
                      &               q(l,i-1,j  ,k  ,4) - q(l,i-1,j  ,k-1,4)+&
                      &               q(l,i-1,j-1,k  ,4) - q(l,i-1,j-1,k-1,4))
              end do
           end if
#endif
           do l=1,ngrid
              div(l,i,j,k) = ux(l) + vy(l) + wz(l)
           end do

        end do
     end do
  end do

end subroutine cmpdivu
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine consup(uin,flux,div,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin 
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::div 

  integer:: i, j, k, l, n
  real(dp)::factor
  real(dp),dimension(1:nvector),save:: div1

  factor=half**(ndim-1)

  ! Add diffusive flux where flow is compressing
  do n = 1, nvar

     do k = kf1, MAX(kf1,ku2-2)
        do j = jf1, MAX(jf1, ju2-2) 
           do i = if1, if2
              div1 = zero
              do l = 1, ngrid
                 div1(l) = factor*div(l,i,j,k)
              end do
#if NDIM>1
              do l = 1, ngrid
                 div1(l)=div1(l)+factor*div(l,i,j+1,k)
              end do
#endif
#if NDIM>2
              do l = 1, ngrid
                 div1(l)=div1(l)+factor*(div(l,i,j,k+1)+div(l,i,j+1,k+1))
              end do
#endif
              do l = 1, ngrid
                 div1(l) = difmag*min(zero,div1(l))
              end do
              do l = 1, ngrid
                 flux(l,i,j,k,n,1) = flux(l,i,j,k,n,1) + &
                      &  dt*div1(l)*(uin(l,i,j,k,n) - uin(l,i-1,j,k,n))
              end do

           end do
        end do
     end do

#if NDIM>1
     do k = kf1, MAX(kf1,ku2-2)
        do j = jf1, jf2
           do i = iu1+2, iu2-2
              div1 = zero
              do l = 1, ngrid
                 div1(l)=div1(l)+factor*(div(l,i,j,k ) + div(l,i+1,j,k))
              end do
#if NDIM>2
              do l = 1, ngrid
                 div1(l)=div1(l)+factor*(div(l,i,j,k+1) + div(l,i+1,j,k+1))
              end do
#endif
              do l = 1, ngrid
                 div1(l) = difmag*min(zero,div1(l))
              end do
              do l = 1, ngrid
                 flux(l,i,j,k,n,2) = flux(l,i,j,k,n,2) + &
                      &  dt*div1(l)*(uin(l,i,j,k,n) - uin(l,i,j-1,k,n))
              end do
           end do
        end do
     end do
#endif

#if NDIM>2
     do k = kf1, kf2
        do j = ju1+2, ju2-2 
           do i = iu1+2, iu2-2 
              do l = 1, ngrid
                 div1(l)=factor*(div(l,i,j  ,k) + div(l,i+1,j  ,k) &
                      &        + div(l,i,j+1,k) + div(l,i+1,j+1,k))
              end do
              do l = 1, ngrid
                 div1(l) = difmag*min(zero,div1(l))
              end do
              do l = 1, ngrid
                 flux(l,i,j,k,n,3) = flux(l,i,j,k,n,3) + &
                      &  dt*div1(l)*(uin(l,i,j,k,n) - uin(l,i,j,k-1,n))
              end do
           end do
        end do
     end do
#endif

  end do

end subroutine consup
