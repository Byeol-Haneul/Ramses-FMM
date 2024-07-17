!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine interpol_hydro(u1,u2,interpol_var,interpol_type,smallr)
  use amr_parameters, only: dp, ndim,twotondim,twondim
  use hydro_parameters, only: nvar, nener
  implicit none
  integer::interpol_var,interpol_type
  real(dp)::smallr
  real(dp),dimension(0:twondim  ,1:nvar)::u1
  real(dp),dimension(1:twotondim,1:nvar)::u2
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  ! The interpolated variables are:
  ! interpol_var=0: rho, rho u and E
  ! interpol_var=1: rho, rho u and rho epsilon
  ! The interpolation method is:
  ! interpol_type=0 straight injection
  ! interpol_type=1 linear interpolation with MinMod slope
  ! interpol_type=2 linear interpolation with Monotonized Central slope
  ! interpol_type=3 linear interpolation with unlimited Central slope
  !----------------------------------------------------------
#if NENER>0
  integer::irad
#endif
  integer::j,ivar,idim,ind,ix,iy,iz
  real(dp),dimension(1:8,1:3)::xc
  real(dp),dimension(0:twondim)::a
  real(dp),dimension(1:ndim)::w
  real(dp)::ekin,erad

#ifdef HYDRO

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! If necessary, convert fathers total energy into internal energy
  if(interpol_var==1)then
     do j=0,twondim
        ekin=0.0d0
        do idim=1,3
           ekin=ekin+0.5d0*u1(j,idim+1)**2/max(u1(j,1),smallr)
        end do
        erad=0.0d0
#if NENER>0
        do irad=1,nener
           erad=erad+u1(j,5+irad)
        end do
#endif
        u1(j,5)=u1(j,5)-ekin-erad
     end do
  end if

  ! Loop over interpolation variables
  do ivar=1,nvar

     ! Load father variable
     do j=0,twondim
        a(j)=u1(j,ivar)
     end do

     ! Reset gradient
     w(1:ndim)=0.0D0

     ! Compute gradient with chosen limiter
     if(interpol_type==1)call compute_limiter_minmod(a,w)
     if(interpol_type==2)call compute_limiter_central(a,w)
     if(interpol_type==3)call compute_central(a,w)

     ! Interpolate over children cells
     do ind=1,twotondim
        u2(ind,ivar)=a(0)
        do idim=1,ndim
           u2(ind,ivar)=u2(ind,ivar)+w(idim)*xc(ind,idim)
        end do
     end do

  end do
  ! End loop over variables
  
  ! If necessary, convert children internal energy into total energy
  if(interpol_var==1)then
     do ind=1,twotondim
        ekin=0.0d0
        do idim=1,3
           ekin=ekin+0.5d0*u2(ind,idim+1)**2/max(u2(ind,1),smallr)
        end do
        erad=0.0d0
#if NENER>0
        do irad=1,nener
           erad=erad+u2(ind,5+irad)
        end do
#endif
        u2(ind,5)=u2(ind,5)+ekin+erad
     end do
  end if

#endif

end subroutine interpol_hydro
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine interpol_mhd(u1,u2,B1,B2,B3,refined,interpol_var,interpol_type,smallr)
  use amr_parameters, only: dp,ndim,twotondim,twondim
  use hydro_parameters, only: nvar, nener
  implicit none
  integer::interpol_var,interpol_type
  real(dp)::smallr
  real(dp),dimension(0:twondim,1:nvar)::u1
  real(dp),dimension(0:twondim,1:6)::B1
  real(dp),dimension(1:twotondim,1:nvar)::u2
  real(dp),dimension(1:twotondim,1:6)::B2
  real(dp),dimension(1:twondim,1:twotondim,1:6)::B3
  logical,dimension(1:twondim)::refined
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  ! The interpolated variables are:
  ! interpol_var=0: rho, rho u and E
  ! interpol_var=1: rho, rho u and rho epsilon
  ! The interpolation method is:
  ! interpol_type=0 straight injection
  ! interpol_type=1 linear interpolation with MinMod slope
  ! interpol_type=2 linear interpolation with Monotonized Central slope
  ! interpol_type=3 linear interpolation without limiters
  !----------------------------------------------------------
#if NENER>0
  integer::irad
#endif
  integer::i,j,ivar,idim,ind,ix,iy,iz
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(0:twondim)::a
  real(dp),dimension(1:ndim)::w
  real(dp)::ekin,emag,erad

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  !---------------------------------------------------------------
  ! If necessary, convert father total energy into internal energy
  !---------------------------------------------------------------
  if(interpol_var==1)then
     do j=0,twondim
        ekin=0.0d0
        do idim=1,3
           ekin=ekin+0.5d0*u1(j,idim+1)**2/max(u1(j,1),smallr)
        end do
        emag=0.0d0
#ifdef MHD
        do idim=1,3
           emag=emag+0.125d0*(B1(j,idim)+B1(j,idim+3))**2
        end do
#endif
        erad=0.0d0
#if NENER>0
        do irad=1,nener
           erad=erad+u1(j,5+irad)
        end do
#endif
        u1(j,5)=u1(j,5)-ekin-emag-erad
     end do
  end if


  !------------------------------------------------
  ! Loop over cell-centered interpolation variables
  !------------------------------------------------
  do ivar=1,nvar+3-ndim

     ! Load father variable
     if(ivar<=nvar)then
        do j=0,twondim
           a(j)=u1(j,ivar)
        end do
     else
#if NDIM==1
        do j=0,twondim
           a(j)=B1(j,ivar-nvar+1)
        end do
#endif
#if NDIM==2
        do j=0,twondim
           a(j)=B1(j,ivar-nvar+2)
        end do
#endif
     endif

     ! Reset gradient
     w(1:ndim)=0.0D0

     ! Compute gradient with chosen limiter
     if(interpol_type==1)call compute_limiter_minmod(a,w)
     if(interpol_type==2)call compute_limiter_central(a,w)
     if(interpol_type==3)call compute_central(a,w)

     ! Interpolate over children cells
     if(ivar<=nvar)then
        do ind=1,twotondim
           u2(ind,ivar)=a(0)
           do idim=1,ndim
              u2(ind,ivar)=u2(ind,ivar)+w(idim)*xc(ind,idim)
           end do
        end do
     else
#if NDIM==1
        do ind=1,twotondim
           B2(ind,ivar-nvar+1)=a(0)
           do idim=1,ndim
              B2(ind,ivar-nvar+1)=B2(ind,ivar-nvar+1)+w(idim)*xc(ind,idim)
           end do
        end do
#endif
#if NDIM==2
        do ind=1,twotondim
           B2(ind,ivar-nvar+2)=a(0)
           do idim=1,ndim
              B2(ind,ivar-nvar+2)=B2(ind,ivar-nvar+2)+w(idim)*xc(ind,idim)
           end do
        end do
#endif
     endif
  end do
  ! End loop over cell-centered variables

  ! Update cell centered magnetic field also in redundant array
#if NDIM<3
  do ind=1,twotondim
     B2(ind,6)=B2(ind,3)
  end do
#endif
#if NDIM<2
  do ind=1,twotondim
     B2(ind,5)=B2(ind,2)
  end do
#endif

  !----------------------------------------
  ! Interpolate face-centered MHD variables
  !----------------------------------------
  call interpol_mag(B1,B2,B3,refined,interpol_type)

  !-----------------------------------------------------------------
  ! If necessary, convert children internal energy into total energy
  !-----------------------------------------------------------------
  if(interpol_var==1)then
     do ind=1,twotondim
        ekin=0.0d0
        do idim=1,3
           ekin=ekin+0.5d0*u2(ind,idim+1)**2/max(u2(ind,1),smallr)
        end do
        emag=0.0d0
#ifdef MHD
        do idim=1,3
           emag=emag+0.125d0*(B2(ind,idim)+B2(ind,idim+3))**2
        end do
#endif
        erad=0.0d0
#if NENER>0
        do irad=1,nener
           erad=erad+u2(ind,5+irad)
        end do
#endif
        u2(ind,5)=u2(ind,5)+ekin+emag+erad
     end do
  end if

end subroutine interpol_mhd
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_mag(B1,B2,B3,refined,interpol_type)
  use amr_parameters, only: dp, twondim, twotondim
  implicit none
  real(dp),dimension(0:twondim,1:6)::B1
  real(dp),dimension(1:twotondim,1:6)::B2
  real(dp),dimension(1:twondim,1:twotondim,1:6)::B3
  logical ,dimension(1:twondim)::refined
  integer::interpol_type
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  ! The interpolated variables are Bx, By and Bz.
  ! Divergence free is garanteed as in Balsara (2001) JCP, 174, 614.
  !----------------------------------------------------------
  integer::i,j,k,ind,imax,jmax,kmax
  real(dp),dimension(-1:1,0:1,0:1),save::u
  real(dp),dimension(0:1,-1:1,0:1),save::v
  real(dp),dimension(0:1,0:1,-1:1),save::w

  imax=1; jmax=0; kmax=0
#if NDIM>1
  jmax=1
#endif
#if NDIM>2
  kmax=1
#endif

  ! Compute interpolated fine B over coarse side faces
  call interpol_faces(B1,u,v,w,interpol_type)

  ! Get fine B from refined faces, if any
  call copy_from_refined_faces(B3,refined,u,v,w)

  ! Compute interpolated fine B inside coarse cell.
  call cmp_central_faces(u,v,w)

  ! Scatter results
  do i=0,imax
  do j=0,jmax
  do k=0,kmax
     ind=1+i+2*j+4*k
     B2(ind,1)=u(i-1,j,k)
     B2(ind,4)=u(i,j,k)
#if NDIM>1
     B2(ind,2)=v(i,j-1,k)
     B2(ind,5)=v(i,j,k)
#endif
#if NDIM>2
     B2(ind,3)=w(i,j,k-1)
     B2(ind,6)=w(i,j,k)
#endif
  end do
  end do
  end do

end subroutine interpol_mag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_faces(B1,u,v,w,interpol_type)
  use amr_parameters, only: dp, twondim
  implicit none
  real(dp),dimension(0:twondim,1:6)::B1
  real(dp),dimension(-1:1,0:1,0:1)::u
  real(dp),dimension(0:1,-1:1,0:1)::v
  real(dp),dimension(0:1,0:1,-1:1)::w
  integer::interpol_type
  ! TVD interpolation from coarse faces
  ! interpol_mag_type=0: straight injection
  ! interpol_mag_type=1: linear interpolation with MinMod slope
  ! interpol_mag_type=2: linear interpolation with Monotonized Central slope
  ! interpol_mag_type=3: linear interpolation without limiters
  integer::i,j,k,imax,jmax,kmax
  real(dp),dimension(0:4)::b
  real(dp),dimension(1:2)::s

  imax=1; jmax=0; kmax=0
#if NDIM>1
  jmax=1
#endif
#if NDIM>2
  kmax=1
#endif

  ! Left face along direction x (interpolate Bx)
  b(0) = B1(0,1)
#if NDIM>1
  b(1) = B1(3,1)
  b(2) = B1(4,1)
#endif
#if NDIM>2
  b(3) = B1(5,1)
  b(4) = B1(6,1)
#endif
  s(1:2) = 0.0
#if NDIM==2
  if(interpol_type>0)call compute_1d_tvd(b,s,interpol_type)
#endif
#if NDIM==3
  if(interpol_type>0)call compute_2d_tvd(b,s,interpol_type)
#endif
  do j = 0,jmax
  do k = 0,kmax
     u(-1,j,k) = b(0) + 0.5*s(1)*(dble(j)-0.5) + 0.5*s(2)*(dble(k)-0.5)
  end do
  end do

  ! Right face along direction x (interpolate Bx)
  b(0) = B1(0,4)
#if NDIM>1
  b(1) = B1(3,4)
  b(2) = B1(4,4)
#endif
#if NDIM>2
  b(3) = B1(5,4)
  b(4) = B1(6,4)
#endif

  s(1:2) = 0.0
#if NDIM==2
  if(interpol_type>0)call compute_1d_tvd(b,s,interpol_type)
#endif
#if NDIM==3
  if(interpol_type>0)call compute_2d_tvd(b,s,interpol_type)
#endif
  do j = 0,jmax
  do k = 0,kmax
     u(+1,j,k) = b(0) + 0.5*s(1)*(dble(j)-0.5) + 0.5*s(2)*(dble(k)-0.5)
  end do
  end do

#if NDIM>1
  ! Left face along direction y (interpolate By)
  b(0) = B1(0,2)
  b(1) = B1(1,2)
  b(2) = B1(2,2)
#if NDIM>2
  b(3) = B1(5,2)
  b(4) = B1(6,2)
#endif

  s(1:2) = 0.0
#if NDIM==2
  if(interpol_type>0)call compute_1d_tvd(b,s,interpol_type)
#endif
#if NDIM==3
  if(interpol_type>0)call compute_2d_tvd(b,s,interpol_type)
#endif
  do i = 0,imax
  do k = 0,kmax
     v(i,-1,k) = b(0) + 0.5*s(1)*(dble(i)-0.5) + 0.5*s(2)*(dble(k)-0.5)
  end do
  end do

  ! Right face along direction y (interpolate By)
  b(0) = B1(0,5)
  b(1) = B1(1,5)
  b(2) = B1(2,5)
#if NDIM>2
  b(3) = B1(5,5)
  b(4) = B1(6,5)
#endif

  s(1:2) = 0.0
#if NDIM==2
  if(interpol_type>0)call compute_1d_tvd(b,s,interpol_type)
#endif
#if NDIM==3
  if(interpol_type>0)call compute_2d_tvd(b,s,interpol_type)
#endif
  do i = 0,imax
  do k = 0,kmax
     v(i,+1,k) = b(0) + 0.5*s(1)*(dble(i)-0.5) + 0.5*s(2)*(dble(k)-0.5)
  end do
  end do
#endif

#if NDIM>2
  ! Left face along direction z (interpolate Bz)
  b(0) = B1(0,3)
  b(1) = B1(1,3)
  b(2) = B1(2,3)
  b(3) = B1(3,3)
  b(4) = B1(4,3)

  s(1:2) = 0.0
  if(interpol_type>0)call compute_2d_tvd(b,s,interpol_type)
  do i = 0,imax
  do j = 0,jmax
     w(i,j,-1) = b(0) + 0.5*s(1)*(dble(i)-0.5) + 0.5*s(2)*(dble(j)-0.5)
  end do
  end do

  ! Right face along direction z (interpolate Bz)
  b(0) = B1(0,6)
  b(1) = B1(1,6)
  b(2) = B1(2,6)
  b(3) = B1(3,6)
  b(4) = B1(4,6)

  s(1:2) = 0.0
  if(interpol_type>0)call compute_2d_tvd(b,s,interpol_type)
  do i = 0,imax
  do j = 0,jmax
     w(i,j,+1) = b(0) + 0.5*s(1)*(dble(i)-0.5) + 0.5*s(2)*(dble(j)-0.5)
  end do
  end do
#endif

end subroutine interpol_faces
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine copy_from_refined_faces(B3,refined,u,v,w)
  use amr_parameters,only:dp,twondim,twotondim
  implicit none
  logical,dimension(1:twondim)::refined
  real(dp),dimension(1:twondim,1:twotondim,1:6)::B3
  real(dp),dimension(-1:1,0:1,0:1)::u
  real(dp),dimension(0:1,-1:1,0:1)::v
  real(dp),dimension(0:1,0:1,-1:1)::w

  ! TVD interpolation from coarse faces
  integer::i,j,k,ind,imax,jmax,kmax

  imax=1; jmax=0; kmax=0
#if NDIM>1
  jmax=1
#endif
#if NDIM>2
  kmax=1
#endif

  ! Left face along direction x (interpolate Bx)
  do j = 0,jmax
  do k = 0,kmax
     ind = 1+1+j*2+k*4
     if(refined(1))then
        u(-1,j,k) = B3(1,ind,4)
     end if
  end do
  end do

  ! Right face along direction x (interpolate Bx)
  do j = 0,jmax
  do k = 0,kmax
     ind = 1+0+j*2+k*4
     if(refined(2))then
        u(+1,j,k) = B3(2,ind,1)
     end if
  end do
  end do

#if NDIM>1
  ! Left face along direction y (interpolate By)
  do i = 0,imax
  do k = 0,kmax
     ind = 1+i+1*2+k*4
     if(refined(3))then
        v(i,-1,k) = B3(3,ind,5)
     end if
  end do
  end do

  ! Right face along direction y (interpolate By)
  do i = 0,imax
  do k = 0,kmax
     ind = 1+i+0*2+k*4
     if(refined(4))then
        v(i,+1,k) = B3(4,ind,2)
     end if
  end do
  end do
#endif

#if NDIM>2
  ! Left face along direction z (interpolate Bz)
  do i = 0,imax
  do j = 0,kmax
     ind = 1+i+j*2+1*4
     if(refined(5))then
        w(i,j,-1) = B3(5,ind,6)
     end if
  end do
  end do

  ! Right face along direction z (interpolate Bz)
  do i = 0,imax
  do j = 0,kmax
     ind = 1+i+j*2+0*4
     if(refined(6))then
        w(i,j,+1) = B3(6,ind,3)
     end if
  end do
  end do
#endif

end subroutine copy_from_refined_faces
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_central_faces(u,v,w)
  use amr_parameters,only:dp
  implicit none
  real(dp),dimension(-1:1,0:1,0:1)::u
  real(dp),dimension(0:1,-1:1,0:1)::v
  real(dp),dimension(0:1,0:1,-1:1)::w

  integer::i,j,k,ii,jj,kk,imax,jmax,kmax
  real(dp)::UXX,VYY,WZZ,UXYZ,VXYZ,WXYZ

  imax=1; jmax=0; kmax=0
#if NDIM>1
  jmax=1
#endif
#if NDIM>2
  kmax=1
#endif

  UXX = 0.0_dp
  VYY = 0.0_dp
  WZZ = 0.0_dp
  UXYZ = 0.0_dp
  VXYZ = 0.0_dp
  WXYZ = 0.0_dp

#if NDIM==2
  do i = 0,imax
  do j = 0,jmax
  do k = 0,kmax
     ii = 2*i-1
     jj = 2*j-1
     UXX = UXX + (ii*jj*v(i,jj,k))*0.25
     VYY = VYY + (ii*jj*u(ii,j,k))*0.25
  enddo
  enddo
  enddo
#endif

#if NDIM==3
  do i = 0,imax
  do j = 0,jmax
  do k = 0,kmax
     ii = 2*i-1
     jj = 2*j-1
     kk = 2*k-1
     UXX = UXX + (ii*jj*v(i,jj,k)+ii*kk*w(i,j,kk))*0.125
     VYY = VYY + (jj*kk*w(i,j,kk)+ii*jj*u(ii,j,k))*0.125
     WZZ = WZZ + (ii*kk*u(ii,j,k)+jj*kk*v(i,jj,k))*0.125
     UXYZ = UXYZ + (ii*jj*kk*u(ii,j,k))*0.125
     VXYZ = VXYZ + (ii*jj*kk*v(i,jj,k))*0.125
     WXYZ = WXYZ + (ii*jj*kk*w(i,j,kk))*0.125
  enddo
  enddo
  enddo
#endif

#if NDIM==1
  ! Bx on central faces
  do j = 0,jmax
  do k = 0,kmax
     u(0,j,k) = 0.5*(u(-1,j,k)+u(+1,j,k))
  enddo
  enddo
#endif

#if NDIM==2
  ! Bx on central faces
  do j = 0,jmax
  do k = 0,kmax
     u(0,j,k) = 0.5*(u(-1,j,k)+u(+1,j,k)) + UXX
  enddo
  enddo
  do i = 0,imax
  do k = 0,kmax
     v(i,0,k) = 0.5*(v(i,-1,k)+v(i,+1,k)) + VYY
  enddo
  enddo
#endif

#if NDIM==3
  ! Bx on central faces
  do j = 0,jmax
  do k = 0,kmax
     u(0,j,k) = 0.5*(u(-1,j,k)+u(+1,j,k)) + UXX     &
          &   + (dble(k)-0.5)*VXYZ + (dble(j)-0.5)*WXYZ
  enddo
  enddo
  do i = 0,imax
  do k = 0,kmax
     v(i,0,k) = 0.5*(v(i,-1,k)+v(i,+1,k)) + VYY     &
          &   + (dble(i)-0.5)*WXYZ + (dble(k)-0.5)*UXYZ
  enddo
  enddo
  do i = 0,imax
  do j = 0,jmax
     w(i,j,0) = 0.5*(w(i,j,-1)+w(i,j,+1)) + WZZ     &
          &   + (dble(j)-0.5)*UXYZ + (dble(i)-0.5)*VXYZ
  enddo
  enddo
#endif

end subroutine cmp_central_faces
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_2d_tvd(b,s,interpol_type)
  use amr_parameters,only:dp
  use const
  implicit none
  real(dp),dimension(0:4)::b
  real(dp),dimension(1:2)::s
  integer::interpol_type

  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop

  if(interpol_type==3)then
     dlft = half*(b(0) - b(1))
     drgt = half*(b(2) - b(0))
     s(1) = dlft+drgt
     dlft = half*(b(0) - b(3))
     drgt = half*(b(4) - b(0))
     s(2) = dlft+drgt
     return
  endif

  dlft = interpol_type*(b(0) - b(1))
  drgt = interpol_type*(b(2) - b(0))
  dcen = half*(dlft+drgt)/interpol_type
  dsgn = sign(one, dcen)
  slop = min(abs(dlft),abs(drgt))
  dlim = slop
  if((dlft*drgt)<=zero)dlim=zero
  s(1) = dsgn*min(dlim,abs(dcen))

  dlft = interpol_type*(b(0) - b(3))
  drgt = interpol_type*(b(4) - b(0))
  dcen = half*(dlft+drgt)/interpol_type
  dsgn = sign(one, dcen)
  slop = min(abs(dlft),abs(drgt))
  dlim = slop
  if((dlft*drgt)<=zero)dlim=zero
  s(2) = dsgn*min(dlim,abs(dcen))

end subroutine compute_2d_tvd
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_1d_tvd(b,s,interpol_type)
  use amr_parameters,only:dp
  use const
  implicit none
  real(dp),dimension(0:4)::b
  real(dp),dimension(1:2)::s
  integer::interpol_type

  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop

  if(interpol_type==3)then
     dlft = half*(b(0) - b(1))
     drgt = half*(b(2) - b(0))
     s(1) = dlft+drgt
     return
  endif
  dlft = interpol_type*(b(0) - b(1))
  drgt = interpol_type*(b(2) - b(0))
  dcen = half*(dlft+drgt)/interpol_type
  dsgn = sign(one, dcen)
  slop = min(abs(dlft),abs(drgt))
  dlim = slop
  if((dlft*drgt)<=zero)dlim=zero
  s(1) = dsgn*min(dlim,abs(dcen))

end subroutine compute_1d_tvd
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_limiter_minmod(a,w)
  use amr_parameters, only: dp, ndim,twondim,twotondim
  implicit none
  real(dp),dimension(0:twondim)::a
  real(dp),dimension(1:ndim)::w
  !---------------
  ! MinMod slope
  !---------------
  integer::idim
  real(dp)::diff_left,diff_right,minmod

  do idim=1,ndim
     diff_left=0.5*(a(2*idim)-a(0))
     diff_right=0.5*(a(0)-a(2*idim-1))
     if(diff_left*diff_right<=0.0)then
        minmod=0.0
     else
        minmod=MIN(ABS(diff_left),ABS(diff_right)) &
             &   *diff_left/ABS(diff_left)
     end if
     w(idim)=minmod
  end do

end subroutine compute_limiter_minmod
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_limiter_central(a,w)
  use amr_parameters, only: dp, ndim,twondim,twotondim
  implicit none
  real(dp),dimension(0:twondim)::a
  real(dp),dimension(1:ndim)::w
  !---------------------------
  ! Monotonized Central slope
  !---------------------------
  integer::j,idim,ind,ix,iy,iz
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp)::xxc
  real(dp),dimension(1:twotondim)::ac
  real(dp)::corner,kernel,diff_corner,diff_kernel
  real(dp)::max_limiter,min_limiter,limiter

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! Second order central slope
  do idim=1,ndim
     w(idim)=0.25D0*(a(2*idim)-a(2*idim-1))
  end do

  ! Compute corner interpolated values
  do ind=1,twotondim
     ac(ind)=a(0)
  end do
  do idim=1,ndim
     do ind=1,twotondim
        xxc = xc(ind,idim)
        corner=ac(ind)+2.D0*w(idim)*xxc
        ac(ind)=corner
     end do
  end do

  ! Compute max of corners
  corner=ac(1)
  do j=2,twotondim
     corner=MAX(corner,ac(j))
  end do

  ! Compute max of gradient kernel
  kernel=a(1)
  do j=2,twondim
     kernel=MAX(kernel,a(j))
  end do

  ! Compute differences
  diff_kernel=a(0)-kernel
  diff_corner=a(0)-corner

  ! Compute max_limiter
  max_limiter=0.0D0
  if(diff_kernel*diff_corner > 0.0D0)then
     max_limiter=MIN(1.0_dp,diff_kernel/diff_corner)
  end if

  ! Compute min of corners
  corner=ac(1)
  do j=2,twotondim
     corner=MIN(corner,ac(j))
  end do

  ! Compute min of gradient kernel
  kernel=a(1)
  do j=2,twondim
     kernel=MIN(kernel,a(j))
  end do

  ! Compute differences
  diff_kernel=a(0)-kernel
  diff_corner=a(0)-corner

  ! Compute min_limiter
  min_limiter=0.0D0
  if(diff_kernel*diff_corner > 0.0D0)then
     min_limiter=MIN(1.0_dp,diff_kernel/diff_corner)
  end if

  ! Compute limiter
  limiter=MIN(min_limiter,max_limiter)

  ! Correct gradient with limiter
  do idim=1,ndim
     w(idim)=w(idim)*limiter
  end do

end subroutine compute_limiter_central
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_central(a,w)
  use amr_parameters, only: dp, ndim,twondim,twotondim
  implicit none
  real(dp),dimension(0:twondim)::a
  real(dp),dimension(1:ndim)::w
  !---------------------------
  ! Unlimited Central slope
  !---------------------------
  integer::idim

  ! Second order central slope
  do idim=1,ndim
     w(idim)=0.25D0*(a(2*idim)-a(2*idim-1))
  end do

end subroutine compute_central

#ifdef RT
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine interpol_rt(u1,u2,interpol_type,smallnp)
  use amr_parameters, only: dp, ndim,twotondim,twondim
  use rt_parameters, only: nrtvar
  implicit none
  integer::interpol_type
  real(dp)::smallnp
  real(dp),dimension(0:twondim  ,1:nrtvar)::u1
  real(dp),dimension(1:twotondim,1:nrtvar)::u2
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  ! The interpolation method is:
  ! interpol_type=0 straight injection
  ! interpol_type=1 linear interpolation with MinMod slope
  ! interpol_type=2 linear interpolation with Monotonized Central slope
  ! interpol_type=3 linear interpolation with unlimited Central slope
  !----------------------------------------------------------
  integer::j,ivar,idim,ind,ix,iy,iz
  real(dp),dimension(1:8,1:3)::xc
  real(dp),dimension(0:twondim)::a
  real(dp),dimension(1:ndim)::w

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! Loop over interpolation variables
  do ivar=1,nrtvar

     ! Load father variable
     do j=0,twondim
        a(j)=u1(j,ivar)
     end do

     ! Reset gradient
     w(1:ndim)=0.0D0

     ! Compute gradient with chosen limiter
     if(interpol_type==1)call compute_limiter_minmod(a,w)
     if(interpol_type==2)call compute_limiter_central(a,w)
     if(interpol_type==3)call compute_central(a,w)

     ! Interpolate over children cells
     do ind=1,twotondim
        u2(ind,ivar)=a(0)
        do idim=1,ndim
           u2(ind,ivar)=u2(ind,ivar)+w(idim)*xc(ind,idim)
        end do
     end do

  end do
  ! End loop over variables
  

end subroutine interpol_rt

#endif