!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_input_hydro_condinit(pst,input_size,output_size,ilevel)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer::ilevel
  
  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_INPUT_HYDRO_CONDINIT,pst%iUpper+1,input_size,output_size,ilevel)
     call r_input_hydro_condinit(pst%pLower,input_size,output_size,ilevel)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     call input_hydro_condinit(pst%s%r,pst%s%g,pst%s%m,ilevel)
  endif
  
end subroutine r_input_hydro_condinit
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_hydro_condinit(r,g,m,ilevel)
  use amr_parameters, only: ndim,twotondim,dp,nvector
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  
  ! Local variables
  integer::igrid,ngrid,ind,idim,nstride,i,ivar
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp)::dx

  if(m%noct(ilevel)==0)return

  !----------------------------------------------------
  ! Compute initial conditions from subroutine condinit
  !----------------------------------------------------
  ! Mesh size at level ilevel in code units
  dx=r%boxlen/2**ilevel
  ! Loop over grids by vector sweeps
  do igrid=m%head(ilevel),m%tail(ilevel),nvector
     ngrid=MIN(nvector,m%tail(ilevel)-igrid+1)
     ! Loop over cells
     do ind=1,twotondim
        ! Compute cell centre position in code units
        do idim=1,ndim
           nstride=2**(idim-1)
           do i=1,ngrid
              xx(i,idim)=(2*m%grid(igrid+i-1)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx
           end do
        end do
        ! Call initial condition routine
        call condinit(r,g,xx,uu,dx,ngrid)
#ifdef HYDRO
        ! Scatter variables to main memory
        do ivar=1,nvar
           do i=1,ngrid
              m%grid(igrid+i-1)%uold(ind,ivar)=uu(i,ivar)
           end do
        end do
#endif
     end do
     ! End loop over cells
  end do
  ! End loop over grids

end subroutine input_hydro_condinit
!################################################################
!################################################################
!################################################################
!################################################################
subroutine region_condinit(r,g,x,q,dx,nn)
  use amr_parameters, only:dp, nvector, ndim
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t, global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  integer ::nn
  real(dp)::dx
  real(dp),dimension(1:nvector,1:nvar)::q
  real(dp),dimension(1:nvector,1:ndim)::x
  !----------------------------------------------------
  ! This routine sets simple pre-defined initial
  ! conditions, like points, squares, etc.
  !----------------------------------------------------
#if NVAR>NDIM+2
  integer::ivar
#endif
  integer::i,k
  real(dp)::vol,rad,weight,xn,yn,zn,en

  ! Set some (tiny) default values in case n_region=0
  q(1:nn,1)=r%smallr
  q(1:nn,2)=0.0d0
#if NDIM>1
  q(1:nn,3)=0.0d0
#endif
#if NDIM>2
  q(1:nn,4)=0.0d0
#endif
  q(1:nn,ndim+2)=r%smallr*r%smallc**2/r%gamma
#if NVAR>NDIM+2
  do ivar=ndim+3,nvar
     q(1:nn,ivar)=0.0d0
  enddo
#endif

  ! Loop over initial conditions regions
  do k=1,r%nregion
     
     ! For "square" regions only:
     if(r%region_type(k) .eq. 'square')then
        ! Exponent of choosen norm
        en=r%exp_region(k)
        do i=1,nn
           ! Compute position in normalized coordinates
           xn=0.0d0; yn=0.0d0; zn=0.0d0
           xn=2.0d0*abs(x(i,1)-r%x_center(k))/r%length_x(k)
#if NDIM>1
           yn=2.0d0*abs(x(i,2)-r%y_center(k))/r%length_y(k)
#endif
#if NDIM>2
           zn=2.0d0*abs(x(i,3)-r%z_center(k))/r%length_z(k)
#endif
           ! Compute cell "radius" relative to region center
           if(r%exp_region(k)<10)then
              rad=(xn**en+yn**en+zn**en)**(1.0/en)
           else
              rad=max(xn,yn,zn)
           end if
           ! If cell lies within region,
           ! REPLACE primitive variables by region values
           if(rad<1.0)then
              q(i,1)=r%d_region(k)
              q(i,2)=r%u_region(k)
#if NDIM>1
              q(i,3)=r%v_region(k)
#endif
#if NDIM>2
              q(i,4)=r%w_region(k)
#endif
              q(i,ndim+2)=r%p_region(k)
#if NENER>0
              do ivar=1,nener
                 q(i,ndim+2+ivar)=r%prad_region(k,ivar)
              enddo
#endif
#if NVAR>NDIM+2+NENER
              do ivar=ndim+3+nener,nvar
                 q(i,ivar)=r%var_region(k,ivar-ndim-2-nener)
              end do
#endif
           end if
        end do
     end if
     
     ! For "point" regions only:
     if(r%region_type(k) .eq. 'point')then
        ! Volume elements
        vol=dx**ndim
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=1.0; yn=1.0; zn=1.0
           xn=max(1.0-abs(x(i,1)-r%x_center(k))/dx,0.0_dp)
#if NDIM>1
           yn=max(1.0-abs(x(i,2)-r%y_center(k))/dx,0.0_dp)
#endif
#if NDIM>2
           zn=max(1.0-abs(x(i,3)-r%z_center(k))/dx,0.0_dp)
#endif
           weight=xn*yn*zn
           ! If cell lies within CIC cloud, 
           ! ADD to primitive variables the region values
           q(i,1)=q(i,1)+r%d_region(k)*weight/vol
           q(i,2)=q(i,2)+r%u_region(k)*weight
#if NDIM>1
           q(i,3)=q(i,3)+r%v_region(k)*weight
#endif
#if NDIM>2
           q(i,4)=q(i,4)+r%w_region(k)*weight
#endif
           q(i,ndim+2)=q(i,ndim+2)+r%p_region(k)*weight/vol
#if NENER>0
           do ivar=1,nener
              q(i,ndim+2+ivar)=q(i,ndim+2+ivar)+r%prad_region(k,ivar)*weight/vol
           enddo
#endif
#if NVAR>NDIM+2+NENER
           do ivar=ndim+3+nener,nvar
              q(i,ivar)=r%var_region(k,ivar-ndim-2-nener)
           end do
#endif
        end do
     end if
  end do

  return
end subroutine region_condinit
!################################################################
!################################################################
!################################################################
!################################################################

