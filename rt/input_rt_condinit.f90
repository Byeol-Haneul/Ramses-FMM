module input_rt_condinit_module
contains
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_input_rt_condinit(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size

  integer::ilevel
  integer::rID
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_RT_CONDINIT,pst%iUpper+1,input_size,0,ilevel)
     call r_input_rt_condinit(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call input_rt_condinit(pst%s%r,pst%s%g,pst%s%m,ilevel)
  endif
  
end subroutine r_input_rt_condinit
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_rt_condinit(r,g,m,ilevel)
  use amr_parameters, only: ndim,twotondim,dp,nvector
  use rt_parameters,only: nrtvar
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  
  ! Local variables
  integer::igrid,ngrid,ind,idim,nstride,i,ivar
  !integer::l
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(dp),dimension(1:nvector,1:nrtvar)::qq
  real(dp)::dx

  if(m%noct(ilevel)==0)return
  !-----------------------------------------------------------
  ! Set RT ICs from user-defined conditions
  !-----------------------------------------------------------
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
              xx(i,idim)=(2*m%grid(igrid+i-1)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx-m%skip(idim)
           end do
        end do
        ! Call initial condition routine
        call rt_condinit(r,g,xx,qq,dx,ngrid)
        ! Scatter primitive variables to main memory
        do ivar=1,nrtvar
           do i=1,ngrid
              m%grid(igrid+i-1)%rtuold(ind,ivar)=qq(i,ivar)
           end do
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

end subroutine input_rt_condinit
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rt_region_condinit(r,g,x,q,dx,nn)
  use amr_parameters, only:dp, nvector, ndim
  use rt_parameters, only: nrtvar, nrtgroups, rt_c  
  use amr_commons, only: run_t, global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  integer ::nn,idim,igrp
  real(dp)::dx
  real(dp),dimension(1:nvector,1:nrtvar)::q
  real(dp),dimension(1:nvector,1:ndim)::x
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::scale_np, scale_fp, dx_cgs
  !----------------------------------------------------
  ! This routine sets simple pre-defined initial
  ! conditions, like points, squares, etc.
  !----------------------------------------------------
  integer::i,k,group_ind
  real(dp)::vol,rad,weight,xn,yn,zn,en
  ! Units are needed for point regions
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(r,g,scale_np, scale_fp)
  dx_cgs=dx*scale_l

  ! Set (tiny) default values in case n_region=0
  do igrp = 1, nrtgroups
     q(1:nn,1+(igrp-1)*(ndim+1)) = r%smallnp  ! photon densities
     do idim = 1, ndim
        q(1:nn,1+idim+(igrp-1)*(ndim+1)) = 0.0d0 ! photon fluxes
     end do
  end do

  ! Loop over initial conditions regions
  do k=1,r%rt_nregion
     if(r%rt_reg_group(k) .le. 0 .or. r%rt_reg_group(k) .gt. nrtgroups) cycle
     if(r%rt_n_region(k).le.0.0) r%rt_n_region(k)=r%smallnp
     group_ind = 1+(r%rt_reg_group(k)-1)*(ndim+1)
     
     ! For "square" regions only:
     if(r%rt_region_type(k) .eq. 'square')then
        ! Exponent of choosen norm
        en=r%rt_exp_region(k)
        do i=1,nn
           ! Compute position in normalized coordinates
           xn=0.0d0; yn=0.0d0; zn=0.0d0
           xn=2.0d0*abs(x(i,1)-r%rt_reg_x_center(k))/r%rt_reg_length_x(k)
#if NDIM>1
           yn=2.0d0*abs(x(i,2)-r%rt_reg_y_center(k))/r%rt_reg_length_y(k)
#endif
#if NDIM>2
           zn=2.0d0*abs(x(i,3)-r%rt_reg_z_center(k))/r%rt_reg_length_z(k)
#endif
           ! Compute cell "radius" relative to region center
           if(r%rt_exp_region(k)<10)then
              rad=(xn**en+yn**en+zn**en)**(1.0/en)
           else
              rad=max(xn,yn,zn)
           end if
           ! If cell lies within region,
           if(rad<1.0)then
              q(i,group_ind)=r%rt_n_region(k)
              q(i,group_ind+1)=r%rt_u_region(k) * rt_c 
#if NDIM>1 
              q(i,group_ind+2)=r%rt_v_region(k) * rt_c
#endif
#if NDIM>2
              q(i,group_ind+3)=r%rt_w_region(k) * rt_c
#endif
           end if
        end do
     end if
     
     ! For "point" regions only:
     if(r%rt_region_type(k) .eq. 'point')then
        ! Volume elements
        vol=dx_cgs**ndim
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=1.0; yn=1.0; zn=1.0
           xn=max(1.0-abs(x(i,1)-r%rt_reg_x_center(k))/dx,0.0_dp)
#if NDIM>1
           yn=max(1.0-abs(x(i,2)-r%rt_reg_x_center(k))/dx,0.0_dp)
#endif
#if NDIM>2
           zn=max(1.0-abs(x(i,3)-r%rt_reg_x_center(k))/dx,0.0_dp)
#endif
           weight=xn*yn*zn
           if(weight.gt.0) then
              ! If cell lies within CIC cloud, 
              ! Convert photon number to photon number density
              q(i,group_ind) = r%rt_n_region(k)/scale_np *weight/vol 
              q(i,group_ind+1) = r%rt_u_region(k)/scale_np*weight/vol &
                               * rt_c
#if NDIM>1
              q(i,group_ind+2) = r%rt_v_region(k)/scale_np*weight/vol   &
                               * rt_c
#endif
#if NDIM>2
              q(i,group_ind+3) = r%rt_w_region(k)/scale_np *weight/vol  &
                               * rt_c
#endif
           endif
        end do
     end if
  end do

  return
end subroutine rt_region_condinit
!################################################################
!################################################################
!################################################################
!################################################################
end module input_rt_condinit_module