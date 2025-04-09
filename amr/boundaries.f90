module boundaries
contains
subroutine get_bound(s,ix,ilevel,ibound)
  use amr_parameters, only: ndim,twotondim
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  integer(kind=8),dimension(1:ndim)::ix
  integer::ibound,ilevel

  logical::in_bound
  integer::i,idim
  integer(kind=4)::bound_ckey_min  ! Min. Cartesian key per level for the boundary region
  integer(kind=4)::bound_ckey_max  ! Max. Cartesian key per level for the boundary region

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

    do i=1,r%nbound
       in_bound = .true.
       if(r%bound_xmin(i).GE.0)then
          bound_ckey_min=r%bound_xmin(i)*2**(ilevel-r%levelmin)
       else
          bound_ckey_min=(m%ckey_max(r%levelmin)+r%bound_xmin(i))*2**(ilevel-r%levelmin)
       endif
       if(r%bound_xmax(i).LE.0)then
          bound_ckey_max=(m%ckey_max(r%levelmin)+r%bound_xmax(i))*2**(ilevel-r%levelmin)-1
       else
          bound_ckey_max=r%bound_xmax(i)*2**(ilevel-r%levelmin)-1
       endif
       in_bound = in_bound .and. ix(1) .ge. bound_ckey_min .and. ix(1) .le. bound_ckey_max
#if NDIM>1
       if(r%bound_ymin(i).GE.0)then
          bound_ckey_min=r%bound_ymin(i)*2**(ilevel-r%levelmin)
       else
          bound_ckey_min=(m%ckey_max(r%levelmin)+r%bound_ymin(i))*2**(ilevel-r%levelmin)
       endif
       if(r%bound_ymax(i).LE.0)then
          bound_ckey_max=(m%ckey_max(r%levelmin)+r%bound_ymax(i))*2**(ilevel-r%levelmin)-1
       else
          bound_ckey_max=r%bound_ymax(i)*2**(ilevel-r%levelmin)-1
       endif
       in_bound = in_bound .and. ix(2) .ge. bound_ckey_min .and. ix(2) .le. bound_ckey_max
#endif
#if NDIM>2
       if(r%bound_zmin(i).GE.0)then
          bound_ckey_min=r%bound_zmin(i)*2**(ilevel-r%levelmin)
       else
          bound_ckey_min=(m%ckey_max(r%levelmin)+r%bound_zmin(i))*2**(ilevel-r%levelmin)
       endif
       if(r%bound_zmax(i).LE.0)then
          bound_ckey_max=(m%ckey_max(r%levelmin)+r%bound_zmax(i))*2**(ilevel-r%levelmin)-1
       else
          bound_ckey_max=r%bound_zmax(i)*2**(ilevel-r%levelmin)-1
       endif
       in_bound = in_bound .and. ix(3) .ge. bound_ckey_min .and. ix(3) .le. bound_ckey_max
#endif
       if(in_bound)then
          ibound=i
          exit
       endif
    end do

  end associate

end subroutine get_bound
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_bound_refine(r,g,m,grid,grid_ref,ibound)
  use amr_parameters, only: ndim, twotondim, dp, nvector
  use hydro_parameters, only: nvar, nener
  use rt_parameters, only: nrtvar, nrtgrp
  use amr_commons, only: run_t, global_t, mesh_t, oct
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(oct)::grid, grid_ref
  integer::ibound
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::scale_np, scale_fp, dx_cgs, dt_cgs

  integer,dimension(1:8,1:3)::ref_right=reshape(&
       & (/ 0,1,0,3,0,5,0,7,&
       &    0,0,1,2,0,0,5,6,&
       &    0,0,0,0,1,2,3,4/),(/8,3/))
  integer,dimension(1:8,1:3)::ref_left=reshape(&
       & (/ 2,0,4,0,6,0,8,0,&
       &    3,4,0,0,7,8,0,0,&
       &    5,6,7,8,0,0,0,0/),(/8,3/))
  integer,dimension(1:8,1:3)::ind1_right=reshape(&
       & (/ 2,1,4,3,6,5,8,7,&
       &    3,4,1,2,7,8,5,6,&
       &    5,6,7,8,1,2,3,4/),(/8,3/))
  integer,dimension(1:8,1:3)::ind1_left=reshape(&
       & (/ 2,1,4,3,6,5,8,7,&
       &    3,4,1,2,7,8,5,6,&
       &    5,6,7,8,1,2,3,4/),(/8,3/))
  integer,dimension(1:8,1:3)::ind2_right=reshape(&
       & (/ 1,1,3,3,5,5,7,7,&
       &    1,2,1,2,5,6,5,6,&
       &    1,2,3,4,1,2,3,4/),(/8,3/))
  integer,dimension(1:8,1:3)::ind2_left=reshape(&
       & (/ 2,2,4,4,6,6,8,8,&
       &    3,4,3,4,7,8,7,8,&
       &    5,6,7,8,5,6,7,8/),(/8,3/))

  integer::idim, ind, ivar, irad
  integer::type, dir, shift, nstride
  real(dp)::reverse, ek_bound
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp),dimension(1:nvector,1:nvar)::qq
  real(dp)::dx,rr,vx,vy,vz,pp,eint,ekin,emag,erad

  type = r%bound_type(ibound)
  dir = r%bound_dir(ibound)
  shift = r%bound_shift(ibound)

  ! Set refinement map
  if(shift==+1)then
     do ind=1,twotondim
        if(ref_right(ind,dir)>0)then
           grid%refined(ind)=grid_ref%refined(ref_right(ind,dir))
        else
           grid%refined(ind)=.false.
        endif
     end do
  end if
  if(shift==-1)then
     do ind=1,twotondim
        if(ref_left(ind,dir)>0)then
           grid%refined(ind)=grid_ref%refined(ref_left(ind,dir))
        else
           grid%refined(ind)=.false.
        endif
     end do
  end if

#ifdef HYDRO

  ! Reflexive BC
  if(type == 1)then

     if(shift==+1)then
        do ivar=1,nvar
           reverse=1
           if(ivar==1+dir)reverse=-1
           do ind=1,twotondim
              grid%uold(ind,ivar)=grid_ref%uold(ind1_right(ind,dir),ivar)*reverse
           end do
        end do
#ifdef MHD
#if NDIM==1
        do ind=1,twotondim
           grid%bold(ind,1)=r%A_ave
           grid%bold(ind,4)=r%A_ave
           grid%bold(ind,2)=grid_ref%bold(ind1_right(ind,dir),2)
           grid%bold(ind,5)=grid_ref%bold(ind1_right(ind,dir),5)
           grid%bold(ind,3)=grid_ref%bold(ind1_right(ind,dir),3)
           grid%bold(ind,6)=grid_ref%bold(ind1_right(ind,dir),6)
        end do
#endif
#if NDIM==2
        do ind=1,twotondim
           grid%bold(ind,3)=grid_ref%bold(ind1_right(ind,dir),3)
           grid%bold(ind,6)=grid_ref%bold(ind1_right(ind,dir),6)
        end do
#endif
#endif
#ifdef RT
        do ivar=1,nrtvar
           reverse=1
           if(mod(ivar-1, ndim+1) == dir) reverse=-1
           do ind=1,twotondim
              grid%rtuold(ind,ivar)=grid_ref%rtuold(ind1_right(ind,dir),ivar)*reverse
           end do
        end do
#endif
     endif

     if(shift==-1)then
        do ivar=1,nvar
           reverse=1
           if(ivar==1+dir)reverse=-1
           do ind=1,twotondim
              grid%uold(ind,ivar)=grid_ref%uold(ind1_left(ind,dir),ivar)*reverse
           end do
        end do
#ifdef MHD
#if NDIM==1
        do ind=1,twotondim
           grid%bold(ind,1)=r%A_ave
           grid%bold(ind,4)=r%A_ave
           grid%bold(ind,2)=grid_ref%bold(ind1_left(ind,dir),2)
           grid%bold(ind,5)=grid_ref%bold(ind1_left(ind,dir),5)
           grid%bold(ind,3)=grid_ref%bold(ind1_left(ind,dir),3)
           grid%bold(ind,6)=grid_ref%bold(ind1_left(ind,dir),6)
        end do
#endif
#if NDIM==2
        do ind=1,twotondim
           grid%bold(ind,3)=grid_ref%bold(ind1_left(ind,dir),3)
           grid%bold(ind,6)=grid_ref%bold(ind1_left(ind,dir),6)
        end do
#endif
#endif
#ifdef RT
        do ivar=1,nrtvar
           reverse=1
           if(mod(ivar-1, ndim+1) == dir) reverse=-1
           do ind=1,twotondim
              grid%rtuold(ind,ivar)=grid_ref%rtuold(ind1_left(ind,dir),ivar)*reverse
           end do
        end do
#endif
     endif

  endif

  ! Zero gradient BC
  if(type == 2)then

     if(shift==+1)then
        do ivar=1,nvar
           do ind=1,twotondim
              grid%uold(ind,ivar)=grid_ref%uold(ind2_right(ind,dir),ivar)
           end do
        end do
#ifdef MHD
#if NDIM==1
        do ind=1,twotondim
           grid%bold(ind,1)=r%A_ave
           grid%bold(ind,4)=r%A_ave
           grid%bold(ind,2)=grid_ref%bold(ind2_right(ind,dir),2)
           grid%bold(ind,5)=grid_ref%bold(ind2_right(ind,dir),5)
           grid%bold(ind,3)=grid_ref%bold(ind2_right(ind,dir),3)
           grid%bold(ind,6)=grid_ref%bold(ind2_right(ind,dir),6)
        end do
#endif
#if NDIM==2
        do ind=1,twotondim
           grid%bold(ind,3)=grid_ref%bold(ind2_right(ind,dir),3)
           grid%bold(ind,6)=grid_ref%bold(ind2_right(ind,dir),6)
        end do
#endif
#endif
#ifdef RT
        do ivar=1,nrtvar
           do ind=1,twotondim
              grid%rtuold(ind,ivar)=grid_ref%rtuold(ind2_right(ind,dir),ivar)
           end do
        end do
#endif
     endif

     if(shift==-1)then
        do ivar=1,nvar
           do ind=1,twotondim
              grid%uold(ind,ivar)=grid_ref%uold(ind2_left(ind,dir),ivar)
           end do
        end do
#ifdef MHD
#if NDIM==1
        do ind=1,twotondim
           grid%bold(ind,1)=r%A_ave
           grid%bold(ind,4)=r%A_ave
           grid%bold(ind,2)=grid_ref%bold(ind2_left(ind,dir),2)
           grid%bold(ind,5)=grid_ref%bold(ind2_left(ind,dir),5)
           grid%bold(ind,3)=grid_ref%bold(ind2_left(ind,dir),3)
           grid%bold(ind,6)=grid_ref%bold(ind2_left(ind,dir),6)
        end do
#endif
#if NDIM==2
        do ind=1,twotondim
           grid%bold(ind,3)=grid_ref%bold(ind2_left(ind,dir),3)
           grid%bold(ind,6)=grid_ref%bold(ind2_left(ind,dir),6)
        end do
#endif
#endif
#ifdef RT
        do ivar=1,nrtvar
           do ind=1,twotondim
              grid%rtuold(ind,ivar)=grid_ref%rtuold(ind2_left(ind,dir),ivar)
           end do
        end do
#endif
     endif

  endif

  ! Imposed BC from namelist constant values
  if(type == 3)then

     do ind=1,twotondim

        grid%uold(ind,1)=r%d_bound(ibound)
        grid%uold(ind,2)=r%d_bound(ibound)*r%u_bound(ibound)
        ek_bound=0.5d0*r%d_bound(ibound)*r%u_bound(ibound)**2
        grid%uold(ind,3)=r%d_bound(ibound)*r%v_bound(ibound)
        ek_bound=ek_bound+0.5d0*r%d_bound(ibound)*r%v_bound(ibound)**2
        grid%uold(ind,4)=r%d_bound(ibound)*r%w_bound(ibound)
        ek_bound=ek_bound+0.5d0*r%d_bound(ibound)*r%w_bound(ibound)**2
#if NENER>0
        do ivar=1,nener
           grid%uold(ind,5+ivar)=r%prad_bound(ibound,ivar)/(r%gamma_rad(ivar)-1.0d0)
           ek_bound=ek_bound+r%prad_bound(ibound,ivar)/(r%gamma_rad(ivar)-1.0d0)
        enddo
#endif
        grid%uold(ind,5)=r%p_bound(ibound)/(r%gamma-1.0d0)+ek_bound
#if NVAR>5+NENER
        do ivar=6+nener,nvar
           grid%uold(ind,ivar)=r%d_bound(ibound)*r%var_bound(ibound,ivar-5-nener)
        end do
#endif
#ifdef MHD
#if NDIM==1
        grid%bold(ind,1)=r%A_ave
        grid%bold(ind,4)=r%A_ave
        grid%bold(ind,2)=r%B_bound(ibound)
        grid%bold(ind,5)=r%B_bound(ibound)
        grid%bold(ind,3)=r%C_bound(ibound)
        grid%bold(ind,6)=r%C_bound(ibound)
#endif
#if NDIM==2
        grid%bold(ind,3)=r%C_bound(ibound)
        grid%bold(ind,6)=r%C_bound(ibound)
#endif
#endif
     end do
#ifdef RT
     call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     call rt_units(r,g,scale_np,scale_fp)
     do ind=1,twotondim
        do ivar=1,nrtgrp
            grid%rtuold(ind,1+(ivar-1)*(ndim+1))= &
              r%rt_n_bound(ibound,ivar)/g%rt_c(grid%lev)/scale_fp
            grid%rtuold(ind,2+(ivar-1)*(ndim+1))= &
              r%rt_n_bound(ibound,ivar)*r%rt_u_bound(ibound,ivar)/scale_fp
            if(ndim>1) grid%rtuold(ind,3+(ivar-1)*(ndim+1))= &
                r%rt_n_bound(ibound,ivar)*r%rt_v_bound(ibound,ivar)/scale_fp
            if(ndim>2) grid%rtuold(ind,4+(ivar-1)*(ndim+1))= &
                r%rt_n_bound(ibound,ivar)*r%rt_w_bound(ibound,ivar)/scale_fp
        end do
     end do
#endif

  endif

  ! Imposed BC from condinit
  if(type == 4)then

     ! Mesh size at level ilevel in code units
     dx=r%boxlen/2**grid%lev

     do ind=1,twotondim
        do idim=1,ndim
           nstride=2**(idim-1)
           xx(1,idim)=(2*grid%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx-m%skip(idim)
        end do
        ! Call initial condition routine
        call condinit(r,g,xx,qq,dx,1)
        ! Compute kinetic and internal energy densities
        rr=qq(1,1)
        vx=qq(1,2)
        vy=qq(1,3)
        vz=qq(1,4)
        pp=qq(1,5)
        ekin=0.5d0*rr*(vx**2+vy**2+vz**2)
        eint=pp/(r%gamma-1.0)
        emag=0.0d0
        erad=0.0d0
#if NENER>0
        ! Compute non-thermal energy densities
        do irad=1,nener
           uu(1,5+irad)=qq(1,5+irad)/(r%gamma_rad(irad)-1.0d0)
           erad=erad+uu(1,5+irad)
        end do
#endif
#if NVAR>5+NENER
        ! Compute passive scalar density
        do ivar=6+nener,nvar
           uu(1,ivar)=rr*qq(1,ivar)
        enddo
#endif
        ! Compute density
        uu(1,1)=qq(1,1)
        ! Compute momentum density
        do idim=1,3
           uu(1,idim+1)=rr*qq(1,idim+1)
        end do
        ! Compute total fluid energy density
        uu(1,5)=eint+ekin+erad+emag
        ! Convert to conservative variables
        do ivar=1,nvar
           grid%uold(ind,ivar)=uu(1,ivar)
        end do
     end do

  endif

  ! Imposed BC from boundana
  if(type == 5)then

     ! Mesh size at level ilevel in code units
     dx=r%boxlen/2**grid%lev

     do ind=1,twotondim
        do idim=1,ndim
           nstride=2**(idim-1)
           xx(1,idim)=(2*grid%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx-m%skip(idim)
        end do
        ! Call initial condition routine
        call boundana(r,g,xx,uu,dx,ibound,1)
        ! Scatter variables to main memory
        do ivar=1,nvar
           grid%uold(ind,ivar)=uu(1,ivar)
        end do
     end do

  endif

#endif

#ifdef GRAV
  ! Isolated BC
  do idim=1,ndim
     do ind=1,twotondim
        grid%f(ind,idim)=0.0d0
     end do
  end do
  do ind=1,twotondim
     grid%phi(ind)=0.0d0
     grid%phi_old(ind)=0.0d0
  end do
#endif

end subroutine init_bound_refine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_bound_flag(r,g,m,grid,grid_ref,ibound)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t, global_t, mesh_t, oct
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(oct)::grid, grid_ref
  integer::ibound

  integer,dimension(1:8,1:3)::ref_right=reshape(&
       & (/ 0,1,0,3,0,5,0,7,&
       &    0,0,1,2,0,0,5,6,&
       &    0,0,0,0,1,2,3,4/),(/8,3/))
  integer,dimension(1:8,1:3)::ref_left=reshape(&
       & (/ 2,0,4,0,6,0,8,0,&
       &    3,4,0,0,7,8,0,0,&
       &    5,6,7,8,0,0,0,0/),(/8,3/))

  integer::idim, ind
  integer::type, dir, shift

  type = r%bound_type(ibound)
  dir = r%bound_dir(ibound)
  shift = r%bound_shift(ibound)

  ! Set refinement map
  if(shift==+1)then
     do ind=1,twotondim
        if(ref_right(ind,dir)>0)then
           grid%flag1(ind)=grid_ref%flag1(ref_right(ind,dir))
        else
           grid%flag1(ind)=0
        endif
     end do
  end if
  if(shift==-1)then
     do ind=1,twotondim
        if(ref_left(ind,dir)>0)then
           grid%flag1(ind)=grid_ref%flag1(ref_left(ind,dir))
        else
           grid%flag1(ind)=0
        endif
     end do
  end if

end subroutine init_bound_flag
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_bound_flag2(r,g,m,grid,grid_ref,ibound)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t, global_t, mesh_t, oct
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(oct)::grid, grid_ref
  integer::ibound

  integer,dimension(1:8,1:3)::ref_right=reshape(&
       & (/ 0,1,0,3,0,5,0,7,&
       &    0,0,1,2,0,0,5,6,&
       &    0,0,0,0,1,2,3,4/),(/8,3/))
  integer,dimension(1:8,1:3)::ref_left=reshape(&
       & (/ 2,0,4,0,6,0,8,0,&
       &    3,4,0,0,7,8,0,0,&
       &    5,6,7,8,0,0,0,0/),(/8,3/))

  integer::idim, ind
  integer::type, dir, shift

  type = r%bound_type(ibound)
  dir = r%bound_dir(ibound)
  shift = r%bound_shift(ibound)

  ! Set refinement map
  if(shift==+1)then
     do ind=1,twotondim
        if(ref_right(ind,dir)>0)then
           grid%flag2(ind)=grid_ref%flag2(ref_right(ind,dir))
        else
           grid%flag2(ind)=0
        endif
     end do
  end if
  if(shift==-1)then
     do ind=1,twotondim
        if(ref_left(ind,dir)>0)then
           grid%flag2(ind)=grid_ref%flag2(ref_left(ind,dir))
        else
           grid%flag2(ind)=0
        endif
     end do
  end if

end subroutine init_bound_flag2
!################################################################
!################################################################
!################################################################
!################################################################
end module boundaries
