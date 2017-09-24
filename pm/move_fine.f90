!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine kick_drift_part_2(r,g,m,p,ilevel,action_part)
  use amr_parameters, only: dp,ndim,twotondim
  use pm_parameters
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use cache_commons
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  integer::ilevel
  integer::action_part
  !
  !
  real(dp),dimension(1:ndim),save::x,dd,dg
  integer,dimension(1:ndim),save::ig,id
  real(dp),dimension(1:twotondim),save::vol
  integer,dimension(1:ndim,1:twotondim),save::ckey
  integer,dimension(1:twotondim),save::igrid,icell
  integer(kind=8),dimension(0:ndim),save::hash_nbor
  integer::ipart,ind,idim
  integer::parent_cell,get_parent_cell_2
  real(kind=8)::dx_loc,vol_loc,dteff
  real(dp),dimension(1:ndim),save::ff
  logical::ok_level
  
  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Open read-only cache
  call open_cache_2(r,g,m,operation_kick,domain_decompos_amr)

  ! Loop over particles
  do ipart=p%headp(ilevel),p%tailp(ilevel)
     
     ! Rescale particle position at level ilevel
     do idim=1,ndim
        x(idim)=p%xp(ipart,idim)/dx_loc
     end do
     
     ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
     do idim=1,ndim
        dd(idim)=x(idim)+0.5D0
        id(idim)=dd(idim)
        dd(idim)=dd(idim)-id(idim)
        dg(idim)=1.0D0-dd(idim)
        ig(idim)=id(idim)-1
     end do
     
     ! Periodic boundary conditions
     do idim=1,ndim
        if(ig(idim)<0)ig(idim)=m%ckey_max(ilevel+1)-1
        if(id(idim)==m%ckey_max(ilevel+1))id(idim)=0
     enddo
     
     ! Compute cells Cartesian key
#if NDIM==1
     ckey(1,1)=ig(1)
     ckey(1,2)=id(1)
#endif
#if NDIM==2
     ckey(1:2,1)=(/ig(1),ig(2)/)
     ckey(1:2,2)=(/id(1),ig(2)/)
     ckey(1:2,3)=(/ig(1),id(2)/)
     ckey(1:2,4)=(/id(1),id(2)/)
#endif
#if NDIM==3
     ckey(1:3,1)=(/ig(1),ig(2),ig(3)/)
     ckey(1:3,2)=(/id(1),ig(2),ig(3)/)
     ckey(1:3,3)=(/ig(1),id(2),ig(3)/)
     ckey(1:3,4)=(/id(1),id(2),ig(3)/)
     ckey(1:3,5)=(/ig(1),ig(2),id(3)/)
     ckey(1:3,6)=(/id(1),ig(2),id(3)/)
     ckey(1:3,7)=(/ig(1),id(2),id(3)/)
     ckey(1:3,8)=(/id(1),id(2),id(3)/)
#endif
     
     ! Get parent cell at level ilevel using read-only cache
     ok_level=.true.
     igrid=0; icell=0
     hash_nbor(0)=ilevel+1
     do ind=1,twotondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        parent_cell=get_parent_cell_2(r,g,m,hash_nbor,m%grid_dict,.false.,.true.)
        if(parent_cell>0)then
           igrid(ind)=(parent_cell-1)/twotondim+1
           icell(ind)=parent_cell-(igrid(ind)-1)*twotondim
           call lock_cache_2(r,g,m,igrid(ind))
        else
           ok_level=.false.
           exit
        end if
     end do
     do ind=1,twotondim
        call unlock_cache_2(r,g,m,igrid(ind))
     end do

     ! If cloud is not fully inside level ilevel, re-do CIC at coarser level
     if(.not. ok_level)then

        ! Rescale particle position at level ilevel
        do idim=1,ndim
           x(idim)=x(idim)/2.0d0
        end do
        
        ! CIC at level ilevel-1 (dd: right cloud boundary; dg: left cloud boundary)
        do idim=1,ndim
           dd(idim)=x(idim)+0.5D0
           id(idim)=dd(idim)
           dd(idim)=dd(idim)-id(idim)
           dg(idim)=1.0D0-dd(idim)
           ig(idim)=id(idim)-1
        end do
        
        ! Periodic boundary conditions
        do idim=1,ndim
           if(ig(idim)<0)ig(idim)=m%ckey_max(ilevel)-1
           if(id(idim)==m%ckey_max(ilevel))id(idim)=0
        enddo
        
        ! Compute cells Cartesian key
#if NDIM==1
        ckey(1,1)=ig(1)
        ckey(1,2)=id(1)
#endif
#if NDIM==2
        ckey(1:2,1)=(/ig(1),ig(2)/)
        ckey(1:2,2)=(/id(1),ig(2)/)
        ckey(1:2,3)=(/ig(1),id(2)/)
        ckey(1:2,4)=(/id(1),id(2)/)
#endif
#if NDIM==3
        ckey(1:3,1)=(/ig(1),ig(2),ig(3)/)
        ckey(1:3,2)=(/id(1),ig(2),ig(3)/)
        ckey(1:3,3)=(/ig(1),id(2),ig(3)/)
        ckey(1:3,4)=(/id(1),id(2),ig(3)/)
        ckey(1:3,5)=(/ig(1),ig(2),id(3)/)
        ckey(1:3,6)=(/id(1),ig(2),id(3)/)
        ckey(1:3,7)=(/ig(1),id(2),id(3)/)
        ckey(1:3,8)=(/id(1),id(2),id(3)/)
#endif
        
        ! Get parent cell at level ilevel-1 using read-only cache
        ok_level=.true.
        hash_nbor(0)=ilevel
        igrid=0; icell=0
        do ind=1,twotondim
           hash_nbor(1:ndim)=ckey(1:ndim,ind)
           parent_cell=get_parent_cell_2(r,g,m,hash_nbor,m%grid_dict,.false.,.true.)
           if(parent_cell>0)then
              igrid(ind)=(parent_cell-1)/twotondim+1
              icell(ind)=parent_cell-(igrid(ind)-1)*twotondim
              call lock_cache_2(r,g,m,igrid(ind))
           else
              ok_level=.false.
              exit
           end if
        end do
        do ind=1,twotondim
           call unlock_cache_2(r,g,m,igrid(ind))
        end do
     end if
        
     ! Compute cloud volumes
#if NDIM==1
     vol(1)=dg(1)
     vol(2)=dd(1)
#endif
#if NDIM==2
     vol(1)=dg(1)*dg(2)
     vol(2)=dd(1)*dg(2)
     vol(3)=dg(1)*dd(2)
     vol(4)=dd(1)*dd(2)
#endif
#if NDIM==3
     vol(1)=dg(1)*dg(2)*dg(3)
     vol(2)=dd(1)*dg(2)*dg(3)
     vol(3)=dg(1)*dd(2)*dg(3)
     vol(4)=dd(1)*dd(2)*dg(3)
     vol(5)=dg(1)*dg(2)*dd(3)
     vol(6)=dd(1)*dg(2)*dd(3)
     vol(7)=dg(1)*dd(2)*dd(3)
     vol(8)=dd(1)*dd(2)*dd(3)
#endif
     
     ! Gather 3-force
     ff(1:ndim)=0.0
     if(ok_level)then
        do ind=1,twotondim
#ifdef GRAV
           ff(1:ndim)=ff(1:ndim)+m%grid(igrid(ind))%f(icell(ind),1:ndim)*vol(ind)
#endif
        end do
     endif

     ! Perform kick, or drift, or both
     if(action_part==action_kick_drift)then

        ! Update velocity
        p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)
        
        ! Update position
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)

     else if(action_part.EQ.action_kick_only)then

        ! Compute proper time step for second kick
        if (p%levelp(ipart)>=ilevel)then
           dteff=g%dtnew(p%levelp(ipart))
        else
           dteff=g%dtold(p%levelp(ipart))
        endif

        ! Update level
        p%levelp(ipart)=ilevel

        ! Update velocity
        p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*dteff
        
     endif

  end do
  ! End loop over particles
  
  call close_cache_2(r,g,m,m%grid_dict)

  ! Periodic boundary conditions
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(p%xp(ipart,idim)>r%boxlen)then
              p%xp(ipart,idim)=p%xp(ipart,idim)-r%boxlen
           end if
           if(p%xp(ipart,idim)<0.d0)then
              p%xp(ipart,idim)=p%xp(ipart,idim)+r%boxlen
           end if
        end do
     end do
  end if

111 format('   Entering kick_and_drift_part for level',i2)

end subroutine kick_drift_part_2
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
