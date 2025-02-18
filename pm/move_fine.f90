module move_fine_module
contains
!################################################################
!################################################################
!################################################################
!################################################################
subroutine m_kick_drift_part(pst,ilevel,action_part)
  use amr_parameters, only: ndim,dp,twotondim
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  integer::ilevel
  integer::action_part
  !--------------------------------------------------------------
  ! Move particles according to kick-drift leap frog scheme.
  !--------------------------------------------------------------
  integer,dimension(1:2)::input_array
  integer::dummy(2)

  if(pst%s%m%noct_tot(ilevel)==0)return
  if(pst%s%r%verbose)write(*,'("   Entering kick_drift_part for level",i2," and action=",i2)')ilevel,action_part

  input_array(1)=ilevel
  input_array(2)=action_part
  call r_kick_drift_part(pst,input_array,2,dummy,0)

end subroutine m_kick_drift_part
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_kick_drift_part(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use amr_parameters, only: dp
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel
  integer::action_part
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_KICK_DRIFT_PART,pst%iUpper+1,input_size,output_size,input_array)
     call r_kick_drift_part(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     ilevel=input_array(1)
     action_part=input_array(2)
                     call kick_drift_part(pst%s,pst%s%p   ,ilevel,action_part)
     if(pst%s%r%star)call kick_drift_part(pst%s,pst%s%star,ilevel,action_part)
     if(pst%s%r%sink)call kick_drift_part(pst%s,pst%s%sink,ilevel,action_part)
     if(pst%s%r%tree)call kick_drift_part(pst%s,pst%s%tree,ilevel,action_part)
  endif

end subroutine r_kick_drift_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kick_drift_part(s,p,ilevel,action_part)
  use amr_parameters, only: dp,ndim,twotondim
  use pm_parameters
  use pm_commons, only: part_t
  use amr_commons, only: nbor
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  !
  !
  real(dp),dimension(1:ndim)::x,dr,dl
  integer,dimension(1:ndim)::ir,il
  real(dp),dimension(1:twotondim)::vol
  integer,dimension(1:ndim,1:twotondim)::ckey
  integer,dimension(1:twotondim)::icell
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,ind,idim
  real(kind=8)::dx_loc,vol_loc,dteff
  real(kind=8)::gamma,norm2,fnorm,delta
  real(dp),dimension(1:ndim)::ff
  logical::ok_level
  type(nbor),dimension(1:twotondim)::gridp
  type(msg_three_realdp)::dummy_three_realdp
  
  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Open read-only cache
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_three_realdp)/32,&
                     pack=pack_fetch_kick,unpack=unpack_fetch_kick)

  ! Loop over particles
  do ipart=p%headp(ilevel),p%tailp(ilevel)
     
     ! Rescale particle position at level ilevel
     do idim=1,ndim
        x(idim)=p%xp(ipart,idim)/dx_loc
     end do
     
     ! CIC at level ilevel (dr: right cloud boundary; dl: left cloud boundary)
     do idim=1,ndim
        dr(idim)=x(idim)+0.5D0
        ir(idim)=int(dr(idim))
        dr(idim)=dr(idim)-ir(idim)
        dl(idim)=1.0D0-dr(idim)
        il(idim)=ir(idim)-1
     end do
     
     ! Periodic boundary conditions
     do idim=1,ndim
        if(il(idim)<0)il(idim)=m%ckey_max(ilevel+1)-1
        if(ir(idim)==m%ckey_max(ilevel+1))ir(idim)=0
     enddo
     
     ! Compute cells Cartesian key
#if NDIM==1
     ckey(1,1)=il(1)
     ckey(1,2)=ir(1)
#endif
#if NDIM==2
     ckey(1:2,1)=(/il(1),il(2)/)
     ckey(1:2,2)=(/ir(1),il(2)/)
     ckey(1:2,3)=(/il(1),ir(2)/)
     ckey(1:2,4)=(/ir(1),ir(2)/)
#endif
#if NDIM==3
     ckey(1:3,1)=(/il(1),il(2),il(3)/)
     ckey(1:3,2)=(/ir(1),il(2),il(3)/)
     ckey(1:3,3)=(/il(1),ir(2),il(3)/)
     ckey(1:3,4)=(/ir(1),ir(2),il(3)/)
     ckey(1:3,5)=(/il(1),il(2),ir(3)/)
     ckey(1:3,6)=(/ir(1),il(2),ir(3)/)
     ckey(1:3,7)=(/il(1),ir(2),ir(3)/)
     ckey(1:3,8)=(/ir(1),ir(2),ir(3)/)
#endif
     
     ! Get parent cell at level ilevel using read-only cache
     ok_level=.true.
     hash_nbor(0)=ilevel+1
     icell=0
     do ind=1,twotondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp(ind)%p,icell(ind),flush_cache=.false.,fetch_cache=.true.,lock=.true.)
        if(.not.associated(gridp(ind)%p))then
           ok_level=.false.
        end if
     end do
     do ind=1,twotondim
        call unlock_cache(s,gridp(ind)%p)
     end do

     ! If cloud is not fully inside level ilevel, re-do CIC at coarser level
     if(.not. ok_level)then

        ! Rescale particle position at level ilevel
        do idim=1,ndim
           x(idim)=x(idim)/2.0d0
        end do
        
        ! CIC at level ilevel-1 (dr: right cloud boundary; dl: left cloud boundary)
        do idim=1,ndim
           dr(idim)=x(idim)+0.5D0
           ir(idim)=int(dr(idim))
           dr(idim)=dr(idim)-ir(idim)
           dl(idim)=1.0D0-dr(idim)
           il(idim)=ir(idim)-1
        end do
        
        ! Periodic boundary conditions
        do idim=1,ndim
           if(il(idim)<0)il(idim)=m%ckey_max(ilevel)-1
           if(ir(idim)==m%ckey_max(ilevel))ir(idim)=0
        enddo
        
        ! Compute cells Cartesian key
#if NDIM==1
        ckey(1,1)=il(1)
        ckey(1,2)=ir(1)
#endif
#if NDIM==2
        ckey(1:2,1)=(/il(1),il(2)/)
        ckey(1:2,2)=(/ir(1),il(2)/)
        ckey(1:2,3)=(/il(1),ir(2)/)
        ckey(1:2,4)=(/ir(1),ir(2)/)
#endif
#if NDIM==3
        ckey(1:3,1)=(/il(1),il(2),il(3)/)
        ckey(1:3,2)=(/ir(1),il(2),il(3)/)
        ckey(1:3,3)=(/il(1),ir(2),il(3)/)
        ckey(1:3,4)=(/ir(1),ir(2),il(3)/)
        ckey(1:3,5)=(/il(1),il(2),ir(3)/)
        ckey(1:3,6)=(/ir(1),il(2),ir(3)/)
        ckey(1:3,7)=(/il(1),ir(2),ir(3)/)
        ckey(1:3,8)=(/ir(1),ir(2),ir(3)/)
#endif
        
        ! Get parent cell at level ilevel-1 using read-only cache
        ok_level=.true.
        hash_nbor(0)=ilevel
        icell=0
        do ind=1,twotondim
           hash_nbor(1:ndim)=ckey(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp(ind)%p,icell(ind),flush_cache=.false.,fetch_cache=.true.,lock=.true.)
           if(.not.associated(gridp(ind)%p))then
              ok_level=.false.
           end if
        end do
        do ind=1,twotondim
           call unlock_cache(s,gridp(ind)%p)
        end do
     end if
        
     ! Compute cloud volumes
#if NDIM==1
     vol(1)=dl(1)
     vol(2)=dr(1)
#endif
#if NDIM==2
     vol(1)=dl(1)*dl(2)
     vol(2)=dr(1)*dl(2)
     vol(3)=dl(1)*dr(2)
     vol(4)=dr(1)*dr(2)
#endif
#if NDIM==3
     vol(1)=dl(1)*dl(2)*dl(3)
     vol(2)=dr(1)*dl(2)*dl(3)
     vol(3)=dl(1)*dr(2)*dl(3)
     vol(4)=dr(1)*dr(2)*dl(3)
     vol(5)=dl(1)*dl(2)*dr(3)
     vol(6)=dr(1)*dl(2)*dr(3)
     vol(7)=dl(1)*dr(2)*dr(3)
     vol(8)=dr(1)*dr(2)*dr(3)
#endif
     
     ! Gather 3-force
     ff(1:ndim)=0.0
     if(ok_level)then
        do ind=1,twotondim
#ifdef GRAV
           ff(1:ndim)=ff(1:ndim)+gridp(ind)%p%f(icell(ind),1:ndim)*vol(ind)
#endif
        end do
     endif

     ! Perform kick, or drift, or both
     if(action_part==action_kick_drift)then

        ! Update velocity
        p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)
        
        ! Update position
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)

        ! For sink particle only
        if(p%type==SINK_TYPE.and.r%sink_descent)then

           ! Compute gradient descent coefficients
           gamma = 0d0
           norm2 = 0d0
           fnorm = 0d0
           do idim=1,ndim
              gamma = gamma + p%vp(ipart,idim)*(ff(idim)-p%fp(ipart,idim))
              norm2 = norm2 + (ff(idim)-p%fp(ipart,idim))**2
              fnorm = fnorm + ff(idim)**2
           enddo
           gamma = gamma*g%dtnew(ilevel) ! cm2 s-2
           fnorm = sqrt(fnorm)
           delta = 0d0
           if(norm2>0)then
              delta = MIN(r%fudge_descent*g%dtnew(ilevel)*sqrt(abs(gamma)/norm2)*fnorm,0.5d0*dx_loc)
           endif
!!$           delta = MIN(r%fudge_descent*g%dtnew(ilevel)*sqrt(abs(gamma)),0.5d0*dx_loc)

           ! Update particle positions
           if(fnorm>0)then
              p%xp(ipart,1:ndim) = p%xp(ipart,1:ndim) + ff(1:ndim)/fnorm*delta
           endif
!!$           if(norm2>0)then
!!$              p%xp(ipart,1:ndim) = p%xp(ipart,1:ndim) + (ff(1:ndim)-p%fp(ipart,1:ndim))/sqrt(norm2)*delta
!!$           endif

           ! Store old force
           p%fp(ipart,1:ndim)=ff(1:ndim)

        endif

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
  
  call close_cache(s,m%grid_dict)

  ! Periodic boundary conditions
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(p%xp(ipart,idim)<   0.0d0 )p%xp(ipart,idim)=p%xp(ipart,idim)+r%boxlen
           if(p%xp(ipart,idim)>=r%boxlen)p%xp(ipart,idim)=p%xp(ipart,idim)-r%boxlen
        end do
     end do
  end if

  end associate
  
end subroutine kick_drift_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine tsc_kick_drift_part(s,p,ilevel,action_part)
  use amr_parameters, only: dp, ndim, threetondim
  use pm_parameters
  use pm_commons, only: part_t
  use oct_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  !
  !
  real(dp),dimension(1:ndim)::x,wl,wc,wr
  integer,dimension(1:ndim)::cl,cc,cr
  real(dp),dimension(1:threetondim)::vol
  integer,dimension(1:ndim,1:threetondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,icell,ind,idim
  real(kind=8)::xl,xc,xr
  real(kind=8)::dx_loc,vol_loc,dteff
  real(kind=8)::gamma,norm2,fnorm,delta
  real(dp),dimension(1:ndim)::ff
  logical::ok_level
  type(oct),pointer::gridp
  type(msg_three_realdp)::dummy_three_realdp

  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel
  vol_loc=dx_loc**ndim

  ! Open read-only cache
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_three_realdp)/32,&
                     pack=pack_fetch_kick,unpack=unpack_fetch_kick)

  ! Loop over particles
  do ipart=p%headp(ilevel),p%tailp(ilevel)

     ! Rescale particle position at level ilevel
     do idim=1,ndim
        x(idim)=p%xp(ipart,idim)/dx_loc
     end do

     ! TSC at level ilevel; a particle contributes to 3 cells in each direction
     do idim=1,ndim
        cl(idim)=int(x(idim))-1 ! cell index
        cc(idim)=int(x(idim))
        cr(idim)=int(x(idim))+1
        xl=dble(cl(idim))+0.5D0 ! cell coordinate
        xc=dble(cc(idim))+0.5D0
        xr=dble(cr(idim))+0.5D0
        wl(idim)=0.5D0*(1.5D0-abs(x(idim)-xl))**2 ! weight
        wc(idim)=0.75D0-         (x(idim)-xc) **2
        wr(idim)=0.5D0*(1.5D0-abs(x(idim)-xr))**2
     end do

     ! Periodic boundary conditions
     do idim=1,ndim
        if(cl(idim)<0)cl(idim)=m%ckey_max(ilevel+1)-1
        if(cr(idim)==m%ckey_max(ilevel+1))cr(idim)=0
     enddo

     ! Compute cells Cartesian key
#if NDIM==1
     ckey(1,1)=cl(1)
     ckey(1,2)=cc(1)
     ckey(1,3)=cr(1)
#endif
#if NDIM==2
     ckey(1:2,1)=(/cl(1),cl(2)/)
     ckey(1:2,2)=(/cc(1),cl(2)/)
     ckey(1:2,3)=(/cr(1),cl(2)/)
     ckey(1:2,4)=(/cl(1),cc(2)/)
     ckey(1:2,5)=(/cc(1),cc(2)/)
     ckey(1:2,6)=(/cr(1),cc(2)/)
     ckey(1:2,7)=(/cl(1),cr(2)/)
     ckey(1:2,8)=(/cc(1),cr(2)/)
     ckey(1:2,9)=(/cr(1),cr(2)/)
#endif
#if NDIM==3
     ckey(1:3,1) =(/cl(1),cl(2),cl(3)/)
     ckey(1:3,2) =(/cc(1),cl(2),cl(3)/)
     ckey(1:3,3) =(/cr(1),cl(2),cl(3)/)
     ckey(1:3,4) =(/cl(1),cc(2),cl(3)/)
     ckey(1:3,5) =(/cc(1),cc(2),cl(3)/)
     ckey(1:3,6) =(/cr(1),cc(2),cl(3)/)
     ckey(1:3,7) =(/cl(1),cr(2),cl(3)/)
     ckey(1:3,8) =(/cc(1),cr(2),cl(3)/)
     ckey(1:3,9) =(/cr(1),cr(2),cl(3)/)
     ckey(1:3,10)=(/cl(1),cl(2),cc(3)/)
     ckey(1:3,11)=(/cc(1),cl(2),cc(3)/)
     ckey(1:3,12)=(/cr(1),cl(2),cc(3)/)
     ckey(1:3,13)=(/cl(1),cc(2),cc(3)/)
     ckey(1:3,14)=(/cc(1),cc(2),cc(3)/)
     ckey(1:3,15)=(/cr(1),cc(2),cc(3)/)
     ckey(1:3,16)=(/cl(1),cr(2),cc(3)/)
     ckey(1:3,17)=(/cc(1),cr(2),cc(3)/)
     ckey(1:3,18)=(/cr(1),cr(2),cc(3)/)
     ckey(1:3,19)=(/cl(1),cl(2),cr(3)/)
     ckey(1:3,20)=(/cc(1),cl(2),cr(3)/)
     ckey(1:3,21)=(/cr(1),cl(2),cr(3)/)
     ckey(1:3,22)=(/cl(1),cc(2),cr(3)/)
     ckey(1:3,23)=(/cc(1),cc(2),cr(3)/)
     ckey(1:3,24)=(/cr(1),cc(2),cr(3)/)
     ckey(1:3,25)=(/cl(1),cr(2),cr(3)/)
     ckey(1:3,26)=(/cc(1),cr(2),cr(3)/)
     ckey(1:3,27)=(/cr(1),cr(2),cr(3)/)
#endif

     ! Compute cloud volumes
#if NDIM==1
     vol(1)=wl(1)
     vol(2)=wc(1)
     vol(3)=wr(1)
#endif
#if NDIM==2
     vol(1)=wl(1)*wl(2)
     vol(2)=wc(1)*wl(2)
     vol(3)=wr(1)*wl(2)
     vol(4)=wl(1)*wc(2)
     vol(5)=wc(1)*wc(2)
     vol(6)=wr(1)*wc(2)
     vol(7)=wl(1)*wr(2)
     vol(8)=wc(1)*wr(2)
     vol(9)=wr(1)*wr(2)
#endif
#if NDIM==3
     vol(1) =wl(1)*wl(2)*wl(3)
     vol(2) =wc(1)*wl(2)*wl(3)
     vol(3) =wr(1)*wl(2)*wl(3)
     vol(4) =wl(1)*wc(2)*wl(3)
     vol(5) =wc(1)*wc(2)*wl(3)
     vol(6) =wr(1)*wc(2)*wl(3)
     vol(7) =wl(1)*wr(2)*wl(3)
     vol(8) =wc(1)*wr(2)*wl(3)
     vol(9) =wr(1)*wr(2)*wl(3)
     vol(10)=wl(1)*wl(2)*wc(3)
     vol(11)=wc(1)*wl(2)*wc(3)
     vol(12)=wr(1)*wl(2)*wc(3)
     vol(13)=wl(1)*wc(2)*wc(3)
     vol(14)=wc(1)*wc(2)*wc(3)
     vol(15)=wr(1)*wc(2)*wc(3)
     vol(16)=wl(1)*wr(2)*wc(3)
     vol(17)=wc(1)*wr(2)*wc(3)
     vol(18)=wr(1)*wr(2)*wc(3)
     vol(19)=wl(1)*wl(2)*wr(3)
     vol(20)=wc(1)*wl(2)*wr(3)
     vol(21)=wr(1)*wl(2)*wr(3)
     vol(22)=wl(1)*wc(2)*wr(3)
     vol(23)=wc(1)*wc(2)*wr(3)
     vol(24)=wr(1)*wc(2)*wr(3)
     vol(25)=wl(1)*wr(2)*wr(3)
     vol(26)=wc(1)*wr(2)*wr(3)
     vol(27)=wr(1)*wr(2)*wr(3)
#endif

     ! Gather 3-force
     hash_nbor(0)=ilevel+1
     ff(1:ndim)=0.0
     do ind=1,threetondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        ! Get parent cell at level ilevel using read-only cache
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef GRAV
        if(associated(gridp))then
           ff(1:ndim)=ff(1:ndim)+gridp%f(icell,1:ndim)*vol(ind)
        end if
#endif
     end do

     ! Perform kick, or drift, or both
     if(action_part==action_kick_drift)then

        ! Update velocity
        p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)

        ! Update position
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)

        ! For sink particle only
        if(p%type==SINK_TYPE.and.r%sink_descent)then

           ! Compute gradient descent coefficients
           gamma = 0d0
           norm2 = 0d0
           fnorm = 0d0
           do idim=1,ndim
              gamma = gamma + p%vp(ipart,idim)*(ff(idim)-p%fp(ipart,idim))
              norm2 = norm2 + (ff(idim)-p%fp(ipart,idim))**2
              fnorm = fnorm + ff(idim)**2
           enddo
           gamma = gamma*g%dtnew(ilevel) ! cm2 s-2
           fnorm = sqrt(fnorm)
           delta = 0d0
           if(norm2>0)then
              delta = MIN(r%fudge_descent*g%dtnew(ilevel)*sqrt(abs(gamma)/norm2)*fnorm,0.5d0*dx_loc)
           endif
!!$           delta = MIN(r%fudge_descent*g%dtnew(ilevel)*sqrt(abs(gamma)),0.5d0*dx_loc)

           ! Update particle positions
           if(fnorm>0)then
              p%xp(ipart,1:ndim) = p%xp(ipart,1:ndim) + ff(1:ndim)/fnorm*delta
           endif
!!$           if(norm2>0)then
!!$              p%xp(ipart,1:ndim) = p%xp(ipart,1:ndim) + (ff(1:ndim)-p%fp(ipart,1:ndim))/sqrt(norm2)*delta
!!$           endif

           ! Store old force
           p%fp(ipart,1:ndim)=ff(1:ndim)

        endif

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

  call close_cache(s,m%grid_dict)

  ! Periodic boundary conditions
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(p%xp(ipart,idim)<   0.0d0 )p%xp(ipart,idim)=p%xp(ipart,idim)+r%boxlen
           if(p%xp(ipart,idim)>=r%boxlen)p%xp(ipart,idim)=p%xp(ipart,idim)-r%boxlen
        end do
     end do
  end if

  end associate

end subroutine tsc_kick_drift_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pcs_kick_drift_part(s,p,ilevel,action_part)
  use amr_parameters, only: dp, ndim, fourtondim
  use pm_parameters
  use pm_commons, only: part_t
  use oct_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  !
  !
  real(dp),dimension(1:ndim)::x,wll,wl,wr,wrr
  integer,dimension(1:ndim)::cll,cl,cr,crr
  real(dp),dimension(1:fourtondim)::vol
  integer,dimension(1:ndim,1:fourtondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,icell,ind,idim
  real(kind=8)::xll,xl,xr,xrr
  real(kind=8)::dx_loc,vol_loc,dteff
  real(kind=8)::gamma,norm2,fnorm,delta
  real(dp),dimension(1:ndim)::ff
  logical::ok_level
  type(oct),pointer::gridp
  type(msg_three_realdp)::dummy_three_realdp

  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel
  vol_loc=dx_loc**ndim

  ! Open read-only cache
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
       hilbert=m%domain, pack_size=storage_size(dummy_three_realdp)/32,&
       pack=pack_fetch_kick,unpack=unpack_fetch_kick)

  ! Loop over particles
  do ipart=p%headp(ilevel),p%tailp(ilevel)

     ! Rescale particle position at level ilevel
     do idim=1,ndim
        x(idim)=p%xp(ipart,idim)/dx_loc
     end do

     ! PCS at level ilevel; a particle contributes to 4 cells in each direction
     do idim=1,ndim
        crr(idim)=int(x(idim)+1.5D0) ! rightermost cell index
        cr (idim)=crr(idim)-1
        cl (idim)=crr(idim)-2
        cll(idim)=crr(idim)-3
        xll=dble(cll(idim))+0.5D0 ! cell coordinate
        xl =dble(cl (idim))+0.5D0
        xr =dble(cr (idim))+0.5D0
	xrr=dble(crr(idim))+0.5D0
        wll(idim)=(2D0                        -abs(x(idim)-xll)**3)/6D0 ! weight
        wl (idim)=(4D0-6D0*(x(idim)-xl)**2+3d0*abs(x(idim)-xl )**3)/6D0
        wr (idim)=(4D0-6D0*(x(idim)-xr)**2+3d0*abs(x(idim)-xr )**3)/6D0
        wrr(idim)=(2D0                        -abs(x(idim)-xrr)**3)/6D0
     end do

     ! Periodic boundary conditions
     do idim=1,ndim
        if(cll(idim)<0)cll(idim)=m%ckey_max(ilevel+1)-1
        if(cl (idim)<0)cl (idim)=m%ckey_max(ilevel+1)-1
        if(cr (idim)==m%ckey_max(ilevel+1))cr (idim)=0
        if(crr(idim)==m%ckey_max(ilevel+1))crr(idim)=0
     enddo

     ! Compute cells Cartesian key
#if NDIM==1
     ckey(1,1)=cll(1)
     ckey(1,2)=cl (1)
     ckey(1,3)=cr (1)
     ckey(1,4)=crr(1)
#endif
#if NDIM==2
     ckey(1:2,1) =(/cll(1),cll(2)/)
     ckey(1:2,2) =(/cl (1),cll(2)/)
     ckey(1:2,3) =(/cr (1),cll(2)/)
     ckey(1:2,4) =(/crr(1),cll(2)/)
     ckey(1:2,5) =(/cll(1),cl (2)/)
     ckey(1:2,6) =(/cl (1),cl (2)/)
     ckey(1:2,7) =(/cr (1),cl (2)/)
     ckey(1:2,8) =(/crr(1),cl (2)/)
     ckey(1:2,9) =(/cll(1),cr (2)/)
     ckey(1:2,10)=(/cl (1),cr (2)/)
     ckey(1:2,11)=(/cr (1),cr (2)/)
     ckey(1:2,12)=(/crr(1),cr (2)/)
     ckey(1:2,13)=(/cll(1),crr(2)/)
     ckey(1:2,14)=(/cl (1),crr(2)/)
     ckey(1:2,15)=(/cr (1),crr(2)/)
     ckey(1:2,16)=(/crr(1),crr(2)/)
#endif
#if NDIM==3
     ckey(1:3,1) =(/cll(1),cll(2),cll(3)/)
     ckey(1:3,2) =(/cl (1),cll(2),cll(3)/)
     ckey(1:3,3) =(/cr (1),cll(2),cll(3)/)
     ckey(1:3,4) =(/crr(1),cll(2),cll(3)/)
     ckey(1:3,5) =(/cll(1),cl (2),cll(3)/)
     ckey(1:3,6) =(/cl (1),cl (2),cll(3)/)
     ckey(1:3,7) =(/cr (1),cl (2),cll(3)/)
     ckey(1:3,8) =(/crr(1),cl (2),cll(3)/)
     ckey(1:3,9) =(/cll(1),cr (2),cll(3)/)
     ckey(1:3,10)=(/cl (1),cr (2),cll(3)/)
     ckey(1:3,11)=(/cr (1),cr (2),cll(3)/)
     ckey(1:3,12)=(/crr(1),cr (2),cll(3)/)
     ckey(1:3,13)=(/cll(1),crr(2),cll(3)/)
     ckey(1:3,14)=(/cl (1),crr(2),cll(3)/)
     ckey(1:3,15)=(/cr (1),crr(2),cll(3)/)
     ckey(1:3,16)=(/crr(1),crr(2),cll(3)/)
     ckey(1:3,17)=(/cll(1),cll(2),cl (3)/)
     ckey(1:3,18)=(/cl (1),cll(2),cl (3)/)
     ckey(1:3,19)=(/cr (1),cll(2),cl (3)/)
     ckey(1:3,20)=(/crr(1),cll(2),cl (3)/)
     ckey(1:3,21)=(/cll(1),cl (2),cl (3)/)
     ckey(1:3,22)=(/cl (1),cl (2),cl (3)/)
     ckey(1:3,23)=(/cr (1),cl (2),cl (3)/)
     ckey(1:3,24)=(/crr(1),cl (2),cl (3)/)
     ckey(1:3,25)=(/cll(1),cr (2),cl (3)/)
     ckey(1:3,26)=(/cl (1),cr (2),cl (3)/)
     ckey(1:3,27)=(/cr (1),cr (2),cl (3)/)
     ckey(1:3,28)=(/crr(1),cr (2),cl (3)/)
     ckey(1:3,29)=(/cll(1),crr(2),cl (3)/)
     ckey(1:3,30)=(/cl (1),crr(2),cl (3)/)
     ckey(1:3,31)=(/cr (1),crr(2),cl (3)/)
     ckey(1:3,32)=(/crr(1),crr(2),cl (3)/)
     ckey(1:3,33)=(/cll(1),cll(2),cr (3)/)
     ckey(1:3,34)=(/cl (1),cll(2),cr (3)/)
     ckey(1:3,35)=(/cr (1),cll(2),cr (3)/)
     ckey(1:3,36)=(/crr(1),cll(2),cr (3)/)
     ckey(1:3,37)=(/cll(1),cl (2),cr (3)/)
     ckey(1:3,38)=(/cl (1),cl (2),cr (3)/)
     ckey(1:3,39)=(/cr (1),cl (2),cr (3)/)
     ckey(1:3,40)=(/crr(1),cl (2),cr (3)/)
     ckey(1:3,41)=(/cll(1),cr (2),cr (3)/)
     ckey(1:3,42)=(/cl (1),cr (2),cr (3)/)
     ckey(1:3,43)=(/cr (1),cr (2),cr (3)/)
     ckey(1:3,44)=(/crr(1),cr (2),cr (3)/)
     ckey(1:3,45)=(/cll(1),crr(2),cr (3)/)
     ckey(1:3,46)=(/cl (1),crr(2),cr (3)/)
     ckey(1:3,47)=(/cr (1),crr(2),cr (3)/)
     ckey(1:3,48)=(/crr(1),crr(2),cr (3)/)
     ckey(1:3,49)=(/cll(1),cll(2),crr(3)/)
     ckey(1:3,50)=(/cl (1),cll(2),crr(3)/)
     ckey(1:3,51)=(/cr (1),cll(2),crr(3)/)
     ckey(1:3,52)=(/crr(1),cll(2),crr(3)/)
     ckey(1:3,53)=(/cll(1),cl (2),crr(3)/)
     ckey(1:3,54)=(/cl (1),cl (2),crr(3)/)
     ckey(1:3,55)=(/cr (1),cl (2),crr(3)/)
     ckey(1:3,56)=(/crr(1),cl (2),crr(3)/)
     ckey(1:3,57)=(/cll(1),cr (2),crr(3)/)
     ckey(1:3,58)=(/cl (1),cr (2),crr(3)/)
     ckey(1:3,59)=(/cr (1),cr (2),crr(3)/)
     ckey(1:3,60)=(/crr(1),cr (2),crr(3)/)
     ckey(1:3,61)=(/cll(1),crr(2),crr(3)/)
     ckey(1:3,62)=(/cl (1),crr(2),crr(3)/)
     ckey(1:3,63)=(/cr (1),crr(2),crr(3)/)
     ckey(1:3,64)=(/crr(1),crr(2),crr(3)/)
#endif

     ! Compute cloud volumes
#if NDIM==1
     vol(1)=wll(1)
     vol(2)=wl (1)
     vol(3)=wr (1)
     vol(4)=wrr(1)
#endif
#if NDIM==2
     vol(1) =wll(1)*wll(2)
     vol(2) =wl (1)*wll(2)
     vol(3) =wr (1)*wll(2)
     vol(4) =wrr(1)*wll(2)
     vol(5) =wll(1)*wl (2)
     vol(6) =wl (1)*wl (2)
     vol(7) =wr (1)*wl (2)
     vol(8) =wrr(1)*wl (2)
     vol(9) =wll(1)*wr (2)
     vol(10)=wl (1)*wr (2)
     vol(11)=wr (1)*wr (2)
     vol(12)=wrr(1)*wr (2)
     vol(13)=wll(1)*wrr(2)
     vol(14)=wl (1)*wrr(2)
     vol(15)=wr (1)*wrr(2)
     vol(16)=wrr(1)*wrr(2)
#endif
#if NDIM==3
     vol(1) =wll(1)*wll(2)*wll(3)
     vol(2) =wl (1)*wll(2)*wll(3)
     vol(3) =wr (1)*wll(2)*wll(3)
     vol(4) =wrr(1)*wll(2)*wll(3)
     vol(5) =wll(1)*wl (2)*wll(3)
     vol(6) =wl (1)*wl (2)*wll(3)
     vol(7) =wr (1)*wl (2)*wll(3)
     vol(8) =wrr(1)*wl (2)*wll(3)
     vol(9) =wll(1)*wr (2)*wll(3)
     vol(10)=wl (1)*wr (2)*wll(3)
     vol(11)=wr (1)*wr (2)*wll(3)
     vol(12)=wrr(1)*wr (2)*wll(3)
     vol(13)=wll(1)*wrr(2)*wll(3)
     vol(14)=wl (1)*wrr(2)*wll(3)
     vol(15)=wr (1)*wrr(2)*wll(3)
     vol(16)=wrr(1)*wrr(2)*wll(3)
     vol(17)=wll(1)*wll(2)*wl (3)
     vol(18)=wl (1)*wll(2)*wl (3)
     vol(19)=wr (1)*wll(2)*wl (3)
     vol(20)=wrr(1)*wll(2)*wl (3)
     vol(21)=wll(1)*wl (2)*wl (3)
     vol(22)=wl (1)*wl (2)*wl (3)
     vol(23)=wr (1)*wl (2)*wl (3)
     vol(24)=wrr(1)*wl (2)*wl (3)
     vol(25)=wll(1)*wr (2)*wl (3)
     vol(26)=wl (1)*wr (2)*wl (3)
     vol(27)=wr (1)*wr (2)*wl (3)
     vol(28)=wrr(1)*wr (2)*wl (3)
     vol(29)=wll(1)*wrr(2)*wl (3)
     vol(30)=wl (1)*wrr(2)*wl (3)
     vol(31)=wr (1)*wrr(2)*wl (3)
     vol(32)=wrr(1)*wrr(2)*wl (3)
     vol(33)=wll(1)*wll(2)*wr (3)
     vol(34)=wl (1)*wll(2)*wr (3)
     vol(35)=wr (1)*wll(2)*wr (3)
     vol(36)=wrr(1)*wll(2)*wr (3)
     vol(37)=wll(1)*wl (2)*wr (3)
     vol(38)=wl (1)*wl (2)*wr (3)
     vol(39)=wr (1)*wl (2)*wr (3)
     vol(40)=wrr(1)*wl (2)*wr (3)
     vol(41)=wll(1)*wr (2)*wr (3)
     vol(42)=wl (1)*wr (2)*wr (3)
     vol(43)=wr (1)*wr (2)*wr (3)
     vol(44)=wrr(1)*wr (2)*wr (3)
     vol(45)=wll(1)*wrr(2)*wr (3)
     vol(46)=wl (1)*wrr(2)*wr (3)
     vol(47)=wr (1)*wrr(2)*wr (3)
     vol(48)=wrr(1)*wrr(2)*wr (3)
     vol(49)=wll(1)*wll(2)*wrr(3)
     vol(50)=wl (1)*wll(2)*wrr(3)
     vol(51)=wr (1)*wll(2)*wrr(3)
     vol(52)=wrr(1)*wll(2)*wrr(3)
     vol(53)=wll(1)*wl (2)*wrr(3)
     vol(54)=wl (1)*wl (2)*wrr(3)
     vol(55)=wr (1)*wl (2)*wrr(3)
     vol(56)=wrr(1)*wl (2)*wrr(3)
     vol(57)=wll(1)*wr (2)*wrr(3)
     vol(58)=wl (1)*wr (2)*wrr(3)
     vol(59)=wr (1)*wr (2)*wrr(3)
     vol(60)=wrr(1)*wr (2)*wrr(3)
     vol(61)=wll(1)*wrr(2)*wrr(3)
     vol(62)=wl (1)*wrr(2)*wrr(3)
     vol(63)=wr (1)*wrr(2)*wrr(3)
     vol(64)=wrr(1)*wrr(2)*wrr(3)
#endif

     ! Gather 3-force
     hash_nbor(0)=ilevel+1
     ff(1:ndim)=0.0
     do ind=1,fourtondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        ! Get parent cell at level ilevel using read-only cache
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef GRAV
        if(associated(gridp))then
           ff(1:ndim)=ff(1:ndim)+gridp%f(icell,1:ndim)*vol(ind)
        end if
#endif
     end do

     ! Perform kick, or drift, or both
     if(action_part==action_kick_drift)then

        ! Update velocity
        p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)

        ! Update position
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)

        ! For sink particle only
        if(p%type==SINK_TYPE.and.r%sink_descent)then

           ! Compute gradient descent coefficients
           gamma = 0d0
           norm2 = 0d0
           fnorm = 0d0
           do idim=1,ndim
              gamma = gamma + p%vp(ipart,idim)*(ff(idim)-p%fp(ipart,idim))
              norm2 = norm2 + (ff(idim)-p%fp(ipart,idim))**2
              fnorm = fnorm + ff(idim)**2
           enddo
           gamma = gamma*g%dtnew(ilevel) ! cm2 s-2
           fnorm = sqrt(fnorm)
           delta = 0d0
           if(norm2>0)then
              delta = MIN(r%fudge_descent*g%dtnew(ilevel)*sqrt(abs(gamma)/norm2)*fnorm,0.5d0*dx_loc)
           endif
!!$           delta = MIN(r%fudge_descent*g%dtnew(ilevel)*sqrt(abs(gamma)),0.5d0*dx_loc)

           ! Update particle positions
           if(fnorm>0)then
              p%xp(ipart,1:ndim) = p%xp(ipart,1:ndim) + ff(1:ndim)/fnorm*delta
           endif
!!$           if(norm2>0)then
!!$              p%xp(ipart,1:ndim) = p%xp(ipart,1:ndim) + (ff(1:ndim)-p%fp(ipart,1:ndim))/sqrt(norm2)*delta
!!$           endif

           ! Store old force
           p%fp(ipart,1:ndim)=ff(1:ndim)

        endif

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

  call close_cache(s,m%grid_dict)

  ! Periodic boundary conditions
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(p%xp(ipart,idim)<   0.0d0 )p%xp(ipart,idim)=p%xp(ipart,idim)+r%boxlen
           if(p%xp(ipart,idim)>=r%boxlen)p%xp(ipart,idim)=p%xp(ipart,idim)-r%boxlen
        end do
     end do
  end if

  end associate

end subroutine pcs_kick_drift_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine pack_fetch_kick(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_three_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_three_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp_phi(ind)=grid%f(ind,1)
     msg%realdp_phi_old(ind)=grid%f(ind,2)
     msg%realdp_dis(ind)=grid%f(ind,3)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_kick
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine unpack_fetch_kick(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_three_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_three_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)
  
#ifdef GRAV
  do ind=1,twotondim
     grid%f(ind,1)=msg%realdp_phi(ind)
     grid%f(ind,2)=msg%realdp_phi_old(ind)
     grid%f(ind,3)=msg%realdp_dis(ind)
  end do
#endif

end subroutine unpack_fetch_kick
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module move_fine_module
