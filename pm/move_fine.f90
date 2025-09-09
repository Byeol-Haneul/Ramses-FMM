module move_fine_module
  use rho_fine_module, only: cic_weight, cic_index, tsc_weight, tsc_index, pcs_weight, pcs_index
contains
!################################################################
!################################################################
!################################################################
!################################################################
subroutine m_kick_drift_part(pst,ilevel,action_part)
  use amr_parameters, only: ndim, twotondim
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
     ! Force interpolation for various components (DM particles, star, sink, tree)
     ! based on their respective deposition schemes (CIC 1, TSC 2 or PCS 3)
     if(pst%s%r%part)then
        if(pst%s%r%part_force_interpolation_scheme==1)then
           call cic_kick_drift_part(pst%s,pst%s%p   ,ilevel,action_part)
        elseif(pst%s%r%part_force_interpolation_scheme==2)then
           call tsc_kick_drift_part(pst%s,pst%s%p   ,ilevel,action_part)
        elseif(pst%s%r%part_force_interpolation_scheme==3)then
           call pcs_kick_drift_part(pst%s,pst%s%p   ,ilevel,action_part)
        endif
     endif
     if(pst%s%r%star)then
        if(pst%s%r%star_force_interpolation_scheme==1)then
           call cic_kick_drift_part(pst%s,pst%s%star,ilevel,action_part)
        elseif(pst%s%r%star_force_interpolation_scheme==2)then
           call tsc_kick_drift_part(pst%s,pst%s%star,ilevel,action_part)
        elseif(pst%s%r%star_force_interpolation_scheme==3)then
           call pcs_kick_drift_part(pst%s,pst%s%star,ilevel,action_part)
        endif
     endif
     if(pst%s%r%sink)then
        if(pst%s%r%sink_force_interpolation_scheme==1)then
           call cic_kick_drift_part(pst%s,pst%s%sink,ilevel,action_part)
        elseif(pst%s%r%sink_force_interpolation_scheme==2)then
           call tsc_kick_drift_part(pst%s,pst%s%sink,ilevel,action_part)
        elseif(pst%s%r%sink_force_interpolation_scheme==3)then
           call pcs_kick_drift_part(pst%s,pst%s%sink,ilevel,action_part)
        endif
     endif
     if(pst%s%r%tree)then
        if(pst%s%r%tree_force_interpolation_scheme==1)then
           call cic_kick_drift_part(pst%s,pst%s%tree,ilevel,action_part)
        elseif(pst%s%r%tree_force_interpolation_scheme==2)then
           call tsc_kick_drift_part(pst%s,pst%s%tree,ilevel,action_part)
        elseif(pst%s%r%tree_force_interpolation_scheme==3)then
           call pcs_kick_drift_part(pst%s,pst%s%tree,ilevel,action_part)
        endif
     endif
     if(pst%s%r%trac)then
        if(pst%s%r%trac_force_interpolation_scheme==1)then
           call cic_kick_drift_trac(pst%s,pst%s%trac,ilevel,action_part)
        elseif(pst%s%r%trac_force_interpolation_scheme==2)then
           call tsc_kick_drift_trac(pst%s,pst%s%trac,ilevel,action_part)
        elseif(pst%s%r%trac_force_interpolation_scheme==3)then
           call pcs_kick_drift_trac(pst%s,pst%s%trac,ilevel,action_part)
        endif
     endif
  endif

end subroutine r_kick_drift_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine cic_kick_drift_part(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, twotondim
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
  real(kind=8),dimension(1:ndim)::x,dr,dl
  integer,dimension(1:ndim)::ir,il
  real(kind=8),dimension(1:twotondim)::vol
  integer,dimension(1:ndim,1:twotondim)::ckey
  integer,dimension(1:twotondim)::icell
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,ind,idim
  real(kind=8)::dx_loc,vol_loc,dteff
  real(kind=8)::gamma,norm2,fnorm,delta
  real(kind=8),dimension(1:ndim)::ff
  logical::ok_level
  type(nbor),dimension(1:twotondim)::gridp
  type(msg_three_realdp)::dummy_three_realdp

  associate(r=>s%r,g=>s%g,m=>s%m)

  if(p%static)then
     ! We still need to set the particle levels correctly, even if they are not moved
     ! Loop over particles
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        ! Update level
        p%levelp(ipart)=ilevel
     end do
     return
  end if

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
     ckey = cic_index(il,ir)

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
        ckey = cic_index(il,ir)

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
     vol = cic_weight(dl,dr)

     ! Gather 3-force
     ff(1:ndim)=0.0
     if(ok_level)then
        do ind=1,twotondim
#ifdef GRAV
           ff(1:ndim)=ff(1:ndim)+gridp(ind)%p%f(icell(ind),1:ndim)*vol(ind)
#else
           continue
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

end subroutine cic_kick_drift_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine tsc_kick_drift_part(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, threetondim
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
  real(kind=8),dimension(1:ndim)::x,wl,wc,wr
  integer,dimension(1:ndim)::cl,cc,cr
  real(kind=8),dimension(1:threetondim)::vol
  integer,dimension(1:ndim,1:threetondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,icell,ind,idim
  real(kind=8)::xl,xc,xr
  real(kind=8)::dx_loc,vol_loc,dteff
  real(kind=8)::gamma,norm2,fnorm,delta
  real(kind=8),dimension(1:ndim)::ff
  logical::ok_level
  type(oct),pointer::gridp
  type(msg_three_realdp)::dummy_three_realdp

  associate(r=>s%r,g=>s%g,m=>s%m)

  if(p%static)then
     ! We still need to set the particle levels correctly, even if they are not moved
     ! Loop over particles
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        ! Update level
        p%levelp(ipart)=ilevel
     end do
     return
  end if

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
     ckey = tsc_index(cl,cc,cr)

     ! Compute cloud volumes
     vol = tsc_weight(wl,wc,wr)

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
  use amr_parameters, only: ndim, fourtondim
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
  real(kind=8),dimension(1:ndim)::x,wll,wl,wr,wrr
  integer,dimension(1:ndim)::cll,cl,cr,crr
  real(kind=8),dimension(1:fourtondim)::vol
  integer,dimension(1:ndim,1:fourtondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,icell,ind,idim
  real(kind=8)::xll,xl,xr,xrr
  real(kind=8)::dx_loc,vol_loc,dteff
  real(kind=8)::gamma,norm2,fnorm,delta
  real(kind=8),dimension(1:ndim)::ff
  logical::ok_level
  type(oct),pointer::gridp
  type(msg_three_realdp)::dummy_three_realdp

  associate(r=>s%r,g=>s%g,m=>s%m)

  if(p%static)then
     ! We still need to set the particle levels correctly, even if they are not moved
     ! Loop over particles
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        ! Update level
        p%levelp(ipart)=ilevel
     end do
     return
  end if

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
        wll(idim)=(2D0                        -abs(x(idim)-xll))**3/6D0 ! weight
        wl (idim)=(4D0-6D0*(x(idim)-xl)**2+3d0*abs(x(idim)-xl )**3)/6D0
        wr (idim)=(4D0-6D0*(x(idim)-xr)**2+3d0*abs(x(idim)-xr )**3)/6D0
        wrr(idim)=(2D0                        -abs(x(idim)-xrr))**3/6D0
     end do

     ! Periodic boundary conditions
     do idim=1,ndim
        if(cll(idim)<0)cll(idim)=cll(idim)+m%ckey_max(ilevel+1)
        if(cl (idim)<0)cl (idim)=cl (idim)+m%ckey_max(ilevel+1)
        if(cr (idim)>=m%ckey_max(ilevel+1))cr (idim)=cr (idim)-m%ckey_max(ilevel+1)
        if(crr(idim)>=m%ckey_max(ilevel+1))crr(idim)=crr(idim)-m%ckey_max(ilevel+1)
     enddo

     ! Compute cells Cartesian key
     ckey = pcs_index(cll,cl,cr,crr)

     ! Compute cloud volumes
     vol = pcs_weight(wll,wl,wr,wrr)
     

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
! Tracer-only cache pack/unpack (hydro + grav) to keep default paths intact
! In general, different particle types may need different cache packs/unpacks
!#########################################################################
!#########################################################################
subroutine pack_fetch_kick_trac(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim, ndim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_nvar_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_nvar_realdp)::msg

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        msg%realdp_hydro(ind,ivar)=grid%uold(ind,ivar)
     end do
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_kick_trac
!#########################################################################
!#########################################################################
subroutine unpack_fetch_kick_trac(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_nvar_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  type(msg_nvar_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)
  
#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        grid%uold(ind,ivar)=msg%realdp_hydro(ind,ivar)
     end do
  end do
#endif

end subroutine unpack_fetch_kick_trac
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cic_kick_drift_trac(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, twotondim
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
  real(kind=8),dimension(1:ndim)::x,x_mid,dr,dl,dr2,dl2
  integer,dimension(1:ndim)::ir,il
  integer,dimension(1:ndim)::ir2,il2
  real(kind=8),dimension(1:twotondim)::vol,vol2
  integer,dimension(1:ndim,1:twotondim)::ckey,ckey2
  integer,dimension(1:twotondim)::icell,icell2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,ind,idim
  real(kind=8)::dx_loc,vol_loc
  real(kind=8),dimension(1:ndim)::ff,v_pred
  logical::ok_level
  type(nbor),dimension(1:twotondim)::gridp
  type(msg_three_realdp)::dummy_three_realdp
  type(msg_nvar_realdp)::dummy_nvar_realdp
  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  dx_loc=r%boxlen/2**ilevel
  vol_loc=dx_loc**ndim
  if (p%type/=TRAC_TYPE) return
  ! Tracer hydro cache (uold only)
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_nvar_realdp)/32,&
                     pack=pack_fetch_kick_trac,unpack=unpack_fetch_kick_trac)
  do ipart=p%headp(ilevel),p%tailp(ilevel)
     ! Position in cell units at current level
     do idim=1,ndim
        x(idim)=p%xp(ipart,idim)/dx_loc
     end do
     do idim=1,ndim
        x(idim)=x(idim)-dble(m%ckey_max(ilevel+1))*floor(x(idim)/dble(m%ckey_max(ilevel+1)))
     end do

     ! Gather velocity v = mom/rho at x^n using CIC
     do idim=1,ndim
        dr(idim)=x(idim)+0.5D0
        ir(idim)=int(dr(idim))
        dr(idim)=dr(idim)-ir(idim)
        dl(idim)=1.0D0-dr(idim)
        il(idim)=ir(idim)-1
     end do
     do idim=1,ndim
        if(il(idim)<0)il(idim)=m%ckey_max(ilevel+1)-1
        if(ir(idim)==m%ckey_max(ilevel+1))ir(idim)=0
     enddo
     ckey = cic_index(il,ir)
     ok_level=.true.
     hash_nbor(0)=ilevel+1
     icell=0
     do ind=1,twotondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp(ind)%p,icell(ind),flush_cache=.false.,fetch_cache=.true.,lock=.true.)
        if(.not.associated(gridp(ind)%p)) ok_level=.false.
     end do
     do ind=1,twotondim
        call unlock_cache(s,gridp(ind)%p)
     end do
     if(.not.ok_level)cycle
     vol = cic_weight(dl,dr)
     ff(1:ndim)=0.0
     do ind=1,twotondim
        ff(1:ndim)=ff(1:ndim)+gridp(ind)%p%uold(icell(ind),2:ndim+1)/max(gridp(ind)%p%uold(icell(ind),1), r%smallr)*vol(ind)
     end do

     if(action_part==action_kick_only)then
        ! RK2 step 1 (early call): stash v^n at x^n, no move
        p%vp(ipart,1:ndim)=ff(1:ndim)
        p%levelp(ipart)=ilevel

     else if(action_part==action_kick_drift)then
        ! RK2 step 2 (late call): predict midpoint and correct
        ! Use stored v^n if available; fallback to current ff at step 0
        if (g%nstep>0) then
           v_pred(1:ndim)=p%vp(ipart,1:ndim)
        else
           v_pred(1:ndim)=ff(1:ndim)
        endif
        ! Predict x_mid in cell units
        do idim=1,ndim
           x_mid(idim)=x(idim)+0.5d0*g%dtnew(ilevel)*v_pred(idim)/dx_loc
        end do
        do idim=1,ndim
           x_mid(idim)=x_mid(idim)-dble(m%ckey_max(ilevel+1))*floor(x_mid(idim)/dble(m%ckey_max(ilevel+1)))
        end do

        ! Gather v^{n+1} at x_mid using CIC
        do idim=1,ndim
           dr2(idim)=x_mid(idim)+0.5D0
           ir2(idim)=int(dr2(idim))
           dr2(idim)=dr2(idim)-ir2(idim)
           dl2(idim)=1.0D0-dr2(idim)
           il2(idim)=ir2(idim)-1
        end do
        do idim=1,ndim
           if(il2(idim)<0)il2(idim)=m%ckey_max(ilevel+1)-1
           if(ir2(idim)==m%ckey_max(ilevel+1))ir2(idim)=0
        enddo
        ckey2 = cic_index(il2,ir2)
        ok_level=.true.
        hash_nbor(0)=ilevel+1
        icell2=0
        do ind=1,twotondim
           hash_nbor(1:ndim)=ckey2(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp(ind)%p,icell2(ind),flush_cache=.false.,fetch_cache=.true.,lock=.true.)
           if(.not.associated(gridp(ind)%p)) ok_level=.false.
        end do
        do ind=1,twotondim
           call unlock_cache(s,gridp(ind)%p)
        end do
        if(.not.ok_level)cycle
        vol2 = cic_weight(dl2,dr2)
        ff(1:ndim)=0.0
        do ind=1,twotondim
           ff(1:ndim)=ff(1:ndim)+gridp(ind)%p%uold(icell2(ind),2:ndim+1)/max(gridp(ind)%p%uold(icell2(ind),1), r%smallr)*vol2(ind)
        end do
        ! Set time-centered velocity and drift
        p%vp(ipart,1:ndim)=ff(1:ndim)
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)
     endif
  end do
  call close_cache(s,m%grid_dict)
  if(action_part==action_kick_drift)then
     ! Periodic wrap in physical units
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(p%xp(ipart,idim)<0.0d0)p%xp(ipart,idim)=p%xp(ipart,idim)+r%boxlen
           if(p%xp(ipart,idim)>=r%boxlen)p%xp(ipart,idim)=p%xp(ipart,idim)-r%boxlen
        end do
     end do
  end if
  end associate
end subroutine cic_kick_drift_trac

subroutine tsc_kick_drift_trac(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, threetondim
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
  real(kind=8),dimension(1:ndim)::x,wl,wc,wr
  real(kind=8),dimension(1:ndim)::x_mid,wl2,wc2,wr2
  integer,dimension(1:ndim)::cl,cc,cr
  integer,dimension(1:ndim)::cl2,cc2,cr2
  real(kind=8),dimension(1:threetondim)::vol,vol2
  integer,dimension(1:ndim,1:threetondim)::ckey,ckey2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,icell,icell2,ind,idim
  real(kind=8)::xl,xc,xr
  real(kind=8)::dx_loc
  real(kind=8),dimension(1:ndim)::ff,v_pred
  type(oct),pointer::gridp
  type(msg_three_realdp)::dummy_three_realdp
  type(msg_nvar_realdp)::dummy_nvar_realdp
  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  if (p%type/=TRAC_TYPE) return
  dx_loc=r%boxlen/2**ilevel
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_nvar_realdp)/32,&
                     pack=pack_fetch_kick_trac,unpack=unpack_fetch_kick_trac)
  do ipart=p%headp(ilevel),p%tailp(ilevel)
     do idim=1,ndim
        x(idim)=p%xp(ipart,idim)/dx_loc
     end do
     do idim=1,ndim
        if(x(idim)<0d0)x(idim)=x(idim)+dble(m%ckey_max(ilevel+1))
        if(x(idim)>=dble(m%ckey_max(ilevel+1)))x(idim)=x(idim)-dble(m%ckey_max(ilevel+1))
     end do
     do idim=1,ndim
        cl(idim)=int(x(idim))-1
        cc(idim)=int(x(idim))
        cr(idim)=int(x(idim))+1
        xl=dble(cl(idim))+0.5D0
        xc=dble(cc(idim))+0.5D0
        xr=dble(cr(idim))+0.5D0
        wl(idim)=0.5D0*(1.5D0-abs(x(idim)-xl))**2
        wc(idim)=0.75D0-         (x(idim)-xc) **2
        wr(idim)=0.5D0*(1.5D0-abs(x(idim)-xr))**2
     end do
     do idim=1,ndim
        if(cl(idim)<0)cl(idim)=m%ckey_max(ilevel+1)-1
        if(cr(idim)==m%ckey_max(ilevel+1))cr(idim)=0
     enddo
     ckey = tsc_index(cl,cc,cr)
     vol = tsc_weight(wl,wc,wr)
     hash_nbor(0)=ilevel+1
     ff(1:ndim)=0.0
     do ind=1,threetondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
        if(associated(gridp))then
           ff(1:ndim)=ff(1:ndim)+gridp%uold(icell,2:ndim+1)/max(gridp%uold(icell,1), r%smallr)*vol(ind)
        end if
     end do

     if(action_part==action_kick_only)then
        ! RK2 step 1: stash v^n at x^n
        p%vp(ipart,1:ndim)=ff(1:ndim)
        p%levelp(ipart)=ilevel

     else if(action_part==action_kick_drift)then
        ! RK2 step 2: predict with v^n, sample v^{n+1} at midpoint
        if (g%nstep>0) then
           v_pred(1:ndim)=p%vp(ipart,1:ndim)
        else
           v_pred(1:ndim)=ff(1:ndim)
        endif
        do idim=1,ndim
           x_mid(idim)=x(idim)+0.5d0*g%dtnew(ilevel)*v_pred(idim)/dx_loc
        end do
        do idim=1,ndim
           if(x_mid(idim)<0d0)x_mid(idim)=x_mid(idim)+dble(m%ckey_max(ilevel+1))
           if(x_mid(idim)>=dble(m%ckey_max(ilevel+1)))x_mid(idim)=x_mid(idim)-dble(m%ckey_max(ilevel+1))
        end do
        do idim=1,ndim
           cl2(idim)=int(x_mid(idim))-1
           cc2(idim)=int(x_mid(idim))
           cr2(idim)=int(x_mid(idim))+1
           xl=dble(cl2(idim))+0.5D0
           xc=dble(cc2(idim))+0.5D0
           xr=dble(cr2(idim))+0.5D0
           wl2(idim)=0.5D0*(1.5D0-abs(x_mid(idim)-xl))**2
           wc2(idim)=0.75D0-         (x_mid(idim)-xc) **2
           wr2(idim)=0.5D0*(1.5D0-abs(x_mid(idim)-xr))**2
        end do
        do idim=1,ndim
           if(cl2(idim)<0)cl2(idim)=m%ckey_max(ilevel+1)-1
           if(cr2(idim)==m%ckey_max(ilevel+1))cr2(idim)=0
        enddo
        ckey2 = tsc_index(cl2,cc2,cr2)
        vol2 = tsc_weight(wl2,wc2,wr2)
        hash_nbor(0)=ilevel+1
        ff(1:ndim)=0.0
        do ind=1,threetondim
           hash_nbor(1:ndim)=ckey2(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell2,flush_cache=.false.,fetch_cache=.true.)
           if(associated(gridp))then
              ff(1:ndim)=ff(1:ndim)+gridp%uold(icell2,2:ndim+1)/max(gridp%uold(icell2,1), r%smallr)*vol2(ind)
           end if
        end do
        p%vp(ipart,1:ndim)=ff(1:ndim)
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)
     endif
  end do
  call close_cache(s,m%grid_dict)
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(p%xp(ipart,idim)<0.0d0)p%xp(ipart,idim)=p%xp(ipart,idim)+r%boxlen
           if(p%xp(ipart,idim)>=r%boxlen)p%xp(ipart,idim)=p%xp(ipart,idim)-r%boxlen
        end do
     end do
  end if
  end associate
end subroutine tsc_kick_drift_trac

subroutine pcs_kick_drift_trac(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, fourtondim
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
  real(kind=8),dimension(1:ndim)::x,wll,wl,wr,wrr
  real(kind=8),dimension(1:ndim)::x_mid,wll2,wl2,wr2,wrr2
  integer,dimension(1:ndim)::cll,cl,cr,crr
  integer,dimension(1:ndim)::cll2,cl2,cr2,crr2
  real(kind=8),dimension(1:fourtondim)::vol,vol2
  integer,dimension(1:ndim,1:fourtondim)::ckey,ckey2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,icell,icell2,ind,idim
  real(kind=8)::xll,xl,xr,xrr
  real(kind=8)::dx_loc
  real(kind=8),dimension(1:ndim)::ff,v_pred
  type(oct),pointer::gridp
  type(msg_large_realdp)::dummy_large_realdp
  type(msg_nvar_realdp)::dummy_nvar_realdp
  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  if (p%type/=TRAC_TYPE) return
  dx_loc=r%boxlen/2**ilevel
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
       hilbert=m%domain, pack_size=storage_size(dummy_nvar_realdp)/32,&
       pack=pack_fetch_kick_trac,unpack=unpack_fetch_kick_trac)
  do ipart=p%headp(ilevel),p%tailp(ilevel)
     do idim=1,ndim
        x(idim)=p%xp(ipart,idim)/dx_loc
     end do
     do idim=1,ndim
        if(x(idim)<0d0)x(idim)=x(idim)+dble(m%ckey_max(ilevel+1))
        if(x(idim)>=dble(m%ckey_max(ilevel+1)))x(idim)=x(idim)-dble(m%ckey_max(ilevel+1))
     end do
     do idim=1,ndim
        crr(idim)=int(x(idim)+1.5D0)
        cr (idim)=crr(idim)-1
        cl (idim)=crr(idim)-2
        cll(idim)=crr(idim)-3
        xll=dble(cll(idim))+0.5D0
        xl =dble(cl (idim))+0.5D0
        xr =dble(cr (idim))+0.5D0
        xrr=dble(crr(idim))+0.5D0
        wll(idim)=(2D0-abs(x(idim)-xll))**3/6D0
        wl (idim)=(4D0-6D0*(x(idim)-xl)**2+3d0*abs(x(idim)-xl )**3)/6D0
        wr (idim)=(4D0-6D0*(x(idim)-xr)**2+3d0*abs(x(idim)-xr )**3)/6D0
        wrr(idim)=(2D0-abs(x(idim)-xrr))**3/6D0
     end do
     do idim=1,ndim
        if(cll(idim)<0)cll(idim)=cll(idim)+m%ckey_max(ilevel+1)
        if(cl (idim)<0)cl (idim)=cl (idim)+m%ckey_max(ilevel+1)
        if(cr (idim)>=m%ckey_max(ilevel+1))cr (idim)=cr (idim)-m%ckey_max(ilevel+1)
        if(crr(idim)>=m%ckey_max(ilevel+1))crr(idim)=crr(idim)-m%ckey_max(ilevel+1)
     enddo
     ckey = pcs_index(cll,cl,cr,crr)
     vol = pcs_weight(wll,wl,wr,wrr)
     hash_nbor(0)=ilevel+1
     ff(1:ndim)=0.0
     do ind=1,fourtondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
        if(associated(gridp))then
           ff(1:ndim)=ff(1:ndim)+gridp%uold(icell,2:ndim+1)/max(gridp%uold(icell,1), r%smallr)*vol(ind)
        end if
     end do

     if(action_part==action_kick_only)then
        ! RK2 step 1: stash v^n at x^n
        p%vp(ipart,1:ndim)=ff(1:ndim)
        p%levelp(ipart)=ilevel

     else if(action_part==action_kick_drift)then
        ! RK2 step 2: predict with v^n, sample v^{n+1} at midpoint
        if (g%nstep>0) then
           v_pred(1:ndim)=p%vp(ipart,1:ndim)
        else
           v_pred(1:ndim)=ff(1:ndim)
        endif
        do idim=1,ndim
           x_mid(idim)=x(idim)+0.5d0*g%dtnew(ilevel)*v_pred(idim)/dx_loc
        end do
        do idim=1,ndim
           if(x_mid(idim)<0d0)x_mid(idim)=x_mid(idim)+dble(m%ckey_max(ilevel+1))
           if(x_mid(idim)>=dble(m%ckey_max(ilevel+1)))x_mid(idim)=x_mid(idim)-dble(m%ckey_max(ilevel+1))
        end do
        do idim=1,ndim
           crr2(idim)=int(x_mid(idim)+1.5D0)
           cr2 (idim)=crr2(idim)-1
           cl2 (idim)=crr2(idim)-2
           cll2(idim)=crr2(idim)-3
           xll=dble(cll2(idim))+0.5D0
           xl =dble(cl2 (idim))+0.5D0
           xr =dble(cr2 (idim))+0.5D0
           xrr=dble(crr2(idim))+0.5D0
           wll2(idim)=(2D0-abs(x_mid(idim)-xll))**3/6D0
           wl2 (idim)=(4D0-6D0*(x_mid(idim)-xl)**2+3d0*abs(x_mid(idim)-xl )**3)/6D0
           wr2 (idim)=(4D0-6D0*(x_mid(idim)-xr)**2+3d0*abs(x_mid(idim)-xr )**3)/6D0
           wrr2(idim)=(2D0-abs(x_mid(idim)-xrr))**3/6D0
        end do
        do idim=1,ndim
           if(cll2(idim)<0)cll2(idim)=cll2(idim)+m%ckey_max(ilevel+1)
           if(cl2 (idim)<0)cl2 (idim)=cl2 (idim)+m%ckey_max(ilevel+1)
           if(cr2 (idim)>=m%ckey_max(ilevel+1))cr2 (idim)=cr2 (idim)-m%ckey_max(ilevel+1)
           if(crr2(idim)>=m%ckey_max(ilevel+1))crr2(idim)=crr2(idim)-m%ckey_max(ilevel+1)
        enddo
        ckey2 = pcs_index(cll2,cl2,cr2,crr2)
        vol2 = pcs_weight(wll2,wl2,wr2,wrr2)
        hash_nbor(0)=ilevel+1
        ff(1:ndim)=0.0
        do ind=1,fourtondim
           hash_nbor(1:ndim)=ckey2(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell2,flush_cache=.false.,fetch_cache=.true.)
           if(associated(gridp))then
              ff(1:ndim)=ff(1:ndim)+gridp%uold(icell2,2:ndim+1)/max(gridp%uold(icell2,1), r%smallr)*vol2(ind)
           end if
        end do
        p%vp(ipart,1:ndim)=ff(1:ndim)
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)
     endif
  end do
  call close_cache(s,m%grid_dict)
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(p%xp(ipart,idim)<0.0d0)p%xp(ipart,idim)=p%xp(ipart,idim)+r%boxlen
           if(p%xp(ipart,idim)>=r%boxlen)p%xp(ipart,idim)=p%xp(ipart,idim)-r%boxlen
        end do
     end do
  end if
  end associate
end subroutine pcs_kick_drift_trac

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module move_fine_module
