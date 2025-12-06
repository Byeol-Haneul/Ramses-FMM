module move_fine_module
  use rho_fine_module, only: cic_weight, cic_index, tsc_weight, tsc_index, pcs_weight, pcs_index
  use rng
  type(RngStream), save :: tracer_rng
  logical, save :: tracer_rng_ready = .false.
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
        if(pst%s%r%trac_interpolation_scheme==0)then
           call cic_trace_gas_part_ito(pst%s,pst%s%trac,ilevel,action_part)
        elseif(pst%s%r%trac_interpolation_scheme==1)then
           call cic_trace_gas_part(pst%s,pst%s%trac,ilevel,action_part)
        elseif(pst%s%r%trac_interpolation_scheme==2)then
           call tsc_trace_gas_part(pst%s,pst%s%trac,ilevel,action_part)
        elseif(pst%s%r%trac_interpolation_scheme==3)then
           call pcs_trace_gas_part(pst%s,pst%s%trac,ilevel,action_part)
        elseif(pst%s%r%trac_interpolation_scheme==4)then
           call cic_trace_gas_part_num(pst%s,pst%s%trac,ilevel,action_part)
        endif
     endif
     if(pst%s%r%dust)then
        if(pst%s%r%dust_force_interpolation_scheme==1)then
           call cic_kick_drift_dust(pst%s,pst%s%dust,ilevel,action_part)
        elseif(pst%s%r%dust_force_interpolation_scheme==2)then
           call tsc_kick_drift_dust(pst%s,pst%s%dust,ilevel,action_part)
        elseif(pst%s%r%dust_force_interpolation_scheme==3)then
           call pcs_kick_drift_dust(pst%s,pst%s%dust,ilevel,action_part)
        endif
     endif
  endif

end subroutine r_kick_drift_part
!################################################################
!################################################################
!################################################################
!################################################################
subroutine cic_kick_drift_part(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, twotondim, nvector
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
  real(kind=8),dimension(1:nvector,1:ndim)::xana
  real(kind=8),dimension(1:nvector,1:ndim)::fana
  type(nbor),dimension(1:twotondim)::gridp
  type(msg_three_realdp)::dummy_three_realdp

  associate(r=>s%r,g=>s%g,m=>s%m)

  if(p%static)then
     ! We still need to set the particle levels correctly, even if they are not moved
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        p%levelp(ipart)=ilevel
     end do
     return
  end if

  ! Deal with particles that left the computational domain
  if(ilevel==r%levelmin.and.ANY(.not.r%periodic(1:ndim)))then

     ! Loop over particles outside the box
     do ipart=p%headp(ilevel-1),p%tailp(ilevel-1)

        ! Get particle position
        do idim=1,ndim
           x(idim)=p%xp(ipart,idim)
        end do

        ! Call analytical acceleration routine
        xana(1,1:ndim)=x(1:ndim)
        call gravana(r,g,xana,fana,1)
        ff(1:ndim)=fana(1,1:ndim)

        ! Perform kick, or drift, or both
        if(action_part==action_kick_drift)then

           ! Update velocity (use levelmin time step)
           p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)

           ! Update position (use levelmin time step)
           p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)

        else if(action_part.EQ.action_kick_only)then

           ! Compute proper time step for second kick
           dteff=g%dtnew(p%levelp(ipart))

           ! Update level to levelmin
           p%levelp(ipart)=ilevel

           ! Update velocity
           p%vp(ipart,1:ndim)=p%vp(ipart,1:ndim)+ff(1:ndim)*0.5d0*dteff

        endif

     end do
     ! End loop over particles

  endif

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Open read-only cache
  call open_cache(s,table=m%grid_dict, data_size=storage_size(m%grid(1))/32,&
       hilbert=m%domain, pack_size=storage_size(dummy_three_realdp)/32,&
       pack=pack_fetch_kick, unpack=unpack_fetch_kick)

  ! Loop over particles
  do ipart=p%headp(ilevel),p%tailp(ilevel)

     ! Rescale particle position at level ilevel
     do idim=1,ndim
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
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
        if(r%periodic(idim))then
           if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=m%box_ckey_max(idim,ilevel+1)-1
           if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=m%box_ckey_min(idim,ilevel+1)
        endif
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
           if(r%periodic(idim))then
              if(il(idim)< m%box_ckey_min(idim,ilevel))il(idim)=m%box_ckey_max(idim,ilevel)-1
              if(ir(idim)>=m%box_ckey_max(idim,ilevel))ir(idim)=m%box_ckey_min(idim,ilevel)
           endif
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
#ifdef GRAV
        do ind=1,twotondim
           ff(1:ndim)=ff(1:ndim)+gridp(ind)%p%f(icell(ind),1:ndim)*vol(ind)
        end do
#ifdef OUTPUT_PARTICLE_POTENTIAL
        p%phip(ipart)=0.0
        do ind=1,twotondim
           p%phip(ipart)=p%phip(ipart)+gridp(ind)%p%phi(icell(ind))*vol(ind)
        end do
#endif
#endif
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
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           end if
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
     do ipart=p%headp(ilevel),p%tailp(ilevel)
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
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
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
        if(r%periodic(idim))then
           if(cl(idim)< m%box_ckey_min(idim,ilevel+1))cl(idim)=m%box_ckey_max(idim,ilevel+1)-1
           if(cr(idim)>=m%box_ckey_max(idim,ilevel+1))cr(idim)=m%box_ckey_min(idim,ilevel+1)
        endif
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
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           end if
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
     do ipart=p%headp(ilevel),p%tailp(ilevel)
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
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
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
        if(r%periodic(idim))then
           if(cll(idim)< m%box_ckey_min(idim,ilevel+1))cll(idim)=cll(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
           if(cl (idim)< m%box_ckey_min(idim,ilevel+1))cl (idim)=cl (idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
           if(cr (idim)>=m%box_ckey_max(idim,ilevel+1))cr (idim)=cr (idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
           if(crr(idim)>=m%box_ckey_max(idim,ilevel+1))crr(idim)=crr(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
        endif
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
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           end if
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
  use oct_commons, only: oct
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
#ifdef HYDRO
  msg%realdp_mflux=grid%mflux
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_kick
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine unpack_fetch_kick(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use oct_commons, only: oct
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
#ifdef HYDRO
  grid%mflux=msg%realdp_mflux
#endif

end subroutine unpack_fetch_kick
!#########################################################################
!#########################################################################
! Tracer-only cache pack/unpack (hydro only) to keep default paths intact
! In general, different particle types may need different cache packs/unpacks
!#########################################################################
!#########################################################################
subroutine pack_fetch_kick_trac(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim, ndim
  use hydro_parameters, only: nvar
  use oct_commons, only: oct
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
  do ind=1,twotondim
     msg%realdp_mflux(ind,1:6)=grid%mflux(ind,1:6)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_kick_trac
!#########################################################################
!#########################################################################
subroutine unpack_fetch_kick_trac(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use oct_commons, only: oct
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
  do ind=1,twotondim
     grid%mflux(ind,1:6)=msg%realdp_mflux(ind,1:6)
  end do
#endif

end subroutine unpack_fetch_kick_trac
!#########################################################################
!#########################################################################
subroutine pack_fetch_kick_dust(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use hydro_parameters, only: nvar
  use oct_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_large_realdp)::msg

#ifdef HYDRO
  do ind=1,twotondim
     do ivar=1,nvar
        msg%realdp_hydro(ind,ivar)=grid%uold(ind,ivar)
     end do
  end do
#endif
#ifdef GRAV
  do ind=1,twotondim
     msg%realdp_phi(ind)=grid%f(ind,1)
     msg%realdp_phi_old(ind)=grid%f(ind,2)
     msg%realdp_dis(ind)=grid%f(ind,3)
  end do
#endif
#ifdef RT
  do ind=1,twotondim
     do ivar=1,nrtvar
        msg%realdp_rt(ind,ivar)=grid%rtuold(ind,ivar)
     end do
  end do
#endif
#ifdef MHD
  do ivar=1,6
     do ind=1,twotondim
        msg%realdp_mhd(ind,ivar)=grid%bold(ind,ivar)
     end do
  end do
#endif
  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_kick_dust
!#########################################################################
!#########################################################################
subroutine unpack_fetch_kick_dust(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use oct_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  type(msg_large_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef HYDRO
  do ind=1,twotondim
     do ivar=1,nvar
        grid%uold(ind,ivar)=msg%realdp_hydro(ind,ivar)
     end do
  end do
#endif
#ifdef GRAV
  do ind=1,twotondim
     grid%f(ind,1)=msg%realdp_phi(ind)
     grid%f(ind,2)=msg%realdp_phi_old(ind)
     grid%f(ind,3)=msg%realdp_dis(ind)
  end do
#endif
#ifdef RT
  do ind=1,twotondim
     do ivar=1,nrtvar
        grid%rtuold(ind,ivar)=msg%realdp_rt(ind,ivar)
     end do
  end do
#endif
#ifdef MHD
  do ivar=1,6
     do ind=1,twotondim
        grid%bold(ind,ivar)=msg%realdp_mhd(ind,ivar)
     end do
  end do
#endif

end subroutine unpack_fetch_kick_dust
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cic_trace_gas_part(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, twotondim
  use pm_parameters
  use pm_commons, only: part_t
  use oct_commons, only: oct
  use ramses_commons, only: ramses_t
  use rng
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  real(kind=8),dimension(1:ndim)::x,x_mid,v_pred,vel,vel_mid,disp,xi
  real(kind=8)::dx_loc,dt_level,kappa_mid,noise_amp
  logical :: use_sgs
  type(msg_nvar_realdp)::dummy_nvar_realdp
  type(RngStream)::RngStream_CreateStream
  real(kind=8)::RngStream_RandUni
  integer(kind=8)::stream_skip
  external :: RngStream_SetPackageSeed, RngStream_AdvanceState, gaussdev
  integer :: ipart,idim

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  if (p%type/=TRAC_TYPE) return

  dx_loc=r%boxlen/2**ilevel
  dt_level=g%dtnew(ilevel)
  use_sgs = r%sgs_turb .and. (r%iturb>0)
  if(use_sgs .and. .not. tracer_rng_ready)then
     call RngStream_SetPackageSeed(r%seed)
     tracer_rng = RngStream_CreateStream('tracer_sgs')
     stream_skip = int(2*g%myid,kind=8)
     call RngStream_AdvanceState(tracer_rng,0_8,stream_skip)
     tracer_rng_ready = .true.
  end if

  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_nvar_realdp)/32,&
                     pack=pack_fetch_kick_trac,unpack=unpack_fetch_kick_trac)

  do ipart=p%headp(ilevel),p%tailp(ilevel)

     do idim=1,ndim
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
     end do
     call wrap_cell_coords(s,x,ilevel+1)

     call gather_cic_state(s,x,ilevel,dx_loc,use_sgs,vel,kappa_mid)

     if(action_part==action_kick_only)then
        p%vp(ipart,1:ndim)=vel(1:ndim)
        p%levelp(ipart)=ilevel
        cycle
     endif

     if (g%nstep>0) then
        v_pred(1:ndim)=p%vp(ipart,1:ndim)
     else
        v_pred(1:ndim)=vel(1:ndim)
     endif

     do idim=1,ndim
        x_mid(idim)=x(idim)+0.5d0*dt_level*v_pred(idim)/dx_loc
     end do
     call wrap_cell_coords(s,x_mid,ilevel+1)

     call gather_cic_state(s,x_mid,ilevel,dx_loc,use_sgs,vel_mid,kappa_mid)

     p%vp(ipart,1:ndim)=vel_mid(1:ndim)
     disp(1:ndim)=p%vp(ipart,1:ndim)*dt_level
     if(use_sgs .and. kappa_mid>0.0d0)then
        call sample_tracer_gaussian(xi)
        noise_amp = sqrt(2.0d0*kappa_mid*dt_level)
        disp(1:ndim)=disp(1:ndim)+noise_amp*xi(1:ndim)
     end if
     p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+disp(1:ndim)
  end do

  call close_cache(s,m%grid_dict)

  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0)p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           endif
        end do
     end do
  end if

  end associate

end subroutine cic_trace_gas_part

subroutine cic_trace_gas_part_ito(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, twotondim
  use pm_parameters
  use pm_commons, only: part_t
  use oct_commons, only: oct
  use ramses_commons, only: ramses_t
  use rng
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  real(kind=8),dimension(1:ndim)::x,x_mid,v_pred,vel,vel_mid,disp,xi
  real(kind=8),dimension(1:twotondim)::kappa_cells
  real(kind=8),dimension(1:ndim,1:twotondim)::grad_kappa_cells
  real(kind=8)::dx_loc,dt_level,kappa_mid,noise_amp
  logical :: use_sgs
  type(msg_nvar_realdp)::dummy_nvar_realdp
  type(RngStream)::RngStream_CreateStream
  real(kind=8)::RngStream_RandUni
  integer(kind=8)::stream_skip
  external :: RngStream_SetPackageSeed, RngStream_AdvanceState, gaussdev
  integer :: ipart,idim,ic

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  if (p%type/=TRAC_TYPE) return

  dx_loc=r%boxlen/2**ilevel
  dt_level=g%dtnew(ilevel)
  use_sgs = r%sgs_turb .and. (r%iturb>0)
  if(use_sgs .and. .not. tracer_rng_ready)then
     call RngStream_SetPackageSeed(r%seed)
     tracer_rng = RngStream_CreateStream('tracer_sgs')
     stream_skip = int(2*g%myid,kind=8)
     call RngStream_AdvanceState(tracer_rng,0_8,stream_skip)
     tracer_rng_ready = .true.
  end if

  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_nvar_realdp)/32,&
                     pack=pack_fetch_kick_trac,unpack=unpack_fetch_kick_trac)

  do ipart=p%headp(ilevel),p%tailp(ilevel)

     do idim=1,ndim
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
     end do
     call wrap_cell_coords(s,x,ilevel+1)

     call gather_cic_state(s,x,ilevel,dx_loc,use_sgs,vel,kappa_mid)

     ! also gather cell-centered kappa for gradient estimates
     do ic=1,twotondim
        kappa_cells(ic)=0.d0
     end do
     call gather_cic_scalar(s,x,ilevel,dx_loc,use_sgs,kappa_cells)
     call compute_cell_gradients(kappa_cells,dx_loc,grad_kappa_cells)

     if(action_part==action_kick_only)then
        p%vp(ipart,1:ndim)=vel(1:ndim)
        p%levelp(ipart)=ilevel
        cycle
     endif

     if (g%nstep>0) then
        v_pred(1:ndim)=p%vp(ipart,1:ndim)
     else
        v_pred(1:ndim)=vel(1:ndim)
     endif

     do idim=1,ndim
        x_mid(idim)=x(idim)+0.5d0*dt_level*v_pred(idim)/dx_loc
     end do
     call wrap_cell_coords(s,x_mid,ilevel+1)

     call gather_cic_state(s,x_mid,ilevel,dx_loc,use_sgs,vel_mid,kappa_mid)

     ! interpolate grad(kappa) at midpoint
     disp(1:ndim)=0.d0
     call interp_grad_at_pos(s,x_mid,ilevel,grad_kappa_cells,disp)

     ! Ito drift: u + grad(kappa)
     do idim=1,ndim
        disp(idim) = (vel_mid(idim) + disp(idim))*dt_level
     end do

     if(use_sgs .and. kappa_mid>0.0d0)then
        call sample_tracer_gaussian(xi)
        noise_amp = sqrt(2.0d0*kappa_mid*dt_level)
        disp(1:ndim)=disp(1:ndim)+noise_amp*xi(1:ndim)
     end if

     p%vp(ipart,1:ndim)=vel_mid(1:ndim)
     p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+disp(1:ndim)
  end do

  call close_cache(s,m%grid_dict)

  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0)p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           endif
        end do
     end do
  end if

  end associate
end subroutine cic_trace_gas_part_ito

subroutine cic_trace_gas_part_num(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, twotondim
  use pm_parameters
  use pm_commons, only: part_t
  use oct_commons, only: oct
  use ramses_commons, only: ramses_t
  use rng
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  real(kind=8),dimension(1:ndim)::x,disp,xi,u_eff,d_eff
  real(kind=8),dimension(1:ndim)::dl,dr
  real(kind=8),dimension(1:ndim)::grad_at_p
  integer,dimension(1:ndim)::il,ir
  real(kind=8),dimension(1:twotondim)::vol,phi_slice
  integer,dimension(1:ndim,1:twotondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  real(kind=8),dimension(1:ndim,1:twotondim)::u_cells,d_cells,grad_phi_cells
  type(oct),pointer::gridp
  integer :: ipart,ind,idim,icell,k
  real(kind=8)::dx_loc,dt_level,rho,denom,fluxL,fluxR,jr,jl,noise_amp
  type(msg_nvar_realdp)::dummy_nvar_realdp
  type(RngStream)::RngStream_CreateStream
  real(kind=8)::RngStream_RandUni
  integer(kind=8)::stream_skip
  external :: RngStream_SetPackageSeed, RngStream_AdvanceState, gaussdev

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  if (p%type/=TRAC_TYPE) return

  dx_loc=r%boxlen/2**ilevel
  dt_level=g%dtnew(ilevel)

  if(.not.tracer_rng_ready)then
     call RngStream_SetPackageSeed(r%seed)
     tracer_rng = RngStream_CreateStream('tracer_sgs')
     stream_skip = int(2*g%myid,kind=8)
     call RngStream_AdvanceState(tracer_rng,0_8,stream_skip)
     tracer_rng_ready = .true.
  end if

  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_nvar_realdp)/32,&
                     pack=pack_fetch_kick_trac,unpack=unpack_fetch_kick_trac)

  do ipart=p%headp(ilevel),p%tailp(ilevel)

     do idim=1,ndim
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
     end do
     call wrap_cell_coords(s,x,ilevel+1)

     do idim=1,ndim
        dr(idim)=x(idim)+0.5d0
        ir(idim)=int(dr(idim))
        dr(idim)=dr(idim)-ir(idim)
        dl(idim)=1.0d0-dr(idim)
        il(idim)=ir(idim)-1
     end do
     do idim=1,ndim
        if(r%periodic(idim))then
           if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=m%box_ckey_max(idim,ilevel+1)-1
           if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=m%box_ckey_min(idim,ilevel+1)
        endif
     enddo

     ckey = cic_index(il,ir)
     vol = cic_weight(dl,dr)

     u_cells=0.d0
     d_cells=0.d0

     hash_nbor(0)=ilevel+1
     do ind=1,twotondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
        if(associated(gridp))then
           rho=gridp%uold(icell,1)
           denom=max(rho,r%smallr)
           do idim=1,ndim
              ! mflux stores time-integrated flux ~ (dt/dx)*F; recover physical flux F with factor dx_loc/dt_level
              fluxL=gridp%mflux(icell,idim     )*dx_loc/dt_level
              fluxR=gridp%mflux(icell,idim+ndim)*dx_loc/dt_level
              jr=max(fluxR,0.d0)
              jl=max(-fluxL,0.d0)
              u_cells(idim,ind)=(jr-jl)/denom
              d_cells(idim,ind)=0.5d0*(jr+jl)/denom*dx_loc
           end do
        end if
#endif
     end do

     u_eff=0.d0
     d_eff=0.d0
     do ind=1,twotondim
        do idim=1,ndim
           u_eff(idim)=u_eff(idim)+u_cells(idim,ind)*vol(ind)
           d_eff(idim)=d_eff(idim)+d_cells(idim,ind)*vol(ind)
        end do
     end do

     do k=1,ndim
        do ind=1,twotondim
           phi_slice(ind)=d_cells(k,ind)
        end do
        call compute_cell_gradients(phi_slice,dx_loc,grad_phi_cells)
        call interp_grad_at_pos(s,x,ilevel,grad_phi_cells,grad_at_p)
        u_eff(k) = u_eff(k) + (r%tracer_schmidt_number - 1.0d0) * grad_at_p(k)
        d_eff(k) = d_eff(k) * r%tracer_schmidt_number
     end do

     if(action_part==action_kick_only)then
        p%vp(ipart,1:ndim)=u_eff(1:ndim)
        p%levelp(ipart)=ilevel
        cycle
     endif

     !call sample_tracer_gaussian(xi)
     call sample_tracer_uniform(xi)
     do idim=1,ndim
        disp(idim)=u_eff(idim)*dt_level
        noise_amp = sqrt(max(0.d0,2.d0*d_eff(idim)*dt_level))
        disp(idim)=disp(idim)+noise_amp*xi(idim)
     end do

     p%vp(ipart,1:ndim)=u_eff(1:ndim)
     p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+disp(1:ndim)
  end do

  call close_cache(s,m%grid_dict)

  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           endif
        end do
     end do
  end if

  end associate
end subroutine cic_trace_gas_part_num

subroutine wrap_cell_coords(st,x_cell,levelp1)
  use amr_parameters, only: ndim
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t),intent(in)::st
  real(kind=8),intent(inout)::x_cell(1:ndim)
  integer,intent(in)::levelp1
  integer::jd
  real(kind=8)::range

  do jd=1,ndim
     if(st%r%periodic(jd))then
        range=dble(st%m%box_ckey_max(jd,levelp1)-st%m%box_ckey_min(jd,levelp1))
        if(range<=0.d0)cycle
        if(x_cell(jd)< dble(st%m%box_ckey_min(jd,levelp1)))x_cell(jd)=x_cell(jd)+range
        if(x_cell(jd)>=dble(st%m%box_ckey_max(jd,levelp1)))x_cell(jd)=x_cell(jd)-range
     endif
  end do
end subroutine wrap_cell_coords

subroutine gather_cic_state(st,x_cell,level_in,dx_cell,use_sgs_in,vel_out,kappa_out)
  use amr_parameters, only: ndim, twotondim
  use oct_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache
  implicit none
  type(ramses_t),intent(in)::st
  real(kind=8),intent(in)::x_cell(1:ndim)
  integer,intent(in)::level_in
  real(kind=8),intent(in)::dx_cell
  logical,intent(in)::use_sgs_in
  real(kind=8),intent(out)::vel_out(1:ndim)
  real(kind=8),intent(out)::kappa_out
  real(kind=8),dimension(1:ndim)::dl,dr
  integer,dimension(1:ndim)::il,ir
  real(kind=8),dimension(1:twotondim)::vol
  integer,dimension(1:ndim,1:twotondim)::ckey
  real(kind=8),dimension(1:ndim)::momentum
  real(kind=8)::rho,kappa_sum
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ind,icell,jd
  type(oct),pointer::gridp

  vel_out=0.d0
  momentum=0.d0
  rho=0.d0
  kappa_sum=0.d0

  do jd=1,ndim
     dr(jd)=x_cell(jd)+0.5d0
     ir(jd)=int(dr(jd))
     dr(jd)=dr(jd)-ir(jd)
     dl(jd)=1.0d0-dr(jd)
     il(jd)=ir(jd)-1
  end do
  do jd=1,ndim
     if(st%r%periodic(jd))then
        if(il(jd)< st%m%box_ckey_min(jd,level_in+1))il(jd)=st%m%box_ckey_max(jd,level_in+1)-1
        if(ir(jd)>=st%m%box_ckey_max(jd,level_in+1))ir(jd)=st%m%box_ckey_min(jd,level_in+1)
     endif
  end do
  ckey = cic_index(il,ir)
  vol = cic_weight(dl,dr)

  hash_nbor(0)=level_in+1
  do ind=1,twotondim
     hash_nbor(1:ndim)=ckey(1:ndim,ind)
     call get_parent_cell(st,hash_nbor,st%m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
     if(associated(gridp))then
        momentum(1:ndim)=momentum(1:ndim)+gridp%uold(icell,2:ndim+1)*vol(ind)
        rho=rho+gridp%uold(icell,1)*vol(ind)
        if(use_sgs_in)then
           kappa_sum=kappa_sum+tracer_cell_kappa(gridp%uold(icell,1),gridp%uold(icell,st%r%iturb),dx_cell,st%r%smallr,st%r%tracer_schmidt_number)*vol(ind)
        end if
     end if
#endif
  end do

  if(rho>0.d0)then
     vel_out(1:ndim)=momentum(1:ndim)/max(rho,st%r%smallr)
  else
     vel_out(1:ndim)=0.d0
  end if

  if(use_sgs_in)then
     kappa_out=max(kappa_sum,0.d0)
  else
     kappa_out=0.d0
  end if
end subroutine gather_cic_state

subroutine gather_cic_scalar(st,x_cell,level_in,dx_cell,use_sgs_in,phi_cells)
  use amr_parameters, only: ndim, twotondim
  use oct_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache
  implicit none
  type(ramses_t),intent(in)::st
  real(kind=8),intent(in)::x_cell(1:ndim)
  integer,intent(in)::level_in
  real(kind=8),intent(in)::dx_cell
  logical,intent(in)::use_sgs_in
  real(kind=8),intent(out)::phi_cells(1:twotondim)
  real(kind=8),dimension(1:ndim)::dl,dr
  integer,dimension(1:ndim)::il,ir
  integer,dimension(1:ndim,1:twotondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ind,jd,icell
  type(oct),pointer::gridp

  phi_cells=0.d0
  if(.not.use_sgs_in)return

  do jd=1,ndim
     dr(jd)=x_cell(jd)+0.5d0
     ir(jd)=int(dr(jd))
     dr(jd)=dr(jd)-ir(jd)
     dl(jd)=1.0d0-dr(jd)
     il(jd)=ir(jd)-1
  end do
  do jd=1,ndim
     if(st%r%periodic(jd))then
        if(il(jd)< st%m%box_ckey_min(jd,level_in+1))il(jd)=st%m%box_ckey_max(jd,level_in+1)-1
        if(ir(jd)>=st%m%box_ckey_max(jd,level_in+1))ir(jd)=st%m%box_ckey_min(jd,level_in+1)
     endif
  end do
  ckey = cic_index(il,ir)

  hash_nbor(0)=level_in+1
  do ind=1,twotondim
     hash_nbor(1:ndim)=ckey(1:ndim,ind)
     call get_parent_cell(st,hash_nbor,st%m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
     if(associated(gridp))then
        phi_cells(ind)=tracer_cell_kappa(gridp%uold(icell,1),gridp%uold(icell,st%r%iturb),dx_cell,st%r%smallr,st%r%tracer_schmidt_number)
     end if
#endif
  end do
end subroutine gather_cic_scalar

subroutine compute_cell_gradients(phi_cell,dx_cell,grad_phi)
  use amr_parameters, only: ndim, twotondim
  implicit none
  real(kind=8),intent(in) :: phi_cell(1:twotondim)
  real(kind=8),intent(in) :: dx_cell
  real(kind=8),intent(out):: grad_phi(1:ndim,1:twotondim)
  integer,parameter :: child_coord(3,8)=reshape(&
       (/0,0,0, 1,0,0, 0,1,0, 1,1,0, 0,0,1, 1,0,1, 0,1,1, 1,1,1/),(/3,8/))
  integer :: c,d,plus_idx,minus_idx
  integer :: cx,cy,cz
  real(kind=8) :: phi_plus,phi_minus

  grad_phi=0.d0

  do c=1,twotondim
     cx=child_coord(1,c)
     cy=child_coord(2,c)
     cz=child_coord(3,c)

      do d=1,ndim
        plus_idx=-1; minus_idx=-1
        select case(d)
        case(1)
           if(cx==0) plus_idx = c+1
           if(cx==1) minus_idx = c-1
        case(2)
           if(cy==0) plus_idx = c+2
           if(cy==1) minus_idx = c-2
        case(3)
           if(cz==0) plus_idx = c+4
           if(cz==1) minus_idx = c-4
        end select

        if(plus_idx>0 .and. plus_idx<=twotondim .and. minus_idx>0 .and. minus_idx<=twotondim)then
           phi_plus = phi_cell(plus_idx)
           phi_minus= phi_cell(minus_idx)
           grad_phi(d,c) = (phi_plus-phi_minus)/(2.d0*dx_cell)
        else
           grad_phi(d,c) = 0.d0
        end if
     end do
  end do
end subroutine compute_cell_gradients

subroutine interp_grad_at_pos(st,x_cell,level_in,grad_cells,grad_out)
  use amr_parameters, only: ndim, twotondim
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t),intent(in)::st
  real(kind=8),intent(in)::x_cell(1:ndim)
  integer,intent(in)::level_in
  real(kind=8),intent(in)::grad_cells(1:ndim,1:twotondim)
  real(kind=8),intent(out)::grad_out(1:ndim)
  real(kind=8),dimension(1:ndim)::dl,dr
  integer,dimension(1:ndim)::il,ir
  real(kind=8),dimension(1:twotondim)::vol
  integer,dimension(1:ndim,1:twotondim)::ckey
  integer :: jd,ind

  grad_out=0.d0

  do jd=1,ndim
     dr(jd)=x_cell(jd)+0.5d0
     ir(jd)=int(dr(jd))
     dr(jd)=dr(jd)-ir(jd)
     dl(jd)=1.0d0-dr(jd)
     il(jd)=ir(jd)-1
  end do
  do jd=1,ndim
     if(st%r%periodic(jd))then
        if(il(jd)< st%m%box_ckey_min(jd,level_in+1))il(jd)=st%m%box_ckey_max(jd,level_in+1)-1
        if(ir(jd)>=st%m%box_ckey_max(jd,level_in+1))ir(jd)=st%m%box_ckey_min(jd,level_in+1)
     endif
  end do
  ckey = cic_index(il,ir)
  vol = cic_weight(dl,dr)

  do ind=1,twotondim
     grad_out(1:ndim)=grad_out(1:ndim)+grad_cells(1:ndim,ind)*vol(ind)
  end do
end subroutine interp_grad_at_pos

real(kind=8) function tracer_cell_kappa(dens_in,eturb_in,dx_in,smallr_in,schmidt_in) result(kappa_val)
  implicit none
  real(kind=8),intent(in)::dens_in,eturb_in,dx_in,smallr_in,schmidt_in
  real(kind=8)::rho_eff,sigma_sq

  rho_eff = max(dens_in,smallr_in)
  sigma_sq = max(2.0d0*max(eturb_in,0.0d0)/rho_eff,0.0d0)
  if(sigma_sq>0.0d0)then
     kappa_val = schmidt_in*dx_in*sqrt(sigma_sq)
  else
     kappa_val = 0.0d0
  end if
end function tracer_cell_kappa

subroutine sample_tracer_gaussian(vec)
  use amr_parameters, only: ndim
  implicit none
  real(kind=8),intent(out)::vec(1:ndim)
  integer :: jd
  real(kind=8)::u_rand,tmp
  real(kind=8), external :: RngStream_RandUni
  external :: gaussdev

  vec=0.0d0
  do jd=1,ndim
     u_rand = RngStream_RandUni(tracer_rng)
     call gaussdev(u_rand,tmp)
     vec(jd)=tmp
  end do
end subroutine sample_tracer_gaussian

subroutine sample_tracer_uniform(vec)
  use amr_parameters, only: ndim
  implicit none
  real(kind=8),intent(out)::vec(1:ndim)
  integer :: jd
  real(kind=8)::u_rand
  real(kind=8), external :: RngStream_RandUni

  vec=0.0d0
  do jd=1,ndim
     u_rand = RngStream_RandUni(tracer_rng)
     ! Uniform distribution [-sqrt(3), sqrt(3)]
     ! variance = (sqrt(3) - (-sqrt(3)))^2 / 12 = 12/12 = 1.
     ! mean = 0
     vec(jd) = (2.0d0*u_rand - 1.0d0) * sqrt(3.0d0)
  end do
end subroutine sample_tracer_uniform

subroutine tsc_trace_gas_part(s,p,ilevel,action_part)
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
  logical::ok_level
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
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
     end do
     do idim=1,ndim
        if(r%periodic(idim))then
           if(x(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x(idim)=x(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           if(x(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x(idim)=x(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
        endif
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
        if(r%periodic(idim))then
           if(cl(idim)< m%box_ckey_min(idim,ilevel+1))cl(idim)=m%box_ckey_max(idim,ilevel+1)-1
           if(cr(idim)>=m%box_ckey_max(idim,ilevel+1))cr(idim)=m%box_ckey_min(idim,ilevel+1)
        endif
     enddo
     ckey = tsc_index(cl,cc,cr)
     vol = tsc_weight(wl,wc,wr)
     hash_nbor(0)=ilevel+1
     ff(1:ndim)=0.0
     do ind=1,threetondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
        if(associated(gridp))then
           ff(1:ndim)=ff(1:ndim)+gridp%uold(icell,2:ndim+1)/max(gridp%uold(icell,1), r%smallr)*vol(ind)
        end if
#endif
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
           if(r%periodic(idim))then
              if(x_mid(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x_mid(idim)=x_mid(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
              if(x_mid(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x_mid(idim)=x_mid(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           endif
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
           if(r%periodic(idim))then
              if(cl2(idim)< m%box_ckey_min(idim,ilevel+1))cl2(idim)=m%box_ckey_max(idim,ilevel+1)-1
              if(cr2(idim)>=m%box_ckey_max(idim,ilevel+1))cr2(idim)=m%box_ckey_min(idim,ilevel+1)
           endif
        enddo
        ckey2 = tsc_index(cl2,cc2,cr2)
        vol2 = tsc_weight(wl2,wc2,wr2)
        hash_nbor(0)=ilevel+1
        ff(1:ndim)=0.0
        do ind=1,threetondim
           hash_nbor(1:ndim)=ckey2(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell2,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
           if(associated(gridp))then
              ff(1:ndim)=ff(1:ndim)+gridp%uold(icell2,2:ndim+1)/max(gridp%uold(icell2,1), r%smallr)*vol2(ind)
           end if
#endif
        end do
        p%vp(ipart,1:ndim)=ff(1:ndim)
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)
     endif
  end do
  call close_cache(s,m%grid_dict)
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(ipart,idim) <0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           endif
        end do
     end do
  end if
  end associate
end subroutine tsc_trace_gas_part

subroutine pcs_trace_gas_part(s,p,ilevel,action_part)
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
        x(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
     end do
     do idim=1,ndim
        if(r%periodic(idim))then
           if(x(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x(idim)=x(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           if(x(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x(idim)=x(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
        endif
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
        if(r%periodic(idim))then
           if(cll(idim)< m%box_ckey_min(idim,ilevel+1))cll(idim)=cll(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
           if(cl (idim)< m%box_ckey_min(idim,ilevel+1))cl (idim)=cl (idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
           if(cr (idim)>=m%box_ckey_max(idim,ilevel+1))cr (idim)=cr (idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
           if(crr(idim)>=m%box_ckey_max(idim,ilevel+1))crr(idim)=crr(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
        endif
     enddo
     ckey = pcs_index(cll,cl,cr,crr)
     vol = pcs_weight(wll,wl,wr,wrr)
     hash_nbor(0)=ilevel+1
     ff(1:ndim)=0.0
     do ind=1,fourtondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
        if(associated(gridp))then
           ff(1:ndim)=ff(1:ndim)+gridp%uold(icell,2:ndim+1)/max(gridp%uold(icell,1), r%smallr)*vol(ind)
        end if
#endif
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
           if(r%periodic(idim))then
              if(x_mid(idim)< dble(m%box_ckey_min(idim,ilevel+1)))x_mid(idim)=x_mid(idim)+dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
              if(x_mid(idim)>=dble(m%box_ckey_max(idim,ilevel+1)))x_mid(idim)=x_mid(idim)-dble(m%box_ckey_max(idim,ilevel+1)-m%box_ckey_min(idim,ilevel+1))
           endif
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
           if(r%periodic(idim))then
              if(cll2(idim)< m%box_ckey_min(idim,ilevel+1))cll2(idim)=cll2(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
              if(cl2 (idim)< m%box_ckey_min(idim,ilevel+1))cl2 (idim)=cl2 (idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
              if(cr2 (idim)>=m%box_ckey_max(idim,ilevel+1))cr2 (idim)=cr2 (idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
              if(crr2(idim)>=m%box_ckey_max(idim,ilevel+1))crr2(idim)=crr2(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
           endif
        enddo
        ckey2 = pcs_index(cll2,cl2,cr2,crr2)
        vol2 = pcs_weight(wll2,wl2,wr2,wrr2)
        hash_nbor(0)=ilevel+1
        ff(1:ndim)=0.0
        do ind=1,fourtondim
           hash_nbor(1:ndim)=ckey2(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell2,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
           if(associated(gridp))then
              ff(1:ndim)=ff(1:ndim)+gridp%uold(icell2,2:ndim+1)/max(gridp%uold(icell2,1), r%smallr)*vol2(ind)
           end if
#endif
        end do
        p%vp(ipart,1:ndim)=ff(1:ndim)
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)
     endif
  end do
  call close_cache(s,m%grid_dict)
  if(action_part==action_kick_drift)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
              if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
           endif
        end do
     end do
  end if
  end associate
end subroutine pcs_trace_gas_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cic_kick_drift_dust(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, twotondim
  use hydro_parameters, only: nener
  use pm_parameters
  use pm_commons, only: part_t
  use oct_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  ! Intrinsic functions for drag calculations
  intrinsic :: sinh, cosh
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  real(kind=8),dimension(1:ndim)::x,x_mid,dr,dl,dr2,dl2
  integer,dimension(1:ndim)::ir,il
  integer,dimension(1:ndim)::ir2,il2
  real(kind=8),dimension(1:twotondim)::vol,vol2
  integer,dimension(1:ndim,1:twotondim)::ckey,ckey2
  integer::icell,icell2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,ind,idim,irad
  real(kind=8)::dx_loc,vol_loc
  real(kind=8),dimension(1:ndim)::ff,v_pred
  real(kind=8),dimension(1:ndim)::uu
#ifdef MHD
  real(kind=8),dimension(1:3)::bb
  real(kind=8)::emag
#endif
  real(kind=8)::rho_gas,c_sound,eint,coeff,wdrift2
  real(kind=8)::nu_stop,dens,etot,ekin,erad,cs2,pi=4.0d0*atan(1.0d0)
  real(kind=8),dimension(1:ndim)::what,wdrift ! drift velocity unit vector
  type(oct),pointer :: gridp
  logical :: ok_level
  type(msg_large_realdp)::dummy_large_realdp

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  dx_loc=r%boxlen/2**ilevel
  vol_loc=dx_loc**ndim
  if (p%type/=DUST_TYPE) return
  coeff=9.0d0*pi*r%gamma/128.0d0
  ! Dust hydro+gravity cache
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_large_realdp)/32,&
                     pack=pack_fetch_kick_dust,unpack=unpack_fetch_kick_dust)

  do ipart=p%headp(ilevel),p%tailp(ilevel)
     ! Position in cell units and wrap
     do idim=1,ndim
        x(idim)=p%xp(ipart,idim)/dx_loc
     end do
     do idim=1,ndim
        if(x(idim)<0d0)x(idim)=x(idim)+dble(m%ckey_max(ilevel+1))
        if(x(idim)>=dble(m%ckey_max(ilevel+1)))x(idim)=x(idim)-dble(m%ckey_max(ilevel+1))
     end do

     ! CIC weights/indices at x
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
     vol = cic_weight(dl,dr)



     if(action_part==action_kick_only)then
        p%levelp(ipart)=ilevel
     else if(action_part==action_kick_drift)then
        ! Leapfrog-style sampling: half-step shift only on the first coarse step
        v_pred(1:ndim)=p%vp(ipart,1:ndim)
        if (g%nstep==0) then
           do idim=1,ndim
              x_mid(idim)=x(idim)+0.5d0*g%dtnew(ilevel)*v_pred(idim)/dx_loc
           end do
        else
           do idim=1,ndim
              x_mid(idim)=x(idim)
           end do
        endif
        do idim=1,ndim
           if(x_mid(idim)<0d0)x_mid(idim)=x_mid(idim)+dble(m%ckey_max(ilevel+1))
           if(x_mid(idim)>=dble(m%ckey_max(ilevel+1)))x_mid(idim)=x_mid(idim)-dble(m%ckey_max(ilevel+1))
        end do

        ! CIC weights/indices at x_mid
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
        vol2 = cic_weight(dl2,dr2)

        ! Gather hydro and forces at x_mid
        ff(1:ndim)=0.0
        uu(1:ndim)=0.0
        rho_gas = 0.0
        eint = 0.0
#ifdef MHD
        bb(1:3)=0.0
        emag=0.0d0
#endif
        hash_nbor(0)=ilevel+1
        do ind=1,twotondim
           hash_nbor(1:ndim)=ckey2(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell2,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
           if(associated(gridp))then
#ifdef GRAV
              ff(1:ndim)=ff(1:ndim)+gridp%f(icell2,1:ndim)*vol2(ind)
#endif
              rho_gas = rho_gas + gridp%uold(icell2,1)*vol2(ind)
              uu(1:ndim)=uu(1:ndim)+gridp%uold(icell2,2:ndim+1)/max(gridp%uold(icell2,1), r%smallr)*vol2(ind)
#ifdef MHD
              bb(1:3)=bb(1:3)+0.5d0*(gridp%bold(icell2,1:3)+gridp%bold(icell2,4:6))*vol2(ind)
#endif
              dens = max(dble(gridp%uold(icell2,1)), r%smallr)
              etot = gridp%uold(icell2,5)
              ekin = 0.0d0
              do idim=1,ndim
                 ekin = ekin + 0.5d0 * gridp%uold(icell2,1+idim)**2 / dens
              end do
              erad = 0.0d0
#if NENER>0
              do irad=1,nener
                 erad = erad + gridp%uold(icell2,5+irad)
              end do
#endif
#ifdef MHD
              do idim=1,3
                 emag = emag + 0.125d0*(gridp%bold(icell2,idim)+gridp%bold(icell2,ndim+idim))**2*vol2(ind)
              end do
              eint = eint - emag*vol2(ind)
#endif
              eint = eint + (etot - ekin - erad)*vol2(ind)
           end if
#endif
        end do

        cs2 = r%gamma * (r%gamma-1.0d0) * max(eint, r%smallc**2) / max(rho_gas, r%smallr)
        c_sound = max(sqrt(cs2), r%smallc)
        nu_stop = coeff*c_sound*rho_gas/p%size(ipart)
        wdrift(1:ndim)= v_pred(1:ndim) - uu(1:ndim)

        ! Operator-split: drag(half) -> forces -> drag(half)
        call compute_drag_step(wdrift, c_sound, 0.5d0*g%dtnew(ilevel), nu_stop, coeff, r%analytic_dust_force)
#ifdef GRAV
        wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)
#endif
#ifdef MHD
        if (ndim==3) then
            call compute_lorentz_step(wdrift, bb(1:ndim), g%dtnew(ilevel), p%charge(ipart), r%analytic_dust_force)
        else
            if(g%myid==1 .and. g%nstep==0)then
               write(*,*) 'Warning: Lorentz force not implemented for NDIM != 3; proceeding without it.'
            endif
        endif
#endif
#ifdef GRAV
        wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)
#endif
        call compute_drag_step(wdrift, c_sound, 0.5d0*g%dtnew(ilevel), nu_stop, coeff, r%analytic_dust_force)

        p%vp(ipart,1:ndim)=uu(1:ndim)+wdrift(1:ndim)
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)
     endif
  end do
  call close_cache(s,m%grid_dict)
  
  end associate
end subroutine cic_kick_drift_dust

subroutine tsc_kick_drift_dust(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, threetondim
  use hydro_parameters, only: nener
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
  integer::ipart,icell,icell2,ind,idim,irad
  real(kind=8)::xl,xc,xr
  real(kind=8)::dx_loc
  real(kind=8),dimension(1:ndim)::ff,uu,v_pred,wdrift

  real(kind=8),dimension(1:3)::bb

  real(kind=8)::rho_gas,c_sound,eint,coeff
  real(kind=8)::nu_stop,dens,etot,ekin,erad,emag,cs2,pi
  integer :: ii
  character(LEN=80)::filename,fileloc
  character(LEN=5)::nchar
  type(oct),pointer::gridp
  type(msg_large_realdp)::dummy_large_realdp

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  if (p%type/=DUST_TYPE) return
  dx_loc=r%boxlen/2**ilevel
  pi=4.0d0*atan(1.0d0)
  coeff=9.0d0*pi*r%gamma/128.0d0
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
       hilbert=m%domain, pack_size=storage_size(dummy_large_realdp)/32,&
       pack=pack_fetch_kick_dust,unpack=unpack_fetch_kick_dust)
  do ipart=p%headp(ilevel),p%tailp(ilevel)
     ! particle position in cell units and periodic wrap
     do idim=1,ndim
        x(idim)=p%xp(ipart,idim)/dx_loc
     end do
     do idim=1,ndim
        if(x(idim)<0d0)x(idim)=x(idim)+dble(m%ckey_max(ilevel+1))
        if(x(idim)>=dble(m%ckey_max(ilevel+1)))x(idim)=x(idim)-dble(m%ckey_max(ilevel+1))
     end do

     if(action_part==action_kick_only)then
        p%levelp(ipart)=ilevel

     else if(action_part==action_kick_drift)then
        ! Leapfrog-style sampling: do a half-step advance only on the first step,
        ! then sample at the current position thereafter (we are already time-centered).
        v_pred(1:ndim)=p%vp(ipart,1:ndim)

        if (g%nstep==0) then
           ! First step: predict position at midpoint t + dt/2 for sampling
           do idim=1,ndim
              x_mid(idim)=x(idim)+0.5d0*g%dtnew(ilevel)*v_pred(idim)/dx_loc
           end do
        else
           ! Subsequent steps: positions are already at mid-time, no shift needed
           do idim=1,ndim
              x_mid(idim)=x(idim)
           end do
        endif
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

        ff(1:ndim)=0.0d0
        uu(1:ndim)=0.0d0
        bb(1:ndim)=0.0d0
        rho_gas=0.0d0
        eint=0.0d0
        emag=0.0d0
        hash_nbor(0)=ilevel+1
        do ind=1,threetondim
           hash_nbor(1:ndim)=ckey2(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell2,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
           if(associated(gridp))then
#ifdef GRAV   
              ff(1:ndim)=ff(1:ndim)+gridp%f(icell2,1:ndim)*vol2(ind)
#endif
              rho_gas = rho_gas + gridp%uold(icell2,1)*vol2(ind)
              uu(1:ndim)=uu(1:ndim)+gridp%uold(icell2,2:ndim+1)/max(gridp%uold(icell2,1), r%smallr)*vol2(ind)
#ifdef MHD
              bb(1:3)=bb(1:3)+0.5d0*(gridp%bold(icell2,1:3)+gridp%bold(icell2,4:6))*vol2(ind)
#endif
              dens = max(dble(gridp%uold(icell2,1)), r%smallr)
              etot = gridp%uold(icell2,5)
              ekin = 0.0d0
              do idim=1,ndim
                 ekin = ekin + 0.5d0 * gridp%uold(icell2,1+idim)**2 / dens
              end do
              erad = 0.0d0
#if NENER>0
              do irad=1,nener
                 erad = erad + gridp%uold(icell2,5+irad)
              end do
#endif
#ifdef MHD
              do idim=1,3
                 emag = emag + 0.125d0*(gridp%bold(icell2,idim)+gridp%bold(icell2,ndim+idim))**2*vol2(ind)
              end do
#endif
              eint = eint + (etot - ekin - erad - emag) * vol2(ind)
           end if

        ! Need to add MHD support here
        end do
        cs2 = r%gamma * (r%gamma-1.0d0) * max(eint, r%smallc**2) / max(rho_gas, r%smallr)
        c_sound = max(sqrt(cs2), r%smallc)
        nu_stop = coeff*c_sound*rho_gas/p%size(ipart)
        wdrift(1:ndim)= v_pred(1:ndim) - uu(1:ndim)
        ! This is where we want to split off the different physics
        ! There will also be an if(gyro_pic) gate, since everything must be computed very differently
        ! in the gyro case.

        call compute_drag_step(wdrift, c_sound, 0.5*g%dtnew(ilevel), nu_stop, coeff, r%analytic_dust_force)
        ! Gather ff, and apply to either side of the Lorentz force as a half-step
        wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel) ! External force half-step (includes gravity)
#ifdef MHD 
        if (ndim==3) then
           call compute_lorentz_step(wdrift, bb(1:ndim), g%dtnew(ilevel), p%charge(ipart), r%analytic_dust_force) !Somehow introducing nonphysical oscillations.
        else
           if(g%myid==1 .and. g%nstep==0)then
              write(*,*) 'Warning: Lorentz force not implemented for NDIM != 3; proceeding without it.'
           endif
        endif
#endif
        wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel) ! Second external force half-step
        call compute_drag_step(wdrift, c_sound, 0.5*g%dtnew(ilevel), nu_stop, coeff, r%analytic_dust_force)
        ! Routine will return an intermediate drift velocity.
        p%vp(ipart,1:ndim)=uu(1:ndim)+wdrift(1:ndim)
        ! Leapfrog drift: advance positions using time-centered velocity
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)*g%dtnew(ilevel)
#endif
                ! Trajectory output for selected particles
        if(s%r%ntrajectories>0)then
         do ii=1,s%r%ntrajectories
            if(s%r%trajectories(ii)==p%idp(ipart))then
               call title(g%myid,nchar)
               filename='trajectory.dat'
               fileloc=TRIM(filename)//TRIM(nchar)
               open(25+g%myid,file=fileloc,status='unknown',access='append')
               write(25+g%myid,'(1PE15.7,1X,I12,1X,6(1PE15.7,1X),1PE15.7,1X,1PE15.7)') &
                    g%t, p%idp(ipart), &
                    p%xp(ipart,1),p%xp(ipart,2),p%xp(ipart,3), &
                    p%vp(ipart,1),p%vp(ipart,2),p%vp(ipart,3), &
                    p%charge(ipart), p%size(ipart)
               close(25+g%myid)
               exit
            endif
         end do
      endif
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
end subroutine tsc_kick_drift_dust

subroutine pcs_kick_drift_dust(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, fourtondim
  use hydro_parameters, only: nener
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
  real(kind=8),dimension(1:ndim)::ff,v_pred,uu,what,wdrift
  real(kind=8),dimension(1:3)::bb
  real(kind=8)::emag
  real(kind=8)::rho_gas,c_sound,eint,coeff,wdrift2
  real(kind=8)::nu_stop,dens,etot,ekin,erad,cs2,pi
  integer::irad
  type(oct),pointer::gridp
  type(msg_large_realdp)::dummy_large_realdp

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  if (p%type/=DUST_TYPE) return
  dx_loc=r%boxlen/2**ilevel
  pi=4.0d0*atan(1.0d0)
  coeff=9.0d0*pi*r%gamma/128.0d0
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
       hilbert=m%domain, pack_size=storage_size(dummy_large_realdp)/32,&
       pack=pack_fetch_kick_dust,unpack=unpack_fetch_kick_dust)
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
#ifdef HYDRO
        if(associated(gridp))then
           ff(1:ndim)=ff(1:ndim)+gridp%uold(icell,2:ndim+1)/max(gridp%uold(icell,1), r%smallr)*vol(ind)
        end if
#endif
     end do

    if(action_part==action_kick_only)then
       ! RK2 step 1: stash state; keep velocity as-is for drag scheme
        p%levelp(ipart)=ilevel

    else if(action_part==action_kick_drift)then
       ! RK2 step 2: predict with v^n, sample gas state at midpoint
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
       uu(1:ndim)=0.0
       rho_gas=0.0d0
       eint=0.0d0
       bb(1:3)=0.0d0
       emag=0.0d0
       do ind=1,fourtondim
          hash_nbor(1:ndim)=ckey2(1:ndim,ind)
          call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell2,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
          if(associated(gridp))then
#ifdef GRAV
             ff(1:ndim)=ff(1:ndim)+gridp%f(icell2,1:ndim)*vol2(ind)
#endif
             rho_gas = rho_gas + gridp%uold(icell2,1)*vol2(ind)
             uu(1:ndim)=uu(1:ndim)+gridp%uold(icell2,2:ndim+1)/max(gridp%uold(icell2,1), r%smallr)*vol2(ind)
#ifdef MHD
             bb(1:3)=bb(1:3)+0.5d0*(gridp%bold(icell2,1:3)+gridp%bold(icell2,4:6))*vol2(ind)
#endif
             dens = max(dble(gridp%uold(icell2,1)), r%smallr)
             etot = gridp%uold(icell2,5)
             ekin = 0.0d0
             do idim=1,ndim
                ekin = ekin + 0.5d0 * gridp%uold(icell2,1+idim)**2 / dens
             end do
             erad = 0.0d0
#if NENER>0
             do irad=1,nener
                erad = erad + gridp%uold(icell2,5+irad)
             end do
#endif
#ifdef MHD
             do idim=1,3
                emag = emag + 0.125d0*(gridp%bold(icell2,idim)+gridp%bold(icell2,ndim+idim))**2*vol2(ind)
             end do
             eint = eint - emag*vol2(ind)
#endif
             eint = eint + (etot - ekin - erad) * vol2(ind)
          end if
#endif
       end do
       cs2 = r%gamma * (r%gamma-1.0d0) * max(eint, r%smallc**2) / max(rho_gas, r%smallr)
       c_sound = max(sqrt(cs2), r%smallc)
       nu_stop = coeff*c_sound*rho_gas/p%size(ipart)
       wdrift(1:ndim)= v_pred(1:ndim) - uu(1:ndim)
       ! Operator-split: drag(half) -> forces -> drag(half)
       call compute_drag_step(wdrift, c_sound, 0.5d0*g%dtnew(ilevel), nu_stop, coeff, r%analytic_dust_force)
#ifdef GRAV
       wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)
#endif
#ifdef MHD
       if (ndim==3) then
           call compute_lorentz_step(wdrift, bb(1:ndim), g%dtnew(ilevel), p%charge(ipart), r%analytic_dust_force)
       else
           if(g%myid==1 .and. g%nstep==0)then
              write(*,*) 'Warning: Lorentz force not implemented for NDIM != 3; proceeding without it.'
           endif
       endif
#endif
#ifdef GRAV
       wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*g%dtnew(ilevel)
#endif
       call compute_drag_step(wdrift, c_sound, 0.5d0*g%dtnew(ilevel), nu_stop, coeff, r%analytic_dust_force)
       p%vp(ipart,1:ndim)=uu(1:ndim)+wdrift(1:ndim)
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
end subroutine pcs_kick_drift_dust

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine compute_drag_step(wdrift, c_sound, dt, nu_stop, coeff, analytic_dust_force)
  use amr_parameters, only: ndim
  implicit none
  real(kind=8), intent(inout) :: wdrift(1:ndim)
  real(kind=8), intent(in)    :: c_sound, dt, nu_stop, coeff
  logical, intent(in)         :: analytic_dust_force
  real(kind=8) :: wdrift2
  real(kind=8) :: what(1:ndim)

  wdrift2 = dot_product(wdrift(1:ndim), wdrift(1:ndim))
  if (wdrift2 > 0.0d0) then
     what(1:ndim) = wdrift(1:ndim) / sqrt(wdrift2)
     if (analytic_dust_force) then
        wdrift(1:ndim) = c_sound * what(1:ndim) * &
        & exact_drag(dt, nu_stop, coeff, sqrt(wdrift2)/c_sound)
     else
        wdrift(1:ndim) = c_sound * what(1:ndim) * &
        & fully_implicit_drag(dt, nu_stop, coeff, sqrt(wdrift2)/c_sound)
     end if
  else
     wdrift(1:ndim) = 0.0d0
  end if

end subroutine compute_drag_step

function exact_drag(dt, nu, eta, w0) result(drag_value)
  ! Computes drag coefficient based on the analytic solution for Epstein-Baines drag:
  ! drag(dt) = sqrt(((sinh(nu*dt) + sqrt(1 + eta*w0^2)*cosh(nu*dt)) / 
  !                  (cosh(nu*dt) + sqrt(1 + nu*w0^2)*sinh(nu*dt)))^2 - 1)
  implicit none
  real(kind=8), intent(in) :: dt, nu, eta, w0
  real(kind=8) :: drag_value
  
  real(kind=8) :: sinh_term, cosh_term, sqrt_term, numerator, denominator, ratio
  
  ! Compute hyperbolic functions
  sinh_term = sinh(nu * dt)
  cosh_term = cosh(nu * dt)
  
  ! Compute sqrt(1 + eta*w0^2) - note: using eta as in the original formula
  sqrt_term = sqrt(1.0d0 + eta * w0**2)
  
  ! Compute numerator: sinh(nu*dt) + sqrt(1 + eta*w0^2)*cosh(nu*dt)
  numerator = sinh_term + sqrt_term * cosh_term
  
  ! Compute denominator: cosh(nu*t) + sqrt(1 + eta*w0^2)*sinh(nu*t)
  denominator = cosh_term + sqrt(1.0d0 + eta * w0**2) * sinh_term
  
  ! Avoid division by zero
  if (abs(denominator) < 1.0d-15) then
     drag_value = 0.0d0
     return
  end if
  
  ! Compute ratio and final result
  ratio = numerator / denominator
  drag_value = sqrt(max(0.0d0, ratio**2 - 1.0d0))/sqrt(eta)
  
end function exact_drag

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
function fully_implicit_drag(dt, nu, eta, w0) result(drag_value)
  ! A second order fully implicit drag update
  implicit none
  real(kind=8), intent(in) :: dt, nu, eta, w0
  real(kind=8) :: drag_value, nu_stop
  
  real(kind=8) :: sqrt_term, ratio
  
  
  ! Compute sqrt(1 + eta*w0^2) w0^2 is actually w0^2/c_sound^2
  sqrt_term = sqrt(1.0d0 + eta * w0**2)
  nu_stop = nu * sqrt_term

  drag_value = w0/(1+nu_stop*dt+0.5d0*nu*dt**2) 
  ! Not a typo. Comes from a series expansion of the analytic solution in dt

  
end function fully_implicit_drag

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine compute_lorentz(driftvel, bfield, dt, charge_parameter)
  use amr_parameters, only: ndim
  implicit none
  real(kind=8), dimension(1:ndim), intent(inout) :: driftvel
  real(kind=8), dimension(1:ndim), intent(in)    :: bfield
  real(kind=8), intent(in)                       :: dt
  real(kind=8), intent(in)                       :: charge_parameter
  real(kind=8) :: det, bsquared, dteff
  real(kind=8) :: v1_new, v2_new, v3_new
  real(kind=8), dimension(1:3,1:3) :: matrix
  
  dteff = -1.0d0 * dt * charge_parameter ! Accidentally flipped sign in the original code.
  bsquared = dot_product(bfield(1:ndim), bfield(1:ndim))

  det = 1 + 0.25d0 * bsquared * dteff**2
  
  ! Matrix components from Mathematica
  ! Inverse[(e - 0.5*dteff*m)].(e+0.5*dteff*m).v^n = v^(n+1)
  ! Row 1
  matrix(1,1) = 1.0d0 + 0.25d0 * (bfield(1)**2 - bfield(2)**2 - bfield(3)**2) * dteff**2
  matrix(1,2) = 0.5d0 * dteff * (2.0d0 * bfield(3) + bfield(1) * bfield(2) * dteff)
  matrix(1,3) = 0.5d0 * dteff * (-2.0d0 * bfield(2) + bfield(1) * bfield(3) * dteff)
  
  ! Row 2
  matrix(2,1) = 0.5d0 * dteff * (-2.0d0 * bfield(3) + bfield(1) * bfield(2) * dteff)
  matrix(2,2) = 1.0d0 - 0.25d0 * (bfield(1)**2 - bfield(2)**2 + bfield(3)**2) * dteff**2
  matrix(2,3) = 0.5d0 * dteff * (2.0d0 * bfield(1) + bfield(2) * bfield(3) * dteff)
  
  ! Row 3
  matrix(3,1) = 0.5d0 * dteff * (2.0d0 * bfield(2) + bfield(1) * bfield(3) * dteff)
  matrix(3,2) = 0.5d0 * dteff * (-2.0d0 * bfield(1) + bfield(2) * bfield(3) * dteff)
  matrix(3,3) = 1.0d0 - 0.25d0 * (bfield(1)**2 + bfield(2)**2 - bfield(3)**2) * dteff**2
  
  ! Apply matrix/det to driftvel (hard-coded for efficiency)
  if (ndim == 3) then
     v1_new = (matrix(1,1)*driftvel(1) + matrix(1,2)*driftvel(2) + matrix(1,3)*driftvel(3)) / det
     v2_new = (matrix(2,1)*driftvel(1) + matrix(2,2)*driftvel(2) + matrix(2,3)*driftvel(3)) / det
     v3_new = (matrix(3,1)*driftvel(1) + matrix(3,2)*driftvel(2) + matrix(3,3)*driftvel(3)) / det
     driftvel(1) = v1_new
     driftvel(2) = v2_new
     driftvel(3) = v3_new
  end if
  
end subroutine compute_lorentz

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine compute_lorentz_step(driftvel, bfield, dt, charge_parameter, analytic_dust_force)
  use amr_parameters, only: ndim
  implicit none
  real(kind=8), dimension(1:ndim), intent(inout) :: driftvel
  real(kind=8), dimension(1:ndim), intent(in)    :: bfield
  real(kind=8), intent(in)                       :: dt
  real(kind=8), intent(in)                       :: charge_parameter
  logical, intent(in)                            :: analytic_dust_force
#ifdef MHD
  if (ndim/=3) return
  if (.not. analytic_dust_force) then
     call compute_lorentz(driftvel, bfield, dt, charge_parameter)
  else
     call compute_lorentz_analytic(driftvel, bfield, dt, charge_parameter)
  end if
#else
  ! No MHD compiled; nothing to do
#endif
end subroutine compute_lorentz_step
!#########################################################################
!#########################################################################
subroutine compute_lorentz_analytic(driftvel, bfield, dt, charge_parameter)
  use amr_parameters, only: ndim
  implicit none
  real(kind=8), dimension(1:ndim), intent(inout) :: driftvel
  real(kind=8), dimension(1:ndim), intent(in)    :: bfield
  real(kind=8), intent(in)                       :: dt
  real(kind=8), intent(in)                       :: charge_parameter
  real(kind=8) :: dteff, bsquared, bnorm, theta, costh, sinth, vdotb, t2
  real(kind=8), dimension(1:ndim) :: bhat, v, kxv


  dteff = dt * charge_parameter
  if (dteff == 0.0d0) return

  bsquared = dot_product(bfield(1:ndim), bfield(1:ndim))
  if (bsquared <= 0.0d0) return

  bnorm = sqrt(bsquared)
  bhat(1:ndim) = bfield(1:ndim) / bnorm
  theta = dteff * bnorm
  ! It's worth noting that our gyro_factor limits theta to be less than 0.63
  t2 = theta*theta
  if (abs(theta) < 1.0d-1) then
     ! Horner-form small-angle approximations for better performance/accuracy
     costh = 1.0d0 + t2*(-0.5d0 + t2*( 1.0d0/24.0d0 + t2*( -1.0d0/720.0d0 )))
     sinth = theta * (1.0d0 + t2*( -1.0d0/6.0d0 + t2*( 1.0d0/120.0d0 + t2*( -1.0d0/5040.0d0 ))))
  else
     costh = cos(theta)
     sinth = sin(theta)
  end if

  v(1:ndim) = driftvel(1:ndim)

  if (ndim == 3) then
     ! Rodrigues' rotation formula: rotate v around unit axis bhat by angle theta
     kxv(1) = bhat(2)*v(3) - bhat(3)*v(2)
     kxv(2) = bhat(3)*v(1) - bhat(1)*v(3)
     kxv(3) = bhat(1)*v(2) - bhat(2)*v(1)
     vdotb  = bhat(1)*v(1) + bhat(2)*v(2) + bhat(3)*v(3)

     driftvel(1) = v(1)*costh + kxv(1)*sinth + bhat(1)*vdotb*(1.0d0 - costh)
     driftvel(2) = v(2)*costh + kxv(2)*sinth + bhat(2)*vdotb*(1.0d0 - costh)
     driftvel(3) = v(3)*costh + kxv(3)*sinth + bhat(3)*vdotb*(1.0d0 - costh)
  end if

end subroutine compute_lorentz_analytic

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module move_fine_module
