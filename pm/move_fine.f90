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
           call mc_trace_gas_part(pst%s,pst%s%trac,ilevel,action_part) ! Classical Monte Carlo (may not work with AMR)
        elseif(pst%s%r%trac_interpolation_scheme==1)then
           call cic_trace_gas_part(pst%s,pst%s%trac,ilevel,action_part)
        elseif(pst%s%r%trac_interpolation_scheme==2)then
           call tsc_trace_gas_part(pst%s,pst%s%trac,ilevel,action_part)
        elseif(pst%s%r%trac_interpolation_scheme==3)then
           call pcs_trace_gas_part(pst%s,pst%s%trac,ilevel,action_part)
        elseif(pst%s%r%trac_interpolation_scheme==4)then
           call cic_trace_gas_part_ito_mc(pst%s,pst%s%trac,ilevel,action_part) ! Ito formulation of the flux-based Monte Carlo tracer
        elseif(pst%s%r%trac_interpolation_scheme==5)then
           call tsc_trace_gas_part_ito_mc(pst%s,pst%s%trac,ilevel,action_part) ! Ito MC tracer with TSC
        elseif(pst%s%r%trac_interpolation_scheme==6)then
           call cic_trace_gas_part_sgs_turb(pst%s,pst%s%trac,ilevel,action_part) ! SGS turbulent diffusion tracer with CIC
        elseif(pst%s%r%trac_interpolation_scheme==7)then
           call tsc_trace_gas_part_sgs_turb(pst%s,pst%s%trac,ilevel,action_part) ! SGS turbulent diffusion tracer with TSC
        elseif(pst%s%r%trac_interpolation_scheme==8)then
           call trace_gas_part_trivial(pst%s,pst%s%trac,ilevel,action_part) ! Trivial tracer: no interpolation, fixed velocity, random diffusion
        endif
     endif
     if(pst%s%r%dust)then
        if(pst%s%r%dust_force_interpolation_scheme==1)then
           call cic_kick_drift_dust(pst%s,pst%s%dust,ilevel,action_part)
        elseif(pst%s%r%dust_force_interpolation_scheme==2)then
           call tsc_kick_drift_dust(pst%s,pst%s%dust,ilevel,action_part)
        elseif(pst%s%r%dust_force_interpolation_scheme==3)then
           call pcs_kick_drift_dust(pst%s,pst%s%dust,ilevel,action_part)
        elseif(pst%s%r%dust_force_interpolation_scheme==4)then
           call cic_kick_drift_dust_ito_mc(pst%s,pst%s%dust,ilevel,action_part) ! Asymptotically approaches Ito MC tracer limit
        elseif(pst%s%r%dust_force_interpolation_scheme==5)then
           call tsc_kick_drift_dust_ito_mc(pst%s,pst%s%dust,ilevel,action_part) ! Asymptotically approaches Ito MC tracer limit
        elseif(pst%s%r%dust_force_interpolation_scheme==6)then
           call cic_kick_drift_dust_guiding_center(pst%s,pst%s%dust,ilevel,action_part) ! CIC guiding center
        elseif(pst%s%r%dust_force_interpolation_scheme==7)then
           call tsc_kick_drift_dust_guiding_center(pst%s,pst%s%dust,ilevel,action_part) ! TSC guiding center
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

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

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
        call gravana(r,g,xana,fana,dx_loc,1)
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
        ! Store potential
        if(allocated(p%phip))then
           p%phip(ipart)=0.0
           do ind=1,twotondim
              p%phip(ipart)=p%phip(ipart)+gridp(ind)%p%phi(icell(ind))*vol(ind)
           end do
        endif
        ! Store old force
        if(allocated(p%fp))then
           p%fp(ipart,1:ndim)=ff(1:ndim)
        endif
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
  use cache_commons, only: msg_hydro_mflux
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_hydro_mflux)::msg

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        msg%realdp_hydro(ind,ivar)=grid%uold(ind,ivar)
     end do
  end do
  do ind=1,twotondim
     msg%realdp_mflux(ind,1:2*ndim+1)=grid%mflux(ind,1:2*ndim+1)
     msg%realdp_upwind_rho(ind,1:2*ndim)=grid%upwind_rho(ind,1:2*ndim)
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
  use cache_commons, only: msg_hydro_mflux
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key
  integer::ind,ivar
  type(msg_hydro_mflux)::msg

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
     grid%mflux(ind,1:2*ndim+1)=msg%realdp_mflux(ind,1:2*ndim+1)
     grid%upwind_rho(ind,1:2*ndim)=msg%realdp_upwind_rho(ind,1:2*ndim)
  end do
#endif

end subroutine unpack_fetch_kick_trac
!#########################################################################
!#########################################################################
subroutine pack_fetch_kick_dust(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim, ndim
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
     msg%realdp_mflux(ind,1:2*ndim+1)=grid%mflux(ind,1:2*ndim+1)
     msg%realdp_upwind_rho(ind,1:2*ndim)=grid%upwind_rho(ind,1:2*ndim)
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
     grid%mflux(ind,1:2*ndim+1)=msg%realdp_mflux(ind,1:2*ndim+1)
     grid%upwind_rho(ind,1:2*ndim)=msg%realdp_upwind_rho(ind,1:2*ndim)
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
  type(msg_hydro_mflux)::dummy_nvar_realdp
  type(RngStream),external::RngStream_CreateStream
  real(kind=8),external::RngStream_RandUni
  integer(kind=8)::stream_skip
  external :: RngStream_SetPackageSeed, RngStream_AdvanceState, gaussdev
  integer :: ipart,idim
  integer :: ii
  character(LEN=80)::filename,fileloc
  character(LEN=5)::nchar

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  if (p%type/=TRAC_TYPE) return

  dx_loc=r%boxlen/2**ilevel
  dt_level=g%dtnew(ilevel)
  use_sgs = .false.

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

     ! Trajectory output for selected particles
     if(s%r%ntrajectories>0)then
        do ii=1,s%r%ntrajectories
           if(s%r%trajectories(ii)==p%idp(ipart))then
              call title(g%myid,nchar)
              filename='trajectory.dat'
              fileloc=TRIM(filename)//TRIM(nchar)
              open(25+g%myid,file=fileloc,status='unknown',access='append')
              write(25+g%myid,*) g%t, p%idp(ipart), &
                   (p%xp(ipart,idim), idim=1,ndim), &
                   (p%vp(ipart,idim), idim=1,ndim)
              close(25+g%myid)
              exit
           endif
        end do
     endif

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

subroutine cic_trace_gas_part_sgs_turb(s,p,ilevel,action_part)
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
  real(kind=8),dimension(1:ndim)::x,x_mid,v_pred,u_eff,u_mid,disp,xi,momentum
  real(kind=8),dimension(1:ndim)::grad_at_part
  real(kind=8),dimension(1:twotondim)::vol,kappa_cells
  integer,dimension(1:ndim)::il,ir
  integer,dimension(1:ndim,1:twotondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  real(kind=8),dimension(1:ndim,1:2)::w1d,dw1d
  real(kind=8)::dx_loc,dt_level,kappa_mid,noise_amp,rho
  logical :: use_sgs
  type(oct),pointer::gridp
  type(msg_hydro_mflux)::dummy_nvar_realdp
  type(RngStream),external::RngStream_CreateStream
  real(kind=8),external::RngStream_RandUni
  integer(kind=8)::stream_skip
  external :: RngStream_SetPackageSeed, RngStream_AdvanceState, gaussdev
  integer :: ipart,idim,ind,icell
  integer :: ii
  character(LEN=80)::filename,fileloc
  character(LEN=5)::nchar

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  if (p%type/=TRAC_TYPE) return

  dx_loc=r%boxlen/2**ilevel
  dt_level=g%dtnew(ilevel)
  use_sgs = r%sgs_turb .and. (r%iturb>0)

  ! OPTIMIZATION: Handle kick-only pass immediately without cache operations
  if(action_part==action_kick_only)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        p%levelp(ipart)=ilevel
     end do
     return
  endif

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

     ! Predictor velocity: use stored velocity if available, otherwise gather from grid
     if (g%nstep>0) then
        v_pred(1:ndim)=p%vp(ipart,1:ndim)
     else
        call gather_cic_state(s,x,ilevel,dx_loc,.false.,u_eff,kappa_mid)
        v_pred(1:ndim)=u_eff(1:ndim)
     endif

     do idim=1,ndim
        x_mid(idim)=x(idim)+0.5d0*dt_level*v_pred(idim)/dx_loc
     end do
     call wrap_cell_coords(s,x_mid,ilevel+1)

     ! Build CIC weights, derivatives, indices, and volume weights ONCE for x_mid
     call cic_weights_and_derivs(x_mid, w1d, dw1d, il, ir, vol)
     do idim=1,ndim
        if(r%periodic(idim))then
           if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=m%box_ckey_max(idim,ilevel+1)-1
           if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=m%box_ckey_min(idim,ilevel+1)
        endif
     end do
     ckey = cic_index(il,ir)

     ! Single pass over cells: gather velocity, density, and per-cell kappa
     momentum(1:ndim)=0.d0
     rho=0.d0
     kappa_mid=0.d0
     kappa_cells(1:twotondim)=0.d0
     hash_nbor(0)=ilevel+1
     do ind=1,twotondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
        if(associated(gridp))then
           momentum(1:ndim)=momentum(1:ndim)+gridp%uold(icell,2:ndim+1)*vol(ind)
           rho=rho+gridp%uold(icell,1)*vol(ind)
           if(use_sgs)then
              kappa_cells(ind)=tracer_cell_kappa(gridp%uold(icell,1),gridp%uold(icell,r%iturb),dx_loc,r%smallr)
              kappa_mid=kappa_mid+kappa_cells(ind)*vol(ind)
           end if
        end if
#endif
     end do
     if(rho>r%smallr)then
        u_mid(1:ndim)=momentum(1:ndim)/rho
     else
        u_mid(1:ndim)=0.d0
     end if

     ! Compute gradient of kappa using the already-computed weights and per-cell values
     grad_at_part(1:ndim)=0.d0
     if(use_sgs)then
        call compute_gradient_cic_scalar(w1d, dw1d, kappa_cells, grad_at_part)
     end if

     ! Ito drift: u + grad(kappa)
     disp(1:ndim)=(u_mid(1:ndim) + grad_at_part(1:ndim)/dx_loc)*dt_level

     if(use_sgs .and. kappa_mid>0.0d0)then
        call sample_tracer_gaussian(xi)
        noise_amp = sqrt(2.0d0*kappa_mid*dt_level)
        disp(1:ndim)=disp(1:ndim)+noise_amp*xi(1:ndim)
     end if

     p%vp(ipart,1:ndim)=u_mid(1:ndim)
     p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+disp(1:ndim)

     ! Trajectory output for selected particles
     if(s%r%ntrajectories>0)then
        do ii=1,s%r%ntrajectories
           if(s%r%trajectories(ii)==p%idp(ipart))then
              call title(g%myid,nchar)
              filename='trajectory.dat'
              fileloc=TRIM(filename)//TRIM(nchar)
              open(25+g%myid,file=fileloc,status='unknown',access='append')
              write(25+g%myid,*) g%t, p%idp(ipart), &
                   (p%xp(ipart,idim), idim=1,ndim), &
                   (p%vp(ipart,idim), idim=1,ndim)
              close(25+g%myid)
              exit
           endif
        end do
     endif

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
end subroutine cic_trace_gas_part_sgs_turb

subroutine tsc_trace_gas_part_sgs_turb(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, threetondim
  use pm_parameters
  use pm_commons, only: part_t
  use oct_commons, only: oct
  use ramses_commons, only: ramses_t
  use rng
  use nbors_utils
  use cache_commons
  use cache
  use rho_fine_module, only: tsc_index
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  real(kind=8),dimension(1:ndim)::x,x_mid,v_pred,u_eff,u_mid,disp,xi,momentum
  real(kind=8),dimension(1:ndim)::grad_at_part
  real(kind=8),dimension(1:threetondim)::vol,kappa_cells
  integer,dimension(1:ndim)::il,ic,ir
  integer,dimension(1:ndim,1:threetondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  real(kind=8),dimension(1:ndim,1:3)::w1d,dw1d
  real(kind=8)::dx_loc,dt_level,kappa_mid,noise_amp,rho
  logical :: use_sgs
  type(oct),pointer::gridp
  type(msg_hydro_mflux)::dummy_nvar_realdp
  type(RngStream),external::RngStream_CreateStream
  real(kind=8),external::RngStream_RandUni
  integer(kind=8)::stream_skip
  external :: RngStream_SetPackageSeed, RngStream_AdvanceState, gaussdev
  integer :: ipart,idim,ind,icell

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  if (p%type/=TRAC_TYPE) return

  dx_loc=r%boxlen/2**ilevel
  dt_level=g%dtnew(ilevel)
  use_sgs = r%sgs_turb .and. (r%iturb>0)

  ! OPTIMIZATION: Handle kick-only pass immediately without cache operations
  if(action_part==action_kick_only)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        p%levelp(ipart)=ilevel
     end do
     return
  endif

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

     ! Predictor velocity: use stored velocity if available, otherwise gather from grid
     if (g%nstep>0) then
        v_pred(1:ndim)=p%vp(ipart,1:ndim)
     else
        call gather_tsc_state(s,x,ilevel,dx_loc,.false.,u_eff,kappa_mid)
        v_pred(1:ndim)=u_eff(1:ndim)
     endif

     do idim=1,ndim
        x_mid(idim)=x(idim)+0.5d0*dt_level*v_pred(idim)/dx_loc
     end do
     call wrap_cell_coords(s,x_mid,ilevel+1)

     ! Build TSC weights, derivatives, indices, and volume weights ONCE for x_mid
     call tsc_weights_and_derivs(x_mid, w1d, dw1d, il, ic, ir, vol)
     do idim=1,ndim
        if(r%periodic(idim))then
           if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=il(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
           if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=ir(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
        endif
     end do
     ckey = tsc_index(il,ic,ir)

     ! Single pass over cells: gather velocity, density, and per-cell kappa
     momentum(1:ndim)=0.d0
     rho=0.d0
     kappa_mid=0.d0
     kappa_cells(1:threetondim)=0.d0
     hash_nbor(0)=ilevel+1
     do ind=1,threetondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
        if(associated(gridp))then
           momentum(1:ndim)=momentum(1:ndim)+gridp%uold(icell,2:ndim+1)*vol(ind)
           rho=rho+gridp%uold(icell,1)*vol(ind)
           if(use_sgs)then
              kappa_cells(ind)=tracer_cell_kappa(gridp%uold(icell,1),gridp%uold(icell,r%iturb),dx_loc,r%smallr)
              kappa_mid=kappa_mid+kappa_cells(ind)*vol(ind)
           end if
        end if
#endif
     end do
     if(rho>r%smallr)then
        u_mid(1:ndim)=momentum(1:ndim)/rho
     else
        u_mid(1:ndim)=0.d0
     end if

     ! Compute gradient of kappa using the already-computed weights and per-cell values
     grad_at_part(1:ndim)=0.d0
     if(use_sgs)then
        call compute_gradient_tsc_scalar(w1d, dw1d, kappa_cells, grad_at_part)
     end if

     ! Ito drift: u + grad(kappa)
     disp(1:ndim)=(u_mid(1:ndim) + grad_at_part(1:ndim)/dx_loc)*dt_level

     if(use_sgs .and. kappa_mid>0.0d0)then
        call sample_tracer_gaussian(xi)
        noise_amp = sqrt(2.0d0*kappa_mid*dt_level)
        disp(1:ndim)=disp(1:ndim)+noise_amp*xi(1:ndim)
     end if

     p%vp(ipart,1:ndim)=u_mid(1:ndim)
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
end subroutine tsc_trace_gas_part_sgs_turb

subroutine cic_trace_gas_part_ito_mc(s,p,ilevel,action_part)
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
  real(kind=8),dimension(1:ndim)::x,disp,xi,u_eff,kappa_num
  real(kind=8),dimension(1:ndim)::grad_at_part
  integer,dimension(1:ndim)::il,ir
  real(kind=8),dimension(1:twotondim)::vol,phi_slice
  integer,dimension(1:ndim,1:twotondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  real(kind=8),dimension(1:ndim,1:twotondim)::u_cells,kappa_num_cells,grad_phi_cells,skew_cells,kurt_cells
  real(kind=8),dimension(1:ndim,1:twotondim)::fluxL_cells,fluxR_cells
  real(kind=8),dimension(1:ndim)::skewness_eff,kurtosis_eff
  real(kind=8),dimension(1:ndim,1:2)::w1d,dw1d
  type(oct),pointer::gridp
  integer :: ipart,ind,idim,icell,k
  integer :: ii
  character(LEN=80)::filename,fileloc
  character(LEN=5)::nchar
  real(kind=8)::dx_loc,dt_level,rho,denom,fluxL,fluxR,jr,jl,noise_amp,cfl_eff,one_minus_cfl,pr,pl
  type(msg_hydro_mflux)::dummy_nvar_realdp
  type(RngStream),external::RngStream_CreateStream
  real(kind=8),external::RngStream_RandUni
  integer(kind=8)::stream_skip
  external :: RngStream_SetPackageSeed, RngStream_AdvanceState, gaussdev

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  if (p%type/=TRAC_TYPE) return

  dx_loc=r%boxlen/2**ilevel

  ! OPTIMIZATION: Handle kick-only pass immediately without cache operations or dt_level
  if(action_part==action_kick_only)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        p%levelp(ipart)=ilevel
     end do
     return
  endif

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

     ! Build 1D CIC weights, indices, and 3D volume weights
     call cic_weights_and_derivs(x, w1d, dw1d, il, ir, vol)

     ! Periodic wrap of indices for 2x2x2 stencil
     do idim=1,ndim
        if(r%periodic(idim))then
           if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=m%box_ckey_max(idim,ilevel+1)-1
           if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=m%box_ckey_min(idim,ilevel+1)
        endif
     enddo

     ! Build neighbor index list (twotondim=8)
     ckey = cic_index(il,ir)

     u_cells=0.d0
     kappa_num_cells=0.d0
     skew_cells=0.d0

     hash_nbor(0)=ilevel+1
     do ind=1,twotondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
        if(associated(gridp))then
           rho=gridp%mflux(icell,1)
           denom=max(rho,r%smallr)
           do idim=1,ndim
              ! The trick might be to use a cell-centered average of the face-centered downwind fluxes.
              ! So we use the density corresponding to the upwind density for a given face flux.
              ! mflux stores time-integrated flux ~ (dt/dx)*F/rho; recover physical flux F with factor rho*dx_loc/dt_level
              fluxL=gridp%mflux(icell,1+idim     )
              fluxR=gridp%mflux(icell,1+idim+ndim)
              jr=max(fluxR,0.d0)*dx_loc/dt_level
              jl=max(-fluxL,0.d0)*dx_loc/dt_level
              u_cells(idim,ind)=(jr-jl)/denom
              !u_cells(idim,ind)=u_cells(idim,ind) !Amplify the difference.
              !u_cells(idim,ind)=0.5d0*(fluxR+fluxL)/denom!(jr-jl)/denom
              !u_cells(idim,ind)=gridp%uold(icell,1+idim)/gridp%uold(icell,1)
              !cfl_eff=(jr+jl)/denom*dt_level/dx_loc
              cfl_eff=(jr+jl)/denom*dt_level/dx_loc
              one_minus_cfl = max(0.d0,(1.d0-cfl_eff))
              !kappa_num_cells(idim,ind)=0.5d0*(jr+jl)/denom*dx_loc*one_minus_cfl 
              kappa_num_cells(idim,ind)=0.5d0* cfl_eff * one_minus_cfl * dx_loc**2.d0 / dt_level

              ! Skewness of MC discrete kernel for two-piece uniform sampling
              pr=max(fluxR,0.d0)/denom
              pl=max(-fluxL,0.d0)/denom
              skew_cells(idim,ind)=mc_kernel_skewness(pr,pl)
              kurt_cells(idim,ind)=mc_kernel_kurtosis(pr,pl)

              ! Also require that this is not the first timestep.
              if ((idim==3) .and. (g%t.ge.1.d0) .and. ind==1 .and. (mod(p%idp(ipart), 64**3) == 1)) then
                 write(*,*) 'particle id', p%idp(ipart)
                 write(*,*) 'u_cells(idim,ind)', u_cells(idim,ind)
                 write(*,*) 'kappa_num_cells(idim,ind)', kappa_num_cells(idim,ind)
                 write(*,*) 'rho', denom
                 write(*,*) 'fluxR', fluxR
                 write(*,*) 'fluxL', fluxL
                 write(*,*) 'jr', jr
                 write(*,*) 'jl', jl
                 write(*,*) 'dt_level', dt_level
                 write(*,*) 'dx_loc', dx_loc
                 write(*,*) 'cfl_eff', cfl_eff
                 write(*,*) 'one_minus_cfl', one_minus_cfl
                 write(*,*) 'skew_cells', skew_cells(idim,ind)
                 write(*,*) 'kurt_cells', kurt_cells(idim,ind)
              endif
              ! When we get back, we can examine how using cell-centered values (below)
              ! impacts the results. We can also examine the impact of including the gradient
              ! term. We should also crank up the particle count, and examine the impact of 
              ! changing how we define one_minus_cfl.
              ! Average upwind values algorithm.
            !   jr = fluxR/gridp%upwind_rho(icell,1+idim+ndim)
            !   jl = fluxL/gridp%upwind_rho(icell,1+idim)
            !   u_cells(idim,ind)=0.5d0*(jr+jl)
            !   cfl_eff=(jr+jl)/denom*dt_level/dx_loc
            !   one_minus_cfl = max(0.d0,(1.d0-cfl_eff))
            !   cfl_eff=abs(jr)*dt_level/dx_loc
            !   one_minus_cfl = max(0.d0,(1.d0-cfl_eff))
            !   kappa_num_cells(idim,ind)=0.25d0*jr*dx_loc*one_minus_cfl
            !   cfl_eff=abs(jl)*dt_level/dx_loc
            !   one_minus_cfl = max(0.d0,(1.d0-cfl_eff))
            !   kappa_num_cells(idim,ind)=kappa_num_cells(idim,ind)+0.25d0*jl*dx_loc*one_minus_cfl
           end do
        end if
#endif
     end do

     u_eff=0.d0
     kappa_num=0.d0
     skewness_eff=0.d0
     kurtosis_eff=0.d0
     do ind=1,twotondim
        do idim=1,ndim
           u_eff(idim)=u_eff(idim)+u_cells(idim,ind)*vol(ind)
           kappa_num(idim)=kappa_num(idim)+kappa_num_cells(idim,ind)*vol(ind)
           skewness_eff(idim)=skewness_eff(idim)+skew_cells(idim,ind)*vol(ind)
           kurtosis_eff(idim)=kurtosis_eff(idim)+kurt_cells(idim,ind)*vol(ind)
        end do
     end do

     ! Reintroduce Ito drift: add grad(kappa_k) along each axis
     ! (i.e., d(kappa_x)/dx, d(kappa_y)/dy, d(kappa_z)/dz)
     ! We need even higher particles per cell. It is unclear to me which is correct.
     ! It might somehow need (1/2)*grad(kappa) to be present. 
     ! This depends on what sort of algorithm we are using (Stratonovich vs Ito).
   !   grad_at_part(1:ndim)=0.d0
   !   !call compute_gradient_cic(w1d, dw1d, kappa_num_cells, grad_at_part)
   !   do k=1,ndim
   !      u_eff(k) = u_eff(k) !+ grad_at_part(k)/dx_loc
   !   end do

     !call sample_tracer_gaussian(xi)
     call sample_tracer_uniform(xi)
     !call sample_tracer_two_piece_uniform(xi,skewness_eff)
     do idim=1,ndim
        disp(idim)=u_eff(idim)*dt_level
        noise_amp = sqrt(max(0.d0,2.d0*kappa_num(idim)*dt_level))
        disp(idim)=disp(idim)+noise_amp*xi(idim)
     end do

     p%vp(ipart,1:ndim)=u_eff(1:ndim)
     p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+disp(1:ndim)

     p%vp(ipart,1)=kappa_num(3)

     ! Trajectory output for selected particles
     if(s%r%ntrajectories>0)then
        do ii=1,s%r%ntrajectories
           if(s%r%trajectories(ii)==p%idp(ipart))then
              call title(g%myid,nchar)
              filename='trajectory.dat'
              fileloc=TRIM(filename)//TRIM(nchar)
              open(25+g%myid,file=fileloc,status='unknown',access='append')
              write(25+g%myid,*) g%t, p%idp(ipart), &
                   (p%xp(ipart,idim), idim=1,ndim), &
                   (p%vp(ipart,idim), idim=1,ndim)
              close(25+g%myid)
              exit
           endif
        end do
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
end subroutine cic_trace_gas_part_ito_mc

! subroutine cic_trace_gas_part_slope_limit(s,p,ilevel,action_part)
! end subroutine cic_trace_gas_part_slope_limit

subroutine tsc_trace_gas_part_ito_mc(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, threetondim
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
  real(kind=8),dimension(1:ndim)::x,disp,xi,u_eff,kappa_num
  real(kind=8),dimension(1:ndim)::grad_at_part
  integer,dimension(1:ndim)::il,ic,ir
  real(kind=8),dimension(1:threetondim)::vol
  integer,dimension(1:ndim,1:threetondim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  real(kind=8),dimension(1:ndim,1:threetondim)::u_cells,kappa_num_cells
  real(kind=8),dimension(1:ndim,1:3)::w1d,dw1d
  type(oct),pointer::gridp
  integer :: ipart,ind,idim,icell
  real(kind=8)::dx_loc,dt_level,dx_over_dt,dt_over_dx,rho,denom,fluxL,fluxR,noise_amp,cfl_eff,one_minus_cfl
  real(kind=8)::jr,jl
  type(msg_hydro_mflux)::dummy_nvar_realdp
  type(RngStream),external::RngStream_CreateStream
  real(kind=8),external::RngStream_RandUni
  integer(kind=8)::stream_skip
  external :: RngStream_SetPackageSeed, RngStream_AdvanceState, gaussdev

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  if (p%type/=TRAC_TYPE) return

  dx_loc=r%boxlen/2**ilevel
  dt_level=g%dtnew(ilevel)
  dx_over_dt=dx_loc/dt_level
  dt_over_dx=dt_level/dx_loc

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

     ! Build 1D TSC weights, derivatives, indices, and 3D volume weights
     call tsc_weights_and_derivs(x, w1d, dw1d, il, ic, ir, vol)

     ! Periodic wrap of indices for 3x3x3 stencil
     do idim=1,ndim
        if(r%periodic(idim))then
           if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=il(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
           if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=ir(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
        endif
     enddo

     ! Build neighbor index list (threetondim=27)
     ckey = tsc_index(il,ic,ir)

     ! Single merged loop over 27 neighbors
     u_cells=0.d0
     kappa_num_cells=0.d0
     hash_nbor(0)=ilevel+1
     do ind=1,threetondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
        if(associated(gridp))then
           rho=gridp%mflux(icell,1)
           denom=max(rho,r%smallr)
           do idim=1,ndim
              fluxL=gridp%mflux(icell,1+idim)*dx_over_dt
              fluxR=gridp%mflux(icell,1+idim+ndim)*dx_over_dt
              jr=max(fluxR,0.d0)
              jl=max(-fluxL,0.d0)
              u_cells(idim,ind)=(jr-jl)/denom
              cfl_eff=(jr+jl)/denom*dt_level/dx_loc
              one_minus_cfl = max(0.d0,(1.d0-cfl_eff))
              kappa_num_cells(idim,ind)=0.5d0*(jr+jl)*dx_loc/denom*one_minus_cfl
           end do
        end if
#endif
     end do

     u_eff=0.d0
     kappa_num=0.d0
     do ind=1,threetondim
        do idim=1,ndim
           u_eff(idim)=u_eff(idim)+u_cells(idim,ind)*vol(ind)
           kappa_num(idim)=kappa_num(idim)+kappa_num_cells(idim,ind)*vol(ind)
        end do
     end do

     ! Compute gradient using TSC subroutine
     !call compute_gradient_tsc(w1d, dw1d, kappa_num_cells, grad_at_part)
     ! u_eff = u_eff + grad_at_part  ! Uncomment to add noise-induced drift

     if(action_part==action_kick_only)then
        p%vp(ipart,1:ndim)=u_eff(1:ndim)
        p%levelp(ipart)=ilevel
        cycle
     endif

     call sample_tracer_uniform(xi)
     do idim=1,ndim
        disp(idim)=u_eff(idim)*dt_level
        noise_amp = sqrt(max(0.d0,2.d0*kappa_num(idim)*dt_level))
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
end subroutine tsc_trace_gas_part_ito_mc

subroutine mc_trace_gas_part(s,p,ilevel,action_part)
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
  real(kind=8),dimension(1:ndim)::x,vel
  integer,dimension(1:ndim)::icell_idx
  integer(kind=8),dimension(0:ndim)::hash_nbor
  real(kind=8),dimension(1:2*ndim)::prob_face
  real(kind=8)::out_sum,scale,stay_prob,u,cum
  integer::ipart,idim,iface,selected,icell
  real(kind=8)::rho_cell,denom,dx_loc,dist_to_face
  real(kind=8)::tol_corner,dist_to_corner
  logical :: near_corner
  integer,dimension(1:ndim)::corner_idx
  integer,dimension(1:ndim,1:twotondim)::corner_nbor_idx
  real(kind=8),dimension(1:twotondim)::corner_weight
  real(kind=8)::weight_sum
  integer :: ind,bit,selected_corner,corner_associated
  integer :: ii
  character(LEN=80)::filename,fileloc
  character(LEN=5)::nchar
  type(oct),pointer::gridp
  type(msg_hydro_mflux)::dummy_nvar_realdp
  type(RngStream),external::RngStream_CreateStream
  real(kind=8),external::RngStream_RandUni
  integer(kind=8)::stream_skip
  external :: RngStream_SetPackageSeed, RngStream_AdvanceState

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  if (p%type/=TRAC_TYPE) return

  dx_loc=r%boxlen/2**ilevel
  tol_corner = 0.05d0
  
  if(action_part==action_kick_only)then
    do ipart=p%headp(ilevel),p%tailp(ilevel)
       p%levelp(ipart)=ilevel
    end do
    return
  endif

  if(.not.tracer_rng_ready)then
     call RngStream_SetPackageSeed(r%seed)
     tracer_rng = RngStream_CreateStream('tracer_mc')
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

     near_corner = .true.
     do idim=1,ndim
        dist_to_corner = abs(x(idim) - dble(nint(x(idim))))
        if(dist_to_corner >= tol_corner) near_corner = .false.
     end do

     if(near_corner)then
        do idim=1,ndim
           corner_idx(idim)=nint(x(idim))
        end do
        corner_weight=0.d0
        corner_associated=0
        do ind=1,twotondim
           do idim=1,ndim
              bit = merge(1,0,btest(ind-1,idim-1))
              corner_nbor_idx(idim,ind)=corner_idx(idim)+bit-1
              if(r%periodic(idim))then
                 if(corner_nbor_idx(idim,ind)< m%box_ckey_min(idim,ilevel+1))corner_nbor_idx(idim,ind)=m%box_ckey_max(idim,ilevel+1)-1
                 if(corner_nbor_idx(idim,ind)>=m%box_ckey_max(idim,ilevel+1))corner_nbor_idx(idim,ind)=m%box_ckey_min(idim,ilevel+1)
              end if
           end do
           hash_nbor(0)=ilevel+1
           hash_nbor(1:ndim)=corner_nbor_idx(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
           if(associated(gridp))then
              corner_associated = corner_associated + 1
              corner_weight(ind)=max(gridp%uold(icell,1), r%smallr)
           else
              corner_weight(ind)=0.d0
           end if
#else
           if(associated(gridp))then
              corner_associated = corner_associated + 1
              corner_weight(ind)=1.d0
           else
              corner_weight(ind)=0.d0
           end if
#endif
        end do
        if(corner_associated<twotondim) near_corner = .false.
        weight_sum = 0.d0
        do ind=1,twotondim
           weight_sum = weight_sum + corner_weight(ind)
        end do
        if(near_corner .and. weight_sum>0.d0)then
           u = RngStream_RandUni(tracer_rng)
           cum = 0.d0
           selected_corner = 1
           do ind=1,twotondim
              cum = cum + corner_weight(ind)/weight_sum
              if(u<=cum)then
                 selected_corner = ind
                 exit
              end if
           end do
           do idim=1,ndim
              x(idim)=dble(corner_nbor_idx(idim,selected_corner))+0.5d0
           end do
           call wrap_cell_coords(s,x,ilevel+1)
           do idim=1,ndim
              p%xp(ipart,idim)=x(idim)*dx_loc-m%skip(idim)
           end do
        else
           near_corner = .false.
        end if
     end if

     if(.not.near_corner)then
        do idim=1,ndim
           dist_to_face = abs(x(idim) - dble(nint(x(idim))))
           if(dist_to_face < tol_corner)then
              if(RngStream_RandUni(tracer_rng) < 0.5d0)then
                 x(idim)=dble(nint(x(idim)))-0.5d0
              else
                 x(idim)=dble(nint(x(idim)))+0.5d0
              endif
           else
              x(idim)=dble(int(x(idim)))+0.5d0
           endif
        end do
        call wrap_cell_coords(s,x,ilevel+1)
        do idim=1,ndim
           p%xp(ipart,idim)=x(idim)*dx_loc-m%skip(idim)
        end do
     end if

     do idim=1,ndim
        icell_idx(idim)=int(x(idim))
        if(r%periodic(idim))then
           if(icell_idx(idim)< m%box_ckey_min(idim,ilevel+1))icell_idx(idim)=m%box_ckey_max(idim,ilevel+1)-1
           if(icell_idx(idim)>=m%box_ckey_max(idim,ilevel+1))icell_idx(idim)=m%box_ckey_min(idim,ilevel+1)
        endif
     enddo

     hash_nbor(0)=ilevel+1
     hash_nbor(1:ndim)=icell_idx(1:ndim)
     call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
     if(.not.associated(gridp))cycle

     rho_cell=gridp%mflux(icell,1)
     denom=max(rho_cell,r%smallr)
     vel=0.d0
     vel(1:ndim)=gridp%uold(icell,2:ndim+1)/max(gridp%uold(icell,1),r%smallr)

    prob_face=0.d0
    do idim=1,ndim
       prob_face(2*idim-1)=max(-gridp%mflux(icell,1+idim),0.d0)/denom
       prob_face(2*idim  )=max( gridp%mflux(icell,1+idim+ndim),0.d0)/denom
    end do
#else
     if(.not.associated(gridp))cycle
     prob_face=0.d0
     vel=0.d0
     denom=1.d0
#endif

     out_sum=0.d0
     do iface=1,2*ndim
        if(prob_face(iface)>0.d0)out_sum=out_sum+prob_face(iface)
     end do

     if(out_sum>1.d0)then ! Just for safety, as out_sum SHOULD be the probability of exiting the cell
        scale=1.d0/out_sum
        prob_face=prob_face*scale
        out_sum=1.d0
     end if

     stay_prob=max(0.d0,1.d0-out_sum)

     selected=0
     u=RngStream_RandUni(tracer_rng)
     cum=0.d0
     do iface=1,2*ndim
        cum=cum+prob_face(iface)
        if(u<cum)then
           selected=iface
           exit
        endif
     end do
        ! Remaining probability corresponds to staying in the host cell


     p%vp(ipart,1:ndim)=vel(1:ndim)
     p%levelp(ipart)=ilevel

     if(selected>0)then
        idim = (selected + 1) / 2 ! integer division.
        p%xp(ipart,idim)=p%xp(ipart,idim)+ (-1)**selected * dx_loc
     endif

     ! Trajectory output for selected particles
     if(s%r%ntrajectories>0)then
        do ii=1,s%r%ntrajectories
           if(s%r%trajectories(ii)==p%idp(ipart))then
              call title(g%myid,nchar)
              filename='trajectory.dat'
              fileloc=TRIM(filename)//TRIM(nchar)
              open(25+g%myid,file=fileloc,status='unknown',access='append')
              write(25+g%myid,*) g%t, p%idp(ipart), &
                   (p%xp(ipart,idim), idim=1,ndim), &
                   (p%vp(ipart,idim), idim=1,ndim)
              close(25+g%myid)
              exit
           endif
        end do
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
end subroutine mc_trace_gas_part

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
  integer,dimension(1:ndim)::il,ir
  real(kind=8),dimension(1:twotondim)::vol
  integer,dimension(1:ndim,1:twotondim)::ckey
  real(kind=8),dimension(1:ndim,1:2)::w1d,dw1d
  real(kind=8),dimension(1:ndim)::momentum
  real(kind=8)::rho,kappa_sum
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ind,icell,jd
  type(oct),pointer::gridp

  vel_out=0.d0
  momentum=0.d0
  rho=0.d0
  kappa_sum=0.d0

  ! Build CIC weights, indices, and volume weights
  call cic_weights_and_derivs(x_cell, w1d, dw1d, il, ir, vol)
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
        momentum(1:ndim)=momentum(1:ndim)+gridp%uold(icell,2:ndim+1)*vol(ind)
        rho=rho+gridp%uold(icell,1)*vol(ind)
        if(use_sgs_in)then
           kappa_sum=kappa_sum+tracer_cell_kappa(gridp%uold(icell,1),gridp%uold(icell,st%r%iturb),dx_cell,st%r%smallr)*vol(ind)
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

subroutine gather_tsc_state(st,x_cell,level_in,dx_cell,use_sgs_in,vel_out,kappa_out)
  use amr_parameters, only: ndim, threetondim
  use oct_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache
  use rho_fine_module, only: tsc_index
  implicit none
  type(ramses_t),intent(in)::st
  real(kind=8),intent(in)::x_cell(1:ndim)
  integer,intent(in)::level_in
  real(kind=8),intent(in)::dx_cell
  logical,intent(in)::use_sgs_in
  real(kind=8),intent(out)::vel_out(1:ndim)
  real(kind=8),intent(out)::kappa_out
  integer,dimension(1:ndim)::il,ic,ir
  real(kind=8),dimension(1:threetondim)::vol
  integer,dimension(1:ndim,1:threetondim)::ckey
  real(kind=8),dimension(1:ndim,1:3)::w1d,dw1d
  real(kind=8),dimension(1:ndim)::momentum
  real(kind=8)::rho,kappa_sum
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ind,icell,jd
  type(oct),pointer::gridp

  vel_out=0.d0
  momentum=0.d0
  rho=0.d0
  kappa_sum=0.d0

  ! Build TSC weights, indices, and volume weights
  call tsc_weights_and_derivs(x_cell, w1d, dw1d, il, ic, ir, vol)
  do jd=1,ndim
     if(st%r%periodic(jd))then
        if(il(jd)< st%m%box_ckey_min(jd,level_in+1))il(jd)=il(jd)-st%m%box_ckey_min(jd,level_in+1)+st%m%box_ckey_max(jd,level_in+1)
        if(ir(jd)>=st%m%box_ckey_max(jd,level_in+1))ir(jd)=ir(jd)+st%m%box_ckey_min(jd,level_in+1)-st%m%box_ckey_max(jd,level_in+1)
     endif
  end do
  ckey = tsc_index(il,ic,ir)

  hash_nbor(0)=level_in+1
  do ind=1,threetondim
     hash_nbor(1:ndim)=ckey(1:ndim,ind)
     call get_parent_cell(st,hash_nbor,st%m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
     if(associated(gridp))then
        momentum(1:ndim)=momentum(1:ndim)+gridp%uold(icell,2:ndim+1)*vol(ind)
        rho=rho+gridp%uold(icell,1)*vol(ind)
        if(use_sgs_in)then
           kappa_sum=kappa_sum+tracer_cell_kappa(gridp%uold(icell,1),gridp%uold(icell,st%r%iturb),dx_cell,st%r%smallr)*vol(ind)
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
end subroutine gather_tsc_state

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
  integer,dimension(1:ndim)::il,ir
  integer,dimension(1:ndim,1:twotondim)::ckey
  real(kind=8),dimension(1:ndim,1:2)::w1d,dw1d
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ind,jd,icell
  type(oct),pointer::gridp

  phi_cells=0.d0
  if(.not.use_sgs_in)return

  ! Build CIC weights and indices (only indices used here)
  call cic_weights_and_derivs(x_cell, w1d, dw1d, il, ir)
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
        phi_cells(ind)=tracer_cell_kappa(gridp%uold(icell,1),gridp%uold(icell,st%r%iturb),dx_cell,st%r%smallr)
     end if
#endif
  end do
end subroutine gather_cic_scalar

real(kind=8) function tracer_cell_kappa(dens_in,eturb_in,dx_in,smallr_in) result(kappa_val)
  implicit none
  real(kind=8),intent(in)::dens_in,eturb_in,dx_in,smallr_in
  real(kind=8)::rho_eff,sigma_sq

  rho_eff = max(dens_in,smallr_in)
  sigma_sq = max(2.0d0*max(eturb_in,0.0d0)/rho_eff,0.0d0)
  if(sigma_sq>0.0d0)then
     kappa_val = dx_in*sqrt(sigma_sq)
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

subroutine sample_tracer_two_piece_uniform(vec, gamma1_vec)
  !------------------------------------------------------------------------
  ! Sample from a two-piece uniform distribution with mean 0, variance 1,
  ! and skewness gamma1 (one value per dimension).
  !
  ! The PDF is piecewise constant on [-a, 0] and [0, b], with heights chosen
  ! so that mean=0, variance=1, and third central moment = gamma1.
  !
  ! Given pr = max(mflux(1+idim+ndim), 0) / mflux(1)  (right exit prob),
  !       pl = max(-mflux(1+idim), 0) / mflux(1)      (left exit prob),
  ! the skewness of the MC discrete kernel is:
  !   gamma1 = (pr - pl)*(1 - 3*pr + 4*pr**2 - 2*pr**3 - 3*pl - 8*pr*pl
  !            + 2*pr**2*pl + 4*pl**2 + 2*pr*pl**2 - 2*pl**3)
  !------------------------------------------------------------------------
  use amr_parameters, only: ndim
  implicit none
  real(kind=8),intent(out)::vec(1:ndim)
  real(kind=8),intent(in)::gamma1_vec(1:ndim)
  integer :: jd
  real(kind=8)::u_rand,gamma1,k,a,b,p_left,apb
  real(kind=8), external :: RngStream_RandUni

  vec=0.0d0
  do jd=1,ndim
     gamma1 = gamma1_vec(jd)

     ! Compute two-piece uniform parameters from skewness
     ! Constraints: ab = 3, b - a = 4*gamma1/3
     k = 4.0d0 * gamma1 / 3.0d0
     a = 0.5d0 * (-k + sqrt(k*k + 12.0d0))
     b = a + k
     apb = a + b
     p_left = b / apb   ! probability mass on the left piece [-a, 0]

     ! Sample using inverse CDF
     u_rand = RngStream_RandUni(tracer_rng)
     if (u_rand < p_left) then
        ! Left piece: uniform on [-a, 0]
        vec(jd) = -a + a * u_rand * apb / b
     else
        ! Right piece: uniform on [0, b]
        vec(jd) = b * (u_rand * apb - b) / a
     endif
  end do
end subroutine sample_tracer_two_piece_uniform

real(kind=8) function mc_kernel_skewness(pr, pl) result(gamma1)
  !------------------------------------------------------------------------
  ! Compute the standardized skewness of the MC discrete kernel for one
  ! dimension: gamma1 = mu3 / variance^(3/2).
  !
  ! The MC kernel has X in {-dx, 0, +dx} with probabilities (pl, 1-pl-pr, pr).
  ! Third central moment (mu3) and variance have closed forms in (pl, pr);
  ! standardized skewness = mu3 / variance^1.5.
  !
  ! Input: pr = max(mflux(1+idim+ndim), 0) / mflux(1)  (right exit prob)
  !        pl = max(-mflux(1+idim), 0) / mflux(1)      (left exit prob)
  !------------------------------------------------------------------------
  implicit none
  real(kind=8),intent(in)::pr, pl
  real(kind=8)::mu3, variance

  mu3 = (pr - pl) * (1.0d0 - 3.0d0*pl + 4.0d0*pl**2 - 2.0d0*pl**3 &
       & - 3.0d0*pr - 8.0d0*pl*pr + 2.0d0*pl**2*pr &
       & + 4.0d0*pr**2 + 2.0d0*pl*pr**2 - 2.0d0*pr**3)
  variance = pl - pl**2 + pr + 2.0d0*pl*pr - pr**2
  if (variance > 0.0d0) then
     gamma1 = mu3 / (variance**1.5d0)
  else
     gamma1 = 0.0d0
  endif
end function mc_kernel_skewness

real(kind=8) function mc_kernel_kurtosis(pr, pl) result(gamma2)
  !------------------------------------------------------------------------
  ! Compute the excess kurtosis of the MC discrete kernel for one dimension.
  ! Excess kurtosis = kurtosis - 3 (zero for Gaussian).
  !
  ! The MC kernel has X in {-dx, 0, +dx} with probabilities (pl, 1-pl-pr, pr).
  ! Closed form: gamma2 = -6 + (pl - pl^2 + pr + 14*pl*pr - pr^2) / variance^2.
  !
  ! Input: pr = max(mflux(1+idim+ndim), 0) / mflux(1)  (right exit prob)
  !        pl = max(-mflux(1+idim), 0) / mflux(1)      (left exit prob)
  !------------------------------------------------------------------------
  implicit none
  real(kind=8),intent(in)::pr, pl
  real(kind=8)::variance

  variance = pl - pl**2 + pr + 2.0d0*pl*pr - pr**2
  if (variance > 0.0d0) then
     gamma2 = -6.0d0 + (pl - pl**2 + pr + 14.0d0*pl*pr - pr**2) / (variance**2)
  else
     gamma2 = 0.0d0
  endif
end function mc_kernel_kurtosis

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
  real(kind=8),dimension(1:ndim)::x,x_mid
  integer,dimension(1:ndim)::il,ic,ir
  real(kind=8),dimension(1:threetondim)::vol,vol2
  integer,dimension(1:ndim,1:threetondim)::ckey,ckey2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  real(kind=8),dimension(1:ndim,1:3)::w1d,dw1d
  integer::ipart,icell,icell2,ind,idim
  real(kind=8)::dx_loc
  real(kind=8),dimension(1:ndim)::ff,v_pred
  type(oct),pointer::gridp
  type(msg_three_realdp)::dummy_three_realdp
  type(msg_hydro_mflux)::dummy_nvar_realdp
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
     ! Build 1D TSC weights, indices, and 3D volume weights
     call tsc_weights_and_derivs(x, w1d, dw1d, il, ic, ir, vol)
     do idim=1,ndim
        if(r%periodic(idim))then
           if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=m%box_ckey_max(idim,ilevel+1)-1
           if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=m%box_ckey_min(idim,ilevel+1)
        endif
     enddo
     ckey = tsc_index(il,ic,ir)
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
        ! Build 1D TSC weights, indices, and 3D volume weights at x_mid
        call tsc_weights_and_derivs(x_mid, w1d, dw1d, il, ic, ir, vol2)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=m%box_ckey_max(idim,ilevel+1)-1
              if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=m%box_ckey_min(idim,ilevel+1)
           endif
        enddo
        ckey2 = tsc_index(il,ic,ir)
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
  real(kind=8),dimension(1:ndim)::x,x_mid,wll,wl,wr,wrr
  integer,dimension(1:ndim)::cll,cl,cr,crr
  real(kind=8),dimension(1:fourtondim)::vol,vol2
  integer,dimension(1:ndim,1:fourtondim)::ckey,ckey2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,icell,icell2,ind,idim
  real(kind=8)::dx_loc
  real(kind=8),dimension(1:ndim)::ff,v_pred
  type(oct),pointer::gridp
  type(msg_hydro_mflux)::dummy_nvar_realdp
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
     ! Build PCS 1D weights and indices
     call pcs_1d_weights(x, wll, wl, wr, wrr, cll, cl, cr, crr)
     vol = pcs_weight(wll, wl, wr, wrr)
     do idim=1,ndim
        if(r%periodic(idim))then
           if(cll(idim)< m%box_ckey_min(idim,ilevel+1))cll(idim)=cll(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
           if(cl (idim)< m%box_ckey_min(idim,ilevel+1))cl (idim)=cl (idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
           if(cr (idim)>=m%box_ckey_max(idim,ilevel+1))cr (idim)=cr (idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
           if(crr(idim)>=m%box_ckey_max(idim,ilevel+1))crr(idim)=crr(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
        endif
     enddo
     ckey = pcs_index(cll,cl,cr,crr)
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
        ! Build PCS 1D weights and indices at x_mid
        call pcs_1d_weights(x_mid, wll, wl, wr, wrr, cll, cl, cr, crr)
        vol2 = pcs_weight(wll, wl, wr, wrr)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(cll(idim)< m%box_ckey_min(idim,ilevel+1))cll(idim)=cll(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
              if(cl (idim)< m%box_ckey_min(idim,ilevel+1))cl (idim)=cl (idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
              if(cr (idim)>=m%box_ckey_max(idim,ilevel+1))cr (idim)=cr (idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
              if(crr(idim)>=m%box_ckey_max(idim,ilevel+1))crr(idim)=crr(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
           endif
        enddo
        ckey2 = pcs_index(cll,cl,cr,crr)
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
subroutine trace_gas_part_trivial(s,p,ilevel,action_part)
  use amr_parameters, only: ndim
  use pm_parameters
  use pm_commons, only: part_t
  use ramses_commons, only: ramses_t
  use rng
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part
  real(kind=8),dimension(1:ndim)::disp,xi
  real(kind=8)::dx_loc,dt_level,D_diff,noise_amp
  integer :: ipart,idim
  integer :: ii
  integer :: z_dim
  character(LEN=80)::filename,fileloc
  character(LEN=5)::nchar
  type(RngStream),external::RngStream_CreateStream
  real(kind=8),external::RngStream_RandUni
  integer(kind=8)::stream_skip
  external :: RngStream_SetPackageSeed, RngStream_AdvanceState

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  if (p%type/=TRAC_TYPE) return

  dx_loc=r%boxlen/2**ilevel
  dt_level=g%dtnew(ilevel)
  
  ! Determine z-direction (last dimension)
  z_dim = ndim

  ! Initialize RNG if needed
  if(.not.tracer_rng_ready)then
     call RngStream_SetPackageSeed(r%seed)
     tracer_rng = RngStream_CreateStream('tracer_trivial')
     stream_skip = int(2*g%myid,kind=8)
     call RngStream_AdvanceState(tracer_rng,0_8,stream_skip)
     tracer_rng_ready = .true.
  end if

  do ipart=p%headp(ilevel),p%tailp(ilevel)

     ! Set velocity: 1.0 in z-direction only, 0.0 in others
     do idim=1,ndim
        if(idim==z_dim)then
           p%vp(ipart,idim)=1.0d0
        else
           p%vp(ipart,idim)=0.0d0
        endif
     end do
     p%levelp(ipart)=ilevel

     if(action_part==action_kick_only)then
        cycle
     endif

     ! Compute displacement from advection
     disp(1:ndim)=p%vp(ipart,1:ndim)*dt_level

     ! Add random diffusion kick in z-direction only
     ! Diffusion coefficient: D = dx_loc*(1.d0-dt_level/dx_loc)/2.d0
     D_diff = 0.00653!dx_loc*(1.d0-dt_level/dx_loc)/2.d0
     if(D_diff > 0.0d0)then
        call sample_tracer_uniform(xi)
        noise_amp = sqrt(2.0d0*D_diff*dt_level)
        disp(z_dim)=disp(z_dim)+noise_amp*xi(z_dim)
     end if

     ! Update position
     p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+disp(1:ndim)

     ! Trajectory output for selected particles
     if(s%r%ntrajectories>0)then
        do ii=1,s%r%ntrajectories
           if(s%r%trajectories(ii)==p%idp(ipart))then
              call title(g%myid,nchar)
              filename='trajectory.dat'
              fileloc=TRIM(filename)//TRIM(nchar)
              open(25+g%myid,file=fileloc,status='unknown',access='append')
              write(25+g%myid,*) g%t, p%idp(ipart), &
                   (p%xp(ipart,idim), idim=1,ndim), &
                   (p%vp(ipart,idim), idim=1,ndim)
              close(25+g%myid)
              exit
           endif
        end do
     endif

  end do

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

end subroutine trace_gas_part_trivial
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
  real(kind=8),dimension(1:ndim)::x,x_mid
  integer,dimension(1:ndim)::ir,il
  integer,dimension(1:ndim)::ir2,il2
  real(kind=8),dimension(1:twotondim)::vol,vol2
  integer,dimension(1:ndim,1:twotondim)::ckey,ckey2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  real(kind=8),dimension(1:ndim,1:2)::w1d,dw1d
  integer::ipart,ind,idim,irad,icell,icell2
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

  ! OPTIMIZATION: Handle kick-only pass immediately without cache operations
  if(action_part==action_kick_only)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        p%levelp(ipart)=ilevel
     end do
     return
  endif

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
     call cic_weights_and_derivs(x, w1d, dw1d, il, ir, vol)
     do idim=1,ndim
        if(il(idim)<0)il(idim)=m%ckey_max(ilevel+1)-1
        if(ir(idim)==m%ckey_max(ilevel+1))ir(idim)=0
     enddo
     ckey = cic_index(il,ir)

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
     call cic_weights_and_derivs(x_mid, w1d, dw1d, il2, ir2, vol2)
     do idim=1,ndim
        if(il2(idim)<0)il2(idim)=m%ckey_max(ilevel+1)-1
        if(ir2(idim)==m%ckey_max(ilevel+1))ir2(idim)=0
     enddo
     ckey2 = cic_index(il2,ir2)

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
  end do
  call close_cache(s,m%grid_dict)
  
  end associate
end subroutine cic_kick_drift_dust


! This is missing some key things. 
! First, we need dx_ito (which should be renamed to dx_ito)
! to include a grad(kappa).
! Second, we need to have the noise logic be like sqrt(2*kappa*dt).
!subroutine cic_kick_drift_dust_num_diff(s,p,ilevel,action_part)
!end subroutine cic_kick_drift_dust_num_diff

subroutine cic_kick_drift_dust_ito_mc(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, twotondim
  use hydro_parameters, only: nener
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
  real(kind=8),dimension(1:ndim)::x_main,x_mid
  integer,dimension(1:ndim)::cl2,cr2
  real(kind=8),dimension(1:twotondim)::vol2
  integer,dimension(1:ndim,1:twotondim)::ckey2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,ind,idim,irad,icell,icell2
  real(kind=8)::dx_loc,dt_level
  real(kind=8),dimension(1:ndim)::ff,v_pred,uu
#ifdef MHD
  real(kind=8),dimension(1:3)::bb
  real(kind=8)::emag
#endif
  real(kind=8)::rho_gas,c_sound,eint,coeff
  real(kind=8)::nu_stop,dens,etot,ekin,erad,cs2,pi
  real(kind=8),dimension(1:ndim)::what,wdrift,wdrift0
  real(kind=8),dimension(1:ndim)::disp,xi,dx_ito
  real(kind=8),dimension(1:ndim)::x_diff,grad_at_part,u_eff,kappa_num
  integer,dimension(1:ndim)::il,ir
  real(kind=8),dimension(1:twotondim)::vol_diff
  integer,dimension(1:ndim,1:twotondim)::ckey_diff
  real(kind=8),dimension(1:ndim,1:twotondim)::u_cells,kappa_num_cells
  real(kind=8),dimension(1:ndim,1:2)::w1d,dw1d,w1d_mid,dw1d_mid
  real(kind=8)::wdrift2,w0,a,g_par,g_perp
  type(oct),pointer :: gridp
  integer :: k,jdim
  real(kind=8)::rho,denom,fluxL,fluxR
  real(kind=8)::jr,jl
  type(msg_large_realdp)::dummy_large_realdp
  type(RngStream),external::RngStream_CreateStream
  real(kind=8),external::RngStream_RandUni
  integer(kind=8)::stream_skip
  external :: RngStream_SetPackageSeed, RngStream_AdvanceState, gaussdev

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  dx_loc=r%boxlen/2**ilevel
  dt_level=g%dtnew(ilevel)
  if (p%type/=DUST_TYPE) return
  pi=4.0d0*atan(1.0d0)
  coeff=9.0d0*pi*r%gamma/128.0d0

  ! OPTIMIZATION: Handle kick-only pass immediately without cache operations
  if(action_part==action_kick_only)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        p%levelp(ipart)=ilevel
     end do
     return
  endif

  if(.not.tracer_rng_ready)then
     call RngStream_SetPackageSeed(r%seed)
     tracer_rng = RngStream_CreateStream('tracer_sgs')
     stream_skip = int(2*g%myid,kind=8)
     call RngStream_AdvanceState(tracer_rng,0_8,stream_skip)
     tracer_rng_ready = .true.
  end if

  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_large_realdp)/32,&
                     pack=pack_fetch_kick_dust,unpack=unpack_fetch_kick_dust)

  do ipart=p%headp(ilevel),p%tailp(ilevel)
     ! Compute particle position in cell units
     do idim=1,ndim
        x_main(idim)=p%xp(ipart,idim)/dx_loc
     end do
     do idim=1,ndim
        if(x_main(idim)<0d0)x_main(idim)=x_main(idim)+dble(m%ckey_max(ilevel+1))
        if(x_main(idim)>=dble(m%ckey_max(ilevel+1)))x_main(idim)=x_main(idim)-dble(m%ckey_max(ilevel+1))
     end do

     if(action_part==action_kick_drift)then
        v_pred(1:ndim)=p%vp(ipart,1:ndim)
        if (g%nstep==0) then
           do idim=1,ndim
              x_mid(idim)=x_main(idim)+0.5d0*dt_level*v_pred(idim)/dx_loc
           end do
        else
           do idim=1,ndim
              x_mid(idim)=x_main(idim)
           end do
        endif
        do idim=1,ndim
           if(x_mid(idim)<0d0)x_mid(idim)=x_mid(idim)+dble(m%ckey_max(ilevel+1))
           if(x_mid(idim)>=dble(m%ckey_max(ilevel+1)))x_mid(idim)=x_mid(idim)-dble(m%ckey_max(ilevel+1))
        end do

        ! Build weights and indices for x_mid (derivatives not needed but computed anyway)
        call cic_weights_and_derivs(x_mid, w1d_mid, dw1d_mid, cl2, cr2, vol2)
        do idim=1,ndim
           if(cl2(idim)<0)cl2(idim)=m%ckey_max(ilevel+1)-1
           if(cr2(idim)==m%ckey_max(ilevel+1))cr2(idim)=0
        end do
        ckey2 = cic_index(cl2,cr2)

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
#ifdef MHD
              bb(1:3)=bb(1:3)+0.5d0*(gridp%bold(icell2,1:3)+gridp%bold(icell2,4:6))*vol2(ind)
#endif
              dens = max(dble(gridp%mflux(icell2,1)), r%smallr) ! Must use mflux(1) for computing flux-based velocities

              do idim=1,ndim
                 fluxL = gridp%mflux(icell2,1+idim     )*dx_loc/dt_level
                 fluxR = gridp%mflux(icell2,1+idim+ndim)*dx_loc/dt_level
                 jr = max(fluxR,0.d0)
                 jl = max(-fluxL,0.d0)
                 uu(idim)=uu(idim)+((jr-jl)/dens)*vol2(ind)
              end do
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
        wdrift0(1:ndim)= v_pred(1:ndim) - uu(1:ndim)
        wdrift(1:ndim)=wdrift0(1:ndim)
        call compute_drag_step(wdrift, c_sound, 0.5d0*dt_level, nu_stop, coeff, r%analytic_dust_force)
#ifdef GRAV
        wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*dt_level
#endif
#ifdef MHD
        if (ndim==3) then
            call compute_lorentz_step(wdrift, bb(1:ndim), dt_level, p%charge(ipart), r%analytic_dust_force)
        else
            if(g%myid==1 .and. g%nstep==0)then
               write(*,*) 'Warning: Lorentz force not implemented for NDIM != 3; proceeding without it.'
            endif
        endif
#endif
#ifdef GRAV
        wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*dt_level
#endif
        call compute_drag_step(wdrift, c_sound, 0.5d0*dt_level, nu_stop, coeff, r%analytic_dust_force)

        p%vp(ipart,1:ndim)=uu(1:ndim)+wdrift(1:ndim)

        do idim=1,ndim
           x_diff(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
        end do
        call wrap_cell_coords(s,x_diff,ilevel+1)

        ! Build 1D CIC weights, derivatives, indices, and 3D volume weights
        call cic_weights_and_derivs(x_diff, w1d, dw1d, il, ir, vol_diff)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=m%box_ckey_max(idim,ilevel+1)-1
              if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=m%box_ckey_min(idim,ilevel+1)
           endif
        enddo

        ! Build neighbor index list (twotondim=8)
        ckey_diff = cic_index(il,ir)

        ! Compute cell-by-cell velocities and diffusivities
        u_cells=0.d0
        kappa_num_cells=0.d0
        hash_nbor(0)=ilevel+1
        do ind=1,twotondim
           hash_nbor(1:ndim)=ckey_diff(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
           if(associated(gridp))then
              rho=gridp%mflux(icell,1)
              denom=max(rho,r%smallr)
              do idim=1,ndim
                 fluxL=gridp%mflux(icell,1+idim     )*dx_loc/dt_level
                 fluxR=gridp%mflux(icell,1+idim+ndim)*dx_loc/dt_level
                 jr=max(fluxR,0.d0)
                 jl=max(-fluxL,0.d0)
                 u_cells(idim,ind)=(jr-jl)/denom
                 kappa_num_cells(idim,ind)=0.5d0*(jr+jl)*dx_loc/denom
              end do
           end if
#endif
        end do

        u_eff=0.d0
        kappa_num=0.d0
        do ind=1,twotondim
           do idim=1,ndim
              u_eff(idim)=u_eff(idim)+u_cells(idim,ind)*vol_diff(ind)
              kappa_num(idim)=kappa_num(idim)+kappa_num_cells(idim,ind)*vol_diff(ind)
           end do
        end do

        ! Compute gradient using CIC subroutine
        grad_at_part(1:ndim)=0.d0
        !call compute_gradient_cic(w1d, dw1d, kappa_num_cells, grad_at_part)
        !u_eff(1:ndim) = u_eff(1:ndim) + grad_at_part(1:ndim)

        call sample_tracer_uniform(xi)

        disp(1:ndim)=0.d0
        do idim=1,ndim
           disp(idim)=u_eff(idim)*dt_level
        end do

        wdrift2 = dot_product(wdrift0(1:ndim), wdrift0(1:ndim))
        if (wdrift2 > 0.d0) then
           what(1:ndim) = wdrift0(1:ndim)/sqrt(wdrift2)
           w0 = sqrt(wdrift2)
        else
           what(1:ndim) = 0.d0
           w0 = 0.d0
        end if

        if (nu_stop>0.d0 .and. w0 > 0.d0) then
           a = c_sound/sqrt(coeff)
           if(dt_level*nu_stop < 1.d2)then
              call get_g_factors(dt_level, w0, a, nu_stop, g_par, g_perp)
           else
              g_par = sqrt(1.d0-3.d0/(2.d0*dt_level*nu_stop) &
              &- exp(-2.d0*dt_level*nu_stop)/(2.d0*dt_level*nu_stop) &
              &+ 2.d0*exp(-(dt_level*nu_stop))/(dt_level*nu_stop))
              g_perp = g_par
           end if
        else
           g_par = 0.d0
           g_perp = 0.d0
        end if

        dx_ito(1:ndim)=0.d0
        do idim=1,ndim
           do jdim=1,ndim
              if (wdrift2 > 0.d0) then
                 dx_ito(idim) = dx_ito(idim) + &
                      (g_perp * merge(1.d0,0.d0,idim==jdim) + (g_par - g_perp)*what(idim)*what(jdim)) * &
                      (sqrt(2.d0*kappa_num(jdim)*dt_level) * xi(jdim) + grad_at_part(jdim) * dt_level)
              end if
           end do
        end do

        disp(1:ndim)=disp(1:ndim)+dx_ito(1:ndim)

        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+disp(1:ndim)
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
end subroutine cic_kick_drift_dust_ito_mc

subroutine tsc_kick_drift_dust_ito_mc(s,p,ilevel,action_part)
  use amr_parameters, only: ndim, threetondim
  use hydro_parameters, only: nener
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
  real(kind=8),dimension(1:ndim)::x_main,x_mid
  integer,dimension(1:ndim)::cl2,cc2,cr2
  real(kind=8),dimension(1:threetondim)::vol2
  integer,dimension(1:ndim,1:threetondim)::ckey2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,ind,idim,irad,icell,icell2
  real(kind=8)::dx_loc,dt_level
  real(kind=8),dimension(1:ndim)::ff,v_pred,uu
#ifdef MHD
  real(kind=8),dimension(1:3)::bb
  real(kind=8)::emag
#endif
  real(kind=8)::rho_gas,c_sound,eint,coeff
  real(kind=8)::nu_stop,dens,etot,ekin,erad,cs2,pi
  real(kind=8),dimension(1:ndim)::what,wdrift,wdrift0
  real(kind=8),dimension(1:ndim)::disp,xi,dx_ito
  real(kind=8),dimension(1:ndim)::x_diff,grad_at_part,u_eff,kappa_num
  integer,dimension(1:ndim)::il,ic,ir
  real(kind=8),dimension(1:threetondim)::vol_diff
  integer,dimension(1:ndim,1:threetondim)::ckey_diff
  real(kind=8),dimension(1:ndim,1:threetondim)::u_cells,kappa_num_cells
  real(kind=8),dimension(1:ndim,1:3)::w1d,dw1d,w1d_mid,dw1d_mid
  real(kind=8)::wdrift2,w0,a,g_par,g_perp
  type(oct),pointer :: gridp
  integer :: k,jdim
  real(kind=8)::rho,denom,fluxL,fluxR
  real(kind=8)::jr,jl
  type(msg_large_realdp)::dummy_large_realdp
  type(RngStream),external::RngStream_CreateStream
  real(kind=8),external::RngStream_RandUni
  integer(kind=8)::stream_skip
  external :: RngStream_SetPackageSeed, RngStream_AdvanceState, gaussdev

  associate(r=>s%r,g=>s%g,m=>s%m)
  if(p%static)return
  dx_loc=r%boxlen/2**ilevel
  dt_level=g%dtnew(ilevel)
  if (p%type/=DUST_TYPE) return
  pi=4.0d0*atan(1.0d0)
  coeff=9.0d0*pi*r%gamma/128.0d0

  ! OPTIMIZATION: Handle kick-only pass immediately without cache operations
  if(action_part==action_kick_only)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        p%levelp(ipart)=ilevel
     end do
     return
  endif

  if(.not.tracer_rng_ready)then
     call RngStream_SetPackageSeed(r%seed)
     tracer_rng = RngStream_CreateStream('tracer_sgs')
     stream_skip = int(2*g%myid,kind=8)
     call RngStream_AdvanceState(tracer_rng,0_8,stream_skip)
     tracer_rng_ready = .true.
  end if

  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_large_realdp)/32,&
                     pack=pack_fetch_kick_dust,unpack=unpack_fetch_kick_dust)

  do ipart=p%headp(ilevel),p%tailp(ilevel)
     ! Compute particle position in cell units
     do idim=1,ndim
        x_main(idim)=p%xp(ipart,idim)/dx_loc
     end do
     do idim=1,ndim
        if(x_main(idim)<0d0)x_main(idim)=x_main(idim)+dble(m%ckey_max(ilevel+1))
        if(x_main(idim)>=dble(m%ckey_max(ilevel+1)))x_main(idim)=x_main(idim)-dble(m%ckey_max(ilevel+1))
     end do

     ! NOTE: TSC weights at x_main (ckey, vol_main) were removed as they were unused.
     ! The code uses x_mid for gas property interpolation and x_diff for diffusion.

     if(action_part==action_kick_drift)then
        v_pred(1:ndim)=p%vp(ipart,1:ndim)
        if (g%nstep==0) then
           do idim=1,ndim
              x_mid(idim)=x_main(idim)+0.5d0*dt_level*v_pred(idim)/dx_loc
           end do
        else
           do idim=1,ndim
              x_mid(idim)=x_main(idim)
           end do
        endif
        do idim=1,ndim
           if(x_mid(idim)<0d0)x_mid(idim)=x_mid(idim)+dble(m%ckey_max(ilevel+1))
           if(x_mid(idim)>=dble(m%ckey_max(ilevel+1)))x_mid(idim)=x_mid(idim)-dble(m%ckey_max(ilevel+1))
        end do

        ! Build weights and indices for x_mid (derivatives not needed but computed anyway)
        call tsc_weights_and_derivs(x_mid, w1d_mid, dw1d_mid, cl2, cc2, cr2, vol2)
        do idim=1,ndim
           if(cl2(idim)<0)cl2(idim)=m%ckey_max(ilevel+1)-1
           if(cr2(idim)==m%ckey_max(ilevel+1))cr2(idim)=0
        end do
        ckey2 = tsc_index(cl2,cc2,cr2)

        ff(1:ndim)=0.0
        uu(1:ndim)=0.0
        rho_gas = 0.0
        eint = 0.0
#ifdef MHD
        bb(1:3)=0.0
        emag=0.0d0
#endif
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
#ifdef MHD
              bb(1:3)=bb(1:3)+0.5d0*(gridp%bold(icell2,1:3)+gridp%bold(icell2,4:6))*vol2(ind)
#endif
              dens = max(dble(gridp%mflux(icell2,1)), r%smallr) ! Must use mflux(1) for computing flux-based velocities

              do idim=1,ndim
                 fluxL = gridp%mflux(icell2,1+idim     )*dx_loc/dt_level
                 fluxR = gridp%mflux(icell2,1+idim+ndim)*dx_loc/dt_level
                 jr = max(fluxR,0.d0)
                 jl = max(-fluxL,0.d0)
                 uu(idim)=uu(idim)+((jr-jl)/dens)*vol2(ind)
              end do
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
        wdrift0(1:ndim)= v_pred(1:ndim) - uu(1:ndim)
      !   if(dt_level*nu_stop > 1.d1)then
      !      wdrift0(1:ndim) = 0.d0
      !   end if
        wdrift(1:ndim)=wdrift0(1:ndim)
        call compute_drag_step(wdrift, c_sound, 0.5d0*dt_level, nu_stop, coeff, r%analytic_dust_force)
#ifdef GRAV
        wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*dt_level
#endif
#ifdef MHD
        if (ndim==3) then
            call compute_lorentz_step(wdrift, bb(1:ndim), dt_level, p%charge(ipart), r%analytic_dust_force)
        else
            if(g%myid==1 .and. g%nstep==0)then
               write(*,*) 'Warning: Lorentz force not implemented for NDIM != 3; proceeding without it.'
            endif
        endif
#endif
#ifdef GRAV
        wdrift(1:ndim)=wdrift(1:ndim)+ff(1:ndim)*0.5d0*dt_level
#endif
        call compute_drag_step(wdrift, c_sound, 0.5d0*dt_level, nu_stop, coeff, r%analytic_dust_force)

        p%vp(ipart,1:ndim)=uu(1:ndim)+wdrift(1:ndim)

        do idim=1,ndim
           x_diff(idim)=(p%xp(ipart,idim)+m%skip(idim))/dx_loc
        end do
        call wrap_cell_coords(s,x_diff,ilevel+1)

        ! Build 1D TSC weights, derivatives, indices, and 3D volume weights
        call tsc_weights_and_derivs(x_diff, w1d, dw1d, il, ic, ir, vol_diff)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(il(idim)< m%box_ckey_min(idim,ilevel+1))il(idim)=m%box_ckey_max(idim,ilevel+1)-1
              if(ir(idim)>=m%box_ckey_max(idim,ilevel+1))ir(idim)=m%box_ckey_min(idim,ilevel+1)
           endif
        enddo

        ! Build neighbor index list (threetondim=27)
        ckey_diff = tsc_index(il,ic,ir)

        ! Compute cell-by-cell velocities and diffusivities
        u_cells=0.d0
        kappa_num_cells=0.d0
        hash_nbor(0)=ilevel+1
        do ind=1,threetondim
           hash_nbor(1:ndim)=ckey_diff(1:ndim,ind)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.false.,fetch_cache=.true.)
#ifdef HYDRO
           if(associated(gridp))then
              rho=gridp%mflux(icell,1)
              denom=max(rho,r%smallr)
              do idim=1,ndim
                 fluxL=gridp%mflux(icell,1+idim     )*dx_loc/dt_level
                 fluxR=gridp%mflux(icell,1+idim+ndim)*dx_loc/dt_level
                 jr=max(fluxR,0.d0)
                 jl=max(-fluxL,0.d0)
                 u_cells(idim,ind)=(jr-jl)/denom
                 kappa_num_cells(idim,ind)=0.5d0*(jr+jl)*dx_loc/denom
              end do
           end if
#endif
        end do

        u_eff=0.d0
        kappa_num=0.d0
        do ind=1,threetondim
           do idim=1,ndim
              u_eff(idim)=u_eff(idim)+u_cells(idim,ind)*vol_diff(ind)
              kappa_num(idim)=kappa_num(idim)+kappa_num_cells(idim,ind)*vol_diff(ind)
           end do
        end do

        ! Compute gradient using TSC subroutine
        grad_at_part(1:ndim)=0.d0
        !call compute_gradient_tsc(w1d, dw1d, kappa_num_cells, grad_at_part)
        !u_eff(1:ndim) = u_eff(1:ndim) + grad_at_part(1:ndim)

        call sample_tracer_uniform(xi)

        disp(1:ndim)=0.d0
        do idim=1,ndim
           disp(idim)=u_eff(idim)*dt_level
        end do

        wdrift2 = dot_product(wdrift0(1:ndim), wdrift0(1:ndim))
        if (wdrift2 > 0.d0) then
           what(1:ndim) = wdrift0(1:ndim)/sqrt(wdrift2)
           w0 = sqrt(wdrift2)
        else
           what(1:ndim) = 0.d0
           w0 = 0.d0
        end if

        if (nu_stop>0.d0 .and. w0 > 0.d0) then
           a = c_sound/sqrt(coeff)
           if(dt_level*nu_stop < 1.d2)then
              call get_g_factors(dt_level, w0, a, nu_stop, g_par, g_perp)
           else
              g_par = sqrt(1.d0-3.d0/(2.d0*dt_level*nu_stop) &
              &- exp(-2.d0*dt_level*nu_stop)/(2.d0*dt_level*nu_stop) &
              &+ 2.d0*exp(-(dt_level*nu_stop))/(dt_level*nu_stop))
              g_perp = g_par
           end if
        else
           g_par = 0.d0
           g_perp = 0.d0
        end if

        dx_ito(1:ndim)=0.d0
        do idim=1,ndim
           do jdim=1,ndim
              if (wdrift2 > 0.d0) then
                 dx_ito(idim) = dx_ito(idim) + &
                      (g_perp * merge(1.d0,0.d0,idim==jdim) + (g_par - g_perp)*what(idim)*what(jdim)) * &
                      (sqrt(2.d0*kappa_num(jdim)*dt_level) * xi(jdim) + grad_at_part(jdim) * dt_level)
              end if
           end do
        end do

        disp(1:ndim)=disp(1:ndim)+dx_ito(1:ndim)

        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+disp(1:ndim)
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
end subroutine tsc_kick_drift_dust_ito_mc

subroutine cic_kick_drift_dust_guiding_center(s,p,ilevel,action_part)
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use mdl_module
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part

  ! Guiding-center scheme placeholder.
  ! This must not silently do nothing if selected.
  write(*,*) 'Error: dust_force_interpolation_scheme=6 (CIC guiding center) is not implemented.'
  call mdl_abort(s%mdl)
end subroutine cic_kick_drift_dust_guiding_center

subroutine tsc_kick_drift_dust_guiding_center(s,p,ilevel,action_part)
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use mdl_module
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  integer::action_part

  ! Guiding-center scheme placeholder.
  ! This must not silently do nothing if selected.
  write(*,*) 'Error: dust_force_interpolation_scheme=7 (TSC guiding center) is not implemented.'
  call mdl_abort(s%mdl)
end subroutine tsc_kick_drift_dust_guiding_center
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
  real(kind=8),dimension(1:ndim)::x,x_mid
  integer,dimension(1:ndim)::il,ic,ir
  real(kind=8),dimension(1:threetondim)::vol2
  integer,dimension(1:ndim,1:threetondim)::ckey2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  real(kind=8),dimension(1:ndim,1:3)::w1d,dw1d
  integer::ipart,icell2,ind,idim,irad
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

  ! OPTIMIZATION: Handle kick-only pass immediately without cache operations
  if(action_part==action_kick_only)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        p%levelp(ipart)=ilevel
     end do
     return
  endif

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

     if(action_part==action_kick_drift)then
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
        ! Build 1D TSC weights, indices, and 3D volume weights at x_mid
        call tsc_weights_and_derivs(x_mid, w1d, dw1d, il, ic, ir, vol2)
        do idim=1,ndim
           if(il(idim)<0)il(idim)=m%ckey_max(ilevel+1)-1
           if(ir(idim)==m%ckey_max(ilevel+1))ir(idim)=0
        enddo
        ckey2 = tsc_index(il,ic,ir)

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
  real(kind=8),dimension(1:ndim)::x,x_mid,wll,wl,wr,wrr
  integer,dimension(1:ndim)::cll,cl,cr,crr
  real(kind=8),dimension(1:fourtondim)::vol,vol2
  integer,dimension(1:ndim,1:fourtondim)::ckey,ckey2
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ipart,icell,icell2,ind,idim
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

  ! OPTIMIZATION: Handle kick-only pass immediately without cache operations
  if(action_part==action_kick_only)then
     do ipart=p%headp(ilevel),p%tailp(ilevel)
        p%levelp(ipart)=ilevel
     end do
     return
  endif

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
     ! Build PCS 1D weights and indices
     call pcs_1d_weights(x, wll, wl, wr, wrr, cll, cl, cr, crr)
     vol = pcs_weight(wll, wl, wr, wrr)
     do idim=1,ndim
        if(cll(idim)<0)cll(idim)=cll(idim)+m%ckey_max(ilevel+1)
        if(cl (idim)<0)cl (idim)=cl (idim)+m%ckey_max(ilevel+1)
        if(cr (idim)>=m%ckey_max(ilevel+1))cr (idim)=cr (idim)-m%ckey_max(ilevel+1)
        if(crr(idim)>=m%ckey_max(ilevel+1))crr(idim)=crr(idim)-m%ckey_max(ilevel+1)
     enddo
     ckey = pcs_index(cll,cl,cr,crr)
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

    if(action_part==action_kick_drift)then
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
       ! Build PCS 1D weights and indices at x_mid
       call pcs_1d_weights(x_mid, wll, wl, wr, wrr, cll, cl, cr, crr)
       vol2 = pcs_weight(wll, wl, wr, wrr)
       do idim=1,ndim
          if(cll(idim)<0)cll(idim)=cll(idim)+m%ckey_max(ilevel+1)
          if(cl (idim)<0)cl (idim)=cl (idim)+m%ckey_max(ilevel+1)
          if(cr (idim)>=m%ckey_max(ilevel+1))cr (idim)=cr (idim)-m%ckey_max(ilevel+1)
          if(crr(idim)>=m%ckey_max(ilevel+1))crr(idim)=crr(idim)-m%ckey_max(ilevel+1)
       enddo
       ckey2 = pcs_index(cll,cl,cr,crr)
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

!################################################################
! Weight computation subroutines for CIC and TSC interpolation
!################################################################
subroutine cic_weights_and_derivs(x, w1d, dw1d, il_out, ir_out, vol_out)
  !--------------------------------------------------------------
  ! Compute 1D CIC weights, derivatives, and cell indices from position.
  ! x(1:ndim) is the cell-centered coordinate.
  ! w1d(dim,1) = left weight (dl), w1d(dim,2) = right weight (dr)
  ! dw1d(dim,1) = -1, dw1d(dim,2) = +1 (CIC derivatives are constant)
  ! il_out, ir_out = cell indices (optional)
  ! vol_out(1:twotondim) = 3D volume weights (optional)
  !--------------------------------------------------------------
  use amr_parameters, only: ndim, twotondim
  implicit none
  real(kind=8),intent(in)::x(1:ndim)
  real(kind=8),intent(out)::w1d(1:ndim,1:2),dw1d(1:ndim,1:2)
  integer,intent(out),optional::il_out(1:ndim),ir_out(1:ndim)
  real(kind=8),intent(out),optional::vol_out(1:twotondim)
  integer::idim,ind,ix,iy,iz
  real(kind=8)::xd

  do idim=1,ndim
     xd = x(idim) + 0.5d0
     if(present(ir_out)) ir_out(idim) = int(xd)
     xd = xd - int(xd)  ! fractional offset in [0,1)
     w1d(idim,1) = 1.d0 - xd   ! left weight (dl)
     w1d(idim,2) = xd          ! right weight (dr)
     dw1d(idim,1) = -1.d0      ! derivative of left weight
     dw1d(idim,2) = +1.d0      ! derivative of right weight
     if(present(il_out).and.present(ir_out)) il_out(idim) = ir_out(idim) - 1
  end do
  ! 3D volume weights (optional output)
  if(present(vol_out))then
     ind = 0
     do iz=1,2
        do iy=1,2
           do ix=1,2
              ind = ind+1
              vol_out(ind) = w1d(1,ix)*w1d(2,iy)*w1d(3,iz)
           end do
        end do
     end do
  end if
end subroutine cic_weights_and_derivs

subroutine tsc_weights_and_derivs(x, w1d, dw1d, il_out, ic_out, ir_out, vol_out)
  !--------------------------------------------------------------
  ! Compute 1D TSC weights, derivatives, and cell indices from position.
  ! x(1:ndim) is the cell-centered coordinate.
  ! w1d(dim,1:3) = left, center, right weights
  ! dw1d(dim,1:3) = corresponding derivatives
  ! il_out, ic_out, ir_out = cell indices (optional)
  ! vol_out(1:threetondim) = 3D volume weights (optional)
  !--------------------------------------------------------------
  use amr_parameters, only: ndim, threetondim
  implicit none
  real(kind=8),intent(in)::x(1:ndim)
  real(kind=8),intent(out)::w1d(1:ndim,1:3),dw1d(1:ndim,1:3)
  integer,intent(out),optional::il_out(1:ndim),ic_out(1:ndim),ir_out(1:ndim)
  real(kind=8),intent(out),optional::vol_out(1:threetondim)
  integer::idim,il,ic,ir,ind,ix,iy,iz
  real(kind=8)::xd,xl,xc,xr

  do idim=1,ndim
     xd = x(idim)
     il = int(xd)-1
     ic = int(xd)
     ir = int(xd)+1
     xl = dble(il)+0.5d0
     xc = dble(ic)+0.5d0
     xr = dble(ir)+0.5d0
     ! TSC weights
     w1d(idim,1) = 0.5d0*(1.5d0-abs(xd-xl))**2
     w1d(idim,2) = 0.75d0-(xd-xc)**2
     w1d(idim,3) = 0.5d0*(1.5d0-abs(xd-xr))**2
     ! TSC derivatives
     dw1d(idim,1) = -(1.5d0-abs(xd-xl))*sign(1.d0,xd-xl)
     dw1d(idim,2) = -2.d0*(xd-xc)
     dw1d(idim,3) = -(1.5d0-abs(xd-xr))*sign(1.d0,xd-xr)
     ! Cell indices (optional output)
     if(present(il_out)) il_out(idim) = il
     if(present(ic_out)) ic_out(idim) = ic
     if(present(ir_out)) ir_out(idim) = ir
  end do
  ! 3D volume weights (optional output)
  if(present(vol_out))then
     ind = 0
     do iz=1,3
        do iy=1,3
           do ix=1,3
              ind = ind+1
              vol_out(ind) = w1d(1,ix)*w1d(2,iy)*w1d(3,iz)
           end do
        end do
     end do
  end if
end subroutine tsc_weights_and_derivs

subroutine pcs_1d_weights(x, wll_out, wl_out, wr_out, wrr_out, cll_out, cl_out, cr_out, crr_out)
  !--------------------------------------------------------------
  ! Compute 1D PCS weights and cell indices from position.
  ! x(1:ndim) is the cell-centered coordinate.
  ! wll_out, wl_out, wr_out, wrr_out = 1D weights for each dimension
  ! cll_out, cl_out, cr_out, crr_out = cell indices for each dimension
  ! Use pcs_weight(wll,wl,wr,wrr) to get 3D volume weights.
  ! Use pcs_index(cll,cl,cr,crr) to get neighbor index list (after periodic wrap).
  !--------------------------------------------------------------
  use amr_parameters, only: ndim
  implicit none
  real(kind=8),intent(in)::x(1:ndim)
  real(kind=8),intent(out)::wll_out(1:ndim),wl_out(1:ndim),wr_out(1:ndim),wrr_out(1:ndim)
  integer,intent(out)::cll_out(1:ndim),cl_out(1:ndim),cr_out(1:ndim),crr_out(1:ndim)
  integer::idim
  real(kind=8)::xll,xl,xr,xrr

  do idim=1,ndim
     crr_out(idim) = int(x(idim)+1.5d0)
     cr_out(idim)  = crr_out(idim)-1
     cl_out(idim)  = crr_out(idim)-2
     cll_out(idim) = crr_out(idim)-3
     xll = dble(cll_out(idim))+0.5d0
     xl  = dble(cl_out(idim))+0.5d0
     xr  = dble(cr_out(idim))+0.5d0
     xrr = dble(crr_out(idim))+0.5d0
     wll_out(idim) = (2.d0-abs(x(idim)-xll))**3/6.d0
     wl_out(idim)  = (4.d0-6.d0*(x(idim)-xl)**2+3.d0*abs(x(idim)-xl)**3)/6.d0
     wr_out(idim)  = (4.d0-6.d0*(x(idim)-xr)**2+3.d0*abs(x(idim)-xr)**3)/6.d0
     wrr_out(idim) = (2.d0-abs(x(idim)-xrr))**3/6.d0
  end do
end subroutine pcs_1d_weights

!################################################################
! Gradient computation subroutines for CIC and TSC interpolation
!################################################################
subroutine compute_gradient_cic_scalar(w1d, dw1d, field, grad_out)
  !--------------------------------------------------------------
  ! Compute gradient of a scalar field using CIC (2x2x2) weights.
  ! grad_out(k) = sum_ind dw/dx_k * field(ind)
  !--------------------------------------------------------------
  use amr_parameters, only: ndim, twotondim
  implicit none
  real(kind=8),intent(in)::w1d(1:ndim,1:2),dw1d(1:ndim,1:2)
  real(kind=8),intent(in)::field(1:twotondim)
  real(kind=8),intent(out)::grad_out(1:ndim)
  integer::ix,iy,iz,ind

  grad_out = 0.d0
  ind = 0
  do iz=1,2
     do iy=1,2
        do ix=1,2
           ind = ind+1
           grad_out(1) = grad_out(1) + dw1d(1,ix)*w1d(2,iy)*w1d(3,iz)*field(ind)
           grad_out(2) = grad_out(2) + w1d(1,ix)*dw1d(2,iy)*w1d(3,iz)*field(ind)
           grad_out(3) = grad_out(3) + w1d(1,ix)*w1d(2,iy)*dw1d(3,iz)*field(ind)
        end do
     end do
  end do
end subroutine compute_gradient_cic_scalar

subroutine compute_gradient_cic(w1d, dw1d, field, grad_out)
  !--------------------------------------------------------------
  ! Compute gradient of a vector field using CIC (2x2x2) weights.
  ! grad_out(k) = sum_ind dw/dx_k * field(k,ind)
  ! Each dimension k uses its own field values field(k,:).
  !--------------------------------------------------------------
  use amr_parameters, only: ndim, twotondim
  implicit none
  real(kind=8),intent(in)::w1d(1:ndim,1:2),dw1d(1:ndim,1:2)
  real(kind=8),intent(in)::field(1:ndim,1:twotondim)
  real(kind=8),intent(out)::grad_out(1:ndim)
  integer::ix,iy,iz,ind

  grad_out = 0.d0
  ind = 0
  do iz=1,2
     do iy=1,2
        do ix=1,2
           ind = ind+1
           grad_out(1) = grad_out(1) + dw1d(1,ix)*w1d(2,iy)*w1d(3,iz)*field(1,ind)
           grad_out(2) = grad_out(2) + w1d(1,ix)*dw1d(2,iy)*w1d(3,iz)*field(2,ind)
           grad_out(3) = grad_out(3) + w1d(1,ix)*w1d(2,iy)*dw1d(3,iz)*field(3,ind)
        end do
     end do
  end do
end subroutine compute_gradient_cic

subroutine compute_gradient_tsc_scalar(w1d, dw1d, field, grad_out)
  !--------------------------------------------------------------
  ! Compute gradient of a scalar field using TSC (3x3x3) weights.
  ! grad_out(k) = sum_ind dw/dx_k * field(ind)
  !--------------------------------------------------------------
  use amr_parameters, only: ndim, threetondim
  implicit none
  real(kind=8),intent(in)::w1d(1:ndim,1:3),dw1d(1:ndim,1:3)
  real(kind=8),intent(in)::field(1:threetondim)
  real(kind=8),intent(out)::grad_out(1:ndim)
  integer::ix,iy,iz,ind

  grad_out = 0.d0
  ind = 0
  do iz=1,3
     do iy=1,3
        do ix=1,3
           ind = ind+1
           grad_out(1) = grad_out(1) + dw1d(1,ix)*w1d(2,iy)*w1d(3,iz)*field(ind)
           grad_out(2) = grad_out(2) + w1d(1,ix)*dw1d(2,iy)*w1d(3,iz)*field(ind)
           grad_out(3) = grad_out(3) + w1d(1,ix)*w1d(2,iy)*dw1d(3,iz)*field(ind)
        end do
     end do
  end do
end subroutine compute_gradient_tsc_scalar

subroutine compute_gradient_tsc(w1d, dw1d, field, grad_out)
  !--------------------------------------------------------------
  ! Compute gradient of a vector field using TSC (3x3x3) weights.
  ! grad_out(k) = sum_ind dw/dx_k * field(k,ind)
  ! Each dimension k uses its own field values field(k,:).
  !--------------------------------------------------------------
  use amr_parameters, only: ndim, threetondim
  implicit none
  real(kind=8),intent(in)::w1d(1:ndim,1:3),dw1d(1:ndim,1:3)
  real(kind=8),intent(in)::field(1:ndim,1:threetondim)
  real(kind=8),intent(out)::grad_out(1:ndim)
  integer::ix,iy,iz,ind

  grad_out = 0.d0
  ind = 0
  do iz=1,3
     do iy=1,3
        do ix=1,3
           ind = ind+1
           grad_out(1) = grad_out(1) + dw1d(1,ix)*w1d(2,iy)*w1d(3,iz)*field(1,ind)
           grad_out(2) = grad_out(2) + w1d(1,ix)*dw1d(2,iy)*w1d(3,iz)*field(2,ind)
           grad_out(3) = grad_out(3) + w1d(1,ix)*w1d(2,iy)*dw1d(3,iz)*field(3,ind)
        end do
     end do
  end do
end subroutine compute_gradient_tsc

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
! Subroutines for Anisotropic Sub-Grid Noise Factors (Epstein Drag)
!#########################################################################
!#########################################################################
subroutine get_g_factors(dt, w0, a, nu0, g_par, g_perp)
  implicit none
  real(kind=8), intent(in) :: dt, w0, a, nu0
  real(kind=8), intent(out) :: g_par, g_perp
  
  real(kind=8) :: wt, kt, k0, k_inv_a, prefactor, w_thresh
  real(kind=8) :: u_start, u_end
  real(kind=8) :: val_base, val_tail
  real(kind=8) :: g_par_sq, g_perp_sq
  real(kind=8) :: zeta0, t_thresh_calc
  
  ! Constants
  k_inv_a = 1.0d0/a
  prefactor = (a**2)/nu0
  w_thresh = 1.0d-12 * a 
  
  ! Determine wt
  wt = w_det(dt, w0, a, nu0)
  
  ! Parallel Component (Analytic Formula)
  g_par = get_g_parallel_analytic(dt, w0, a, nu0)
  
  ! Perpendicular Component (Log K-space integral)
  if (wt < w_thresh .and. dt > 1.0d0) then
     zeta0 = asinh(a/w0)
     t_thresh_calc = (asinh(a/w_thresh) - zeta0) / nu0
     if (t_thresh_calc < 0.0d0) t_thresh_calc = 0.0d0
     
     kt = 1.0d0/w_thresh
     k0 = 1.0d0/w0
     u_start = log(k0)
     u_end = log(kt)
     
     val_base = integrate_g_perp_log(u_start, u_end, kt, k_inv_a) * prefactor
     val_tail = (dt - t_thresh_calc) ! Slope is 1 in normalized units
     
     g_perp_sq = val_base + val_tail
  else
     if (wt <= 1.0d-20) then
        kt = 1.0d20
     else
        kt = 1.0d0/wt
     end if
     k0 = 1.0d0/w0
     u_start = log(k0)
     u_end = log(kt)
     
     g_perp_sq = integrate_g_perp_log(u_start, u_end, kt, k_inv_a) * prefactor
  end if
  
  if (g_perp_sq < 0.0d0) g_perp_sq = 0.0d0
  g_perp = sqrt(g_perp_sq)/sqrt(dt)

end subroutine get_g_factors

function get_g_parallel_analytic(t, w0, a, nu0) result(val)
  implicit none
  real(kind=8), intent(in) :: t, w0, a, nu0
  real(kind=8) :: val
  real(kind=8) :: Pi_val, dtau, Q_val, Q2
  real(kind=8) :: term_E_dtau, term_E_2dtau
  real(kind=8) :: inner_most, middle_brack, outer_brack
  real(kind=8) :: numerator_main, denom_main
  
  if (t == 0.0d0) then
     val = 0.0d0
     return
  end if
  
  Pi_val = w0/a
  dtau = nu0 * t
  
  Q_val = sqrt(1.0d0 + Pi_val**(-2)) + 1.0d0/Pi_val
  Q2 = Q_val**2
  
  term_E_dtau = exp(dtau)
  term_E_2dtau = exp(2.0d0*dtau)
  
  inner_most = 4.0d0 + (4.0d0 + (-3.0d0 + 2.0d0*dtau)*term_E_dtau) * Q2
  middle_brack = -8.0d0*dtau + term_E_dtau*inner_most - Q2
  
  outer_brack = 1.0d0 + middle_brack * Q2
  
  numerator_main = 3.0d0 + 2.0d0*dtau + term_E_dtau*(-4.0d0*(1.0d0+Q2) + term_E_dtau*outer_brack)
  
  denom_main = dtau * (-1.0d0 + term_E_2dtau * Q2)**2
  
  if (numerator_main < 0.0d0) numerator_main = 0.0d0
  
  if (denom_main == 0.0d0) then
       val = 0.0d0
  else
       val = sqrt(numerator_main / denom_main) / sqrt(2.0d0)
  end if
  
end function get_g_parallel_analytic

function w_det(t, w0, a, nu0) result(w)
  implicit none
  real(kind=8), intent(in) :: t, w0, a, nu0
  real(kind=8) :: w
  real(kind=8) :: arg
  if (t == 0.0d0) then
     w = w0
     return
  end if
  arg = asinh(a/w0) + nu0*t
  if (arg > 700.0d0) then
     w = 0.0d0
  else
     w = a / sinh(arg)
  end if
end function w_det

function integrate_g_par(t_end, w0, a, nu0) result(integral)
  implicit none
  real(kind=8), intent(in) :: t_end, w0, a, nu0
  real(kind=8) :: integral
  real(kind=8) :: t_mid, dt_half, s
  real(kind=8) :: wt, ws, ratio
  integer :: i
  real(kind=8), dimension(10) :: x_gaus = (/ &
       0.076526521133497d0, 0.227785851141645d0, 0.373706088715419d0, &
       0.510867001950827d0, 0.636053680726515d0, 0.746331906460150d0, &
       0.839116971822218d0, 0.912234428251325d0, 0.963971927277913d0, &
       0.993128599185094d0 /)
  real(kind=8), dimension(10) :: w_gaus = (/ &
       0.152753387130725d0, 0.149172986472603d0, 0.142096109318382d0, &
       0.131688638449176d0, 0.118194531961518d0, 0.101930119817240d0, &
       0.083276741576704d0, 0.062672048334109d0, 0.040601429800386d0, &
       0.017614007139152d0 /)

  wt = w_det(t_end, w0, a, nu0)
  t_mid = 0.5d0 * t_end
  dt_half = 0.5d0 * t_end
  
  integral = 0.0d0
  do i = 1, 10
     s = t_mid + dt_half * x_gaus(i)
     ws = w_det(s, w0, a, nu0)
     if (ws > 0.0d0) then
        ratio = wt/ws
     else
        ratio = 0.0d0
     end if
     integral = integral + w_gaus(i) * (1.0d0 - ratio)**2
     
     s = t_mid - dt_half * x_gaus(i)
     ws = w_det(s, w0, a, nu0)
     if (ws > 0.0d0) then
        ratio = wt/ws
     else
        ratio = 0.0d0
     end if
     integral = integral + w_gaus(i) * (1.0d0 - ratio)**2
  end do
  integral = integral * dt_half
end function integrate_g_par

function integrate_g_perp_log(u_start, u_end, kt, k_inv_a) result(integral)
  implicit none
  real(kind=8), intent(in) :: u_start, u_end, kt, k_inv_a
  real(kind=8) :: integral
  real(kind=8) :: u_mid, du_half, u, k_val, term_val
  integer :: i
  real(kind=8), dimension(10) :: x_gaus = (/ &
       0.076526521133497d0, 0.227785851141645d0, 0.373706088715419d0, &
       0.510867001950827d0, 0.636053680726515d0, 0.746331906460150d0, &
       0.839116971822218d0, 0.912234428251325d0, 0.963971927277913d0, &
       0.993128599185094d0 /)
  real(kind=8), dimension(10) :: w_gaus = (/ &
       0.152753387130725d0, 0.149172986472603d0, 0.142096109318382d0, &
       0.131688638449176d0, 0.118194531961518d0, 0.101930119817240d0, &
       0.083276741576704d0, 0.062672048334109d0, 0.040601429800386d0, &
       0.017614007139152d0 /)

  u_mid = 0.5d0 * (u_start + u_end)
  du_half = 0.5d0 * (u_end - u_start)
  
  integral = 0.0d0
  do i = 1, 10
     u = u_mid + du_half * x_gaus(i)
     k_val = exp(u)
     term_val = integrand_perp_log_eval(k_val, kt, k_inv_a)
     integral = integral + w_gaus(i) * term_val
     
     u = u_mid - du_half * x_gaus(i)
     k_val = exp(u)
     term_val = integrand_perp_log_eval(k_val, kt, k_inv_a)
     integral = integral + w_gaus(i) * term_val
  end do
  integral = integral * du_half
end function integrate_g_perp_log

function integrand_perp_log_eval(k_val, kt, k_inv_a) result(val)
  implicit none
  real(kind=8), intent(in) :: k_val, kt, k_inv_a
  real(kind=8) :: val
  real(kind=8) :: term_sqrt_k, term_sqrt_kt
  real(kind=8) :: diff_linear, diff_sqrt_num, diff_sqrt_den, diff_sqrt
  real(kind=8) :: num_minus_den, den, r_minus_1, log_term
  
  term_sqrt_k = sqrt(k_val**2 + k_inv_a**2)
  term_sqrt_kt = sqrt(kt**2 + k_inv_a**2)
  
  diff_linear = k_inv_a * (kt - k_val)
  diff_sqrt_num = (k_inv_a**2) * (kt - k_val) * (kt + k_val)
  diff_sqrt_den = kt * term_sqrt_k + k_val * term_sqrt_kt
  diff_sqrt = diff_sqrt_num / diff_sqrt_den
  
  num_minus_den = diff_linear + diff_sqrt
  den = k_val * (k_inv_a + term_sqrt_kt)
  
  r_minus_1 = num_minus_den / den
  
  if (abs(r_minus_1) < 1.0d-8) then
     log_term = r_minus_1 - 0.5d0*r_minus_1**2 + (1.0d0/3.0d0)*r_minus_1**3
  else
     log_term = log(1.0d0 + r_minus_1)
  end if
  val = term_sqrt_k * (log_term**2) * k_val
end function integrand_perp_log_eval

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module move_fine_module
