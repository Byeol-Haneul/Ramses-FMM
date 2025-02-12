module feedback_module

    type :: out_feedback_t
     real(kind=8)::mass
  end type out_feedback_t

contains
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
recursive subroutine r_thermal_feedback(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(out_feedback_t)::output,next_output

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_THERMAL_FEEDBACK,pst%iUpper+1,input_size,output_size,ilevel)
     call r_thermal_feedback(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
  else
     call thermal_feedback(pst%s,pst%s%star,ilevel,output%mass)
  endif

end subroutine r_thermal_feedback
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine thermal_feedback(s,p,ilevel,msn_loc)
  use amr_parameters, only: ndim,twotondim,dp
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use nbors_utils
  use cache_commons
  use cache
  use marshal, only: pack_fetch_refine,unpack_fetch_refine
  use boundaries, only: init_bound_refine
  use godunov_fine_module, only: init_flush_godunov,pack_flush_godunov,unpack_flush_godunov
  use hilbert
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  real(kind=8)::msn_loc
  !==================================================================
  ! This is the RAMSES routine for supernovae feedback using
  ! a thermal energy dump. The routine is called at every time step
  ! and modify hydro variables unew.
  ! Written by Romain Teyssier (mini-ramses version in June 2023).
  !==================================================================
  ! Local variables
  integer,dimension(1:ndim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_cell
  integer::i,ipart,icellp,ind,idim
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx_loc,vol_loc,vol_cell
  real(dp)::mejecta,dloss,dzloss,zloss,ekinetic,ethermal
  real(dp)::birth_time,t_sn,e_sn,dteff,dold
  type(oct),pointer::gridp
  type(msg_large_realdp)::dummy_large_realdp
  logical::ok_level,ok_leaf

#ifdef HYDRO
#if NDIM==3
  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Supernovae progenitors life time from Myr to proper time in code units
  t_SN=r%t_SNII*1d6*(365.*24.*3600.)/(scale_t/g%aexp**2)

  ! Supernovae specific energy from cgs to code units
  e_SN=r%E_SNII/(r%M_SNII*2d33)/scale_v**2

  ! Total mass of ejecta
  msn_loc=0d0

  ! Open cache for array uold (fetch) and unew (flush)
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_large_realdp)/32,&
                pack=pack_fetch_refine,unpack=unpack_fetch_refine,&
                init=init_flush_godunov, flush=pack_flush_godunov,&
                combine=unpack_flush_godunov, bound=init_bound_refine)

  ! Loop over particles in Hilbert order
  do ipart=p%headp(ilevel),p%tailp(ilevel)

     ! Compute time step for that particle
     dteff=g%dtnew(p%levelp(ipart))*g%aexp**2

     ! Select only recently formed stars
     birth_time=p%tp(ipart) ! Proper time
     if(birth_time.lt.(g%texp-t_SN-dteff))cycle ! Already exploded
     if(birth_time.ge.(g%texp-t_SN))cycle ! Not old enough

     ok_level=.true.

     ! Find parent cell at level ilevel
     do idim=1,ndim
        ckey(idim)=int(p%xp(ipart,idim)/dx_loc)
     end do

     ! Cell volume at level ilevel
     vol_cell=vol_loc

     ! Get parent cell at level ilevel using cache
     hash_cell(0)=ilevel+1
     hash_cell(1:ndim)=ckey(1:ndim)
     call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icellp,flush_cache=.true.,fetch_cache=.true.)

     ! If cell does not exist at current level, then find cell at coarser level
     if(.not.associated(gridp))then

        ! NGP at level ilevel-1
        do idim=1,ndim
           ckey(idim)=int(p%xp(ipart,idim)/dx_loc/2)
        end do

        ! Cell volume at level ilevel-1
        vol_cell=vol_loc*2**ndim

        ! Get parent cell at level ilevel-1 using cache
        hash_cell(0)=ilevel
        hash_cell(1:ndim)=ckey(1:ndim)
        call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icellp,flush_cache=.true.,fetch_cache=.true.)
        if(.not.associated(gridp))ok_level=.false.

     end if

     if(.not. ok_level)then
        write(*,*)"Something went wrong in thermal_feedback"
        write(*,*)"Current level grid and coarser grid both dont exist..."
        stop
     endif

     ok_leaf = .not. gridp%refined(icellp)
     if(.not. ok_leaf)then
        write(*,*)"Something went wrong in thermal_feedback"
        write(*,*)"Cell should be a leaf cell..."
        stop
     endif

     ! Compute supernova properties
     mejecta=r%eta_SNII*p%mp(ipart)
     dloss=mejecta/vol_cell
     ethermal=dloss*e_SN
     ekinetic=dloss*0.5d0*(p%vp(ipart,1)**2+p%vp(ipart,2)**2+p%vp(ipart,3)**2)
     zloss=r%yield_SNII+(1d0-r%yield_SNII)*p%zp(ipart)
     dzloss=dloss*zloss

     ! Update unew
     gridp%unew(icellp,1)=gridp%unew(icellp,1)+dloss
     gridp%unew(icellp,2)=gridp%unew(icellp,2)+dloss*p%vp(ipart,1)
     gridp%unew(icellp,3)=gridp%unew(icellp,3)+dloss*p%vp(ipart,2)
     gridp%unew(icellp,4)=gridp%unew(icellp,4)+dloss*p%vp(ipart,3)
     gridp%unew(icellp,5)=gridp%unew(icellp,5)+ekinetic+ethermal
     if(r%metal)gridp%unew(icellp,r%imetal)=gridp%unew(icellp,r%imetal)+dzloss

     ! If dual energy scheme is activated, update entropy
     if(r%entropy.and.r%dual_energy.GE.0)then
        dold = gridp%uold(icellp,1)
        gridp%unew(icellp,r%ientropy)=gridp%unew(icellp,r%ientropy)+ethermal/dold**(r%gamma-1)*(r%gamma-1)
     endif

     ! Update particle mass
     p%mp(ipart)=p%mp(ipart)-mejecta

     ! Update total mass of ejecta
     msn_loc=msn_loc+mejecta

  end do
  ! End loop over particles

  call close_cache(s,m%grid_dict)

end associate
#endif  
#endif
end subroutine thermal_feedback
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine m_mechanical_feedback(pst,ilevel,mass_fbk)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::ilevel
  real(kind=8)::mass_fbk
  ! This is the master routine calling mechanical feedback.
  ! Note that it loops over all levels that are time-synchronized.
  integer::i
  type(out_feedback_t)::output_fbk

  mass_fbk=0d0
  do i=ilevel,pst%s%r%levelmin,-1
     call r_mechanical_feedback(pst,i,1,output_fbk,2)
     mass_fbk=mass_fbk+output_fbk%mass
     if(i>pst%s%r%levelmin)then
        if(pst%s%r%nsubcycle(i-1)==2.and.pst%s%g%isubcycle(i)==1)exit
     endif
  end do

end subroutine m_mechanical_feedback
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
recursive subroutine r_mechanical_feedback(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(out_feedback_t)::output,next_output

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_MECHANICAL_FEEDBACK,pst%iUpper+1,input_size,output_size,ilevel)
     call r_mechanical_feedback(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
  else
     call mechanical_feedback(pst%s,pst%s%star,ilevel,output%mass)
  endif

end subroutine r_mechanical_feedback
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine mechanical_feedback(s,p,ilevel,msn_loc)
  use amr_parameters, only: ndim,twotondim,dp
  use amr_commons, only: nbor,oct
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use nbors_utils
  use cache_commons
  use cache
  use marshal, only: pack_fetch_refine,unpack_fetch_refine
  use boundaries, only: init_bound_refine
  use godunov_fine_module, only: init_flush_godunov,pack_flush_godunov,unpack_flush_godunov
  use hilbert
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  real(kind=8)::msn_loc
  !==================================================================
  ! This is the RAMSES routine for supernovae feedback using
  ! only mechanical energy, no thermal energy.
  ! The routine modifies hydro variables unew.
  ! Written by Taysun Kimm (mini-ramses version in June 2023).
  !==================================================================
  ! Parameters
  ! Number of neighboring cells to deposit mass/momentum/energy
  integer, parameter::nSNnei=48
  ! Number of cells corresponding to the central cell to deposit mass
  real(dp),parameter::nSNcen=4
  ! Momentum input
  ! p_sn = A_SN*nH**(alpha)*ESN**(beta)*ZpSN**(gamma)
  ! ex) Thornton et al.
  !     A_SN = 3e5, alphaN = -2/17, beta = 16/17, gamma = -0.14
  ! ex) Kim & Ostriker (2015) uniform case
  !     A_SN = 2.17e5, alpha = -0.13, beta = 0.93
  real(dp),parameter::A_SN=3d5
  real(dp),parameter::expN_SN=-2d0/17d0
  real(dp),parameter::expE_SN=+16d0/17d0
  real(dp),parameter::expZ_SN=-0.14

  real(dp)::d,d_nei,dloss,ekloss,dzloss,zloss,dm_ejecta,dm_load,m_SN
  real(dp)::e,ekk,eth,p_solid,ek_solid,f_esn2,f_w_cell,f_w_crit
  real(dp)::nH_nei,u,v,w,up,vp,wp,T2,x,y,z,rr,Z_nei,Zdepen
  real(dp)::vload,vload_rad,vol_nei
  real(dp),dimension(1:3,1:nSNnei)::xSNnei
  real(dp),dimension(1:3,1:nSNnei)::vSNnei
  real(dp)::f_LOAD,f_CANCEL,f_ESN,f_LOAD_CEN
  integer,dimension(1:ndim)::ckey,ckey_ref,ckey_nbor
  integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
  integer::i,j,k,ipart,icellp,icelln,ind,idim,ivar,ipart_ref
  integer,dimension(1:ndim)::ix
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx_loc,vol_loc,vol_cell
  real(dp)::mejecta,ekinetic,ethermal
  real(dp)::birth_time,t_sn,e_sn,dteff,dold,num_SN
  real(dp),dimension(1:3)::xcen,xnei
  real(dp),dimension(1:nvar)::q
  integer,dimension(1:nSNnei)::icell_nbor,level_nbor
  type(nbor),dimension(1:nSNnei)::grid_nbor
  type(oct),pointer::gridp,gridn
  type(msg_large_realdp)::dummy_large_realdp
  logical::ok_level,ok_leaf,ok


#ifdef HYDRO
#if NDIM==3
  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Mechanical feedback parameters
  f_LOAD = nSNnei / dble(nSNcen + nSNnei) ! mass loading factor of the ejecta
  f_LOAD_CEN = 0.0 ! mass loading factor of the central cell
  f_CANCEL = 0.9387 ! correction due to momentum cancellation.
  f_ESN = 0.676 ! Blondin et al. (98) at t=trad

  ! Arrays to define neighbors (center=[0,0,0])
  ! normalized to dx = 1 = size of the central leaf cell in which a SN particle sits
  ! from -0.75 to 0.75
  ind=0
  do k=1,4
     do j=1,4
        do i=1,4
           ok=.true.
           if((i==1.or.i==4).and.(j==1.or.j==4).and.(k==1.or.k==4)) ok=.false. ! edge
           if((i==2.or.i==3).and.(j==2.or.j==3).and.(k==2.or.k==3)) ok=.false. ! centre
           if(ok)then
              ind = ind+1
              x = (i-1)+0.5d0 - 2
              y = (j-1)+0.5d0 - 2
              z = (k-1)+0.5d0 - 2
              rr = sqrt(dble(x*x+y*y+z*z))
              xSNnei(1,ind) = x/2d0
              xSNnei(2,ind) = y/2d0
              xSNnei(3,ind) = z/2d0
              vSNnei(1,ind) = x/rr
              vSNnei(2,ind) = y/rr
              vSNnei(3,ind) = z/rr
           endif
        enddo
     enddo
  enddo

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Supernovae progenitors life time from Myr to proper time in code units
  t_SN=r%t_SNII*1d6*(365.*24.*3600.)/(scale_t/g%aexp**2)

  ! Supernovae specific energy from cgs to code units
  e_SN=r%E_SNII/(r%M_SNII*2d33)/scale_v**2

  ! Supernovae progenitor mass in code units
  m_SN=r%M_SNII*2d33/(scale_d*scale_l**3)

  ! Total mass of ejecta
  msn_loc=0d0

  ! Open cache for array uold (fetch) and unew (flush)
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_large_realdp)/32,&
                pack=pack_fetch_refine,unpack=unpack_fetch_refine,&
                init=init_flush_godunov, flush=pack_flush_godunov,&
                combine=unpack_flush_godunov, bound=init_bound_refine)

  ! Loop over particles in Hilbert order
  do ipart=p%headp(ilevel),p%tailp(ilevel)

     ! Compute time step for that particle
     dteff=g%dtnew(p%levelp(ipart))*g%aexp**2

     ! Select only recently formed stars
     birth_time=p%tp(ipart) ! Proper time
     if(birth_time.lt.(g%texp-t_SN-dteff))cycle ! Already exploded
     if(birth_time.ge.(g%texp-t_SN))cycle ! Not old enough

     ok_level=.true.

     ! Find parent cell at level ilevel
     do idim=1,ndim
        ckey(idim)=int(p%xp(ipart,idim)/dx_loc)
        xcen(idim)=ckey(idim)+0.5
     end do

     ! Cell volume at level ilevel
     vol_cell=vol_loc

     ! Get parent cell at level ilevel using cache
     hash_cell(0)=ilevel+1
     hash_cell(1:ndim)=ckey(1:ndim)
     call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icellp,flush_cache=.true.,fetch_cache=.true.)

     ! If cell does not exist at current level, then find cell at coarser level
     if(.not.associated(gridp))then

        ! NGP at level ilevel-1
        do idim=1,ndim
           ckey(idim)=int(p%xp(ipart,idim)/dx_loc/2)
           xcen(idim)=ckey(idim)+0.5
        end do

        ! Cell volume at level ilevel-1
        vol_cell=vol_loc*2**ndim

        ! Get parent cell at level ilevel-1 using cache
        hash_cell(0)=ilevel
        hash_cell(1:ndim)=ckey(1:ndim)
        call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icellp,flush_cache=.true.,fetch_cache=.true.)
        if(.not.associated(gridp))ok_level=.false.

     end if

     ! Lock grid in cache
     call lock_cache(s,gridp)

     if(.not. ok_level)then
        write(*,*)"Something went wrong in mechanical_feedback"
        write(*,*)"Current level grid and coarser grid both dont exist..."
        stop
     endif

     ok_leaf = .not. gridp%refined(icellp)
     if(.not. ok_leaf)then
        write(*,*)"Something went wrong in mechanical_feedback"
        write(*,*)"Cell should be a leaf cell..."
        stop
     endif

     ! Compute supernova properties
     mejecta=r%eta_SNII*p%mp(ipart)
     up=p%vp(ipart,1)
     vp=p%vp(ipart,2)
     wp=p%vp(ipart,3)
     dloss=mejecta/vol_cell
     ekloss=dloss*0.5d0*(up**2+vp**2+wp**2)
     zloss=r%yield_SNII+(1d0-r%yield_SNII)*p%zp(ipart)
     dzloss=dloss*zloss
     num_sn=mejecta/m_SN

     ! Compute central cell properties
     d=max(gridp%uold(icellp,1),r%smallr)
     u=gridp%uold(icellp,2)/d
     v=gridp%uold(icellp,3)/d
     w=gridp%uold(icellp,4)/d
     e=gridp%uold(icellp,5)
     ekk=0.5*d*(u**2+v**2+w**2)
     eth=e-ekk
     T2=eth/d*scale_T2*(r%gamma-1)
     if(r%metal)z=gridp%uold(icellp,r%imetal)/d

     ! Collect all neighboring cell from hash table
     do j=1,nSNnei

        ! Compute neighboring cell coordinates
        xnei(1:ndim)=xcen(1:ndim)+xSNnei(1:ndim,j)
        ! Periodic boundary conditions
        do idim=1,ndim
           if(xnei(idim)<                    0.0d0)xnei(idim)=xnei(idim)+m%ckey_max(hash_cell(0))
           if(xnei(idim)>=m%ckey_max(hash_cell(0)))xnei(idim)=xnei(idim)-m%ckey_max(hash_cell(0))
        end do

        ! Get neighboring cell at current level
        ckey_nbor(1:ndim)=int(xnei(1:ndim))
        hash_nbor(0)=hash_cell(0)
        hash_nbor(1:ndim)=ckey_nbor(1:ndim)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.)

        ! If missing, get neighboring cell at current level - 1
        if(.not.associated(gridn))then
           ckey_nbor(1:ndim)=int(xnei(1:ndim)/2.0)
           hash_nbor(0)=hash_cell(0)-1
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.)

        ! If refined, get neighboring cell at current level + 1
        else if (gridn%refined(icelln))then
           ckey_nbor(1:ndim)=int(xnei(1:ndim)*2.0)
           hash_nbor(0)=hash_cell(0)+1
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.)
        endif

        ! Lock grid in cache
        call lock_cache(s,gridn)

        grid_nbor(j)%p => gridn
        icell_nbor(j) = icelln
        level_nbor(j) = hash_nbor(0)-1

     end do

     ! Update unew in central cell
     gridp%unew(icellp,1)=gridp%unew(icellp,1)+dloss-dloss*f_LOAD-d*f_LOAD_CEN
     gridp%unew(icellp,2)=gridp%unew(icellp,2)+dloss*up-dloss*up*f_LOAD-d*u*f_LOAD_CEN
     gridp%unew(icellp,3)=gridp%unew(icellp,3)+dloss*vp-dloss*vp*f_LOAD-d*v*f_LOAD_CEN
     gridp%unew(icellp,4)=gridp%unew(icellp,4)+dloss*wp-dloss*wp*f_LOAD-d*w*f_LOAD_CEN
     gridp%unew(icellp,5)=gridp%unew(icellp,5)+ekloss-ekloss*f_LOAD-(ekk+eth)*f_LOAD_CEN

     ! Update metals
     if(r%metal)gridp%unew(icellp,r%imetal)=gridp%unew(icellp,r%imetal)+dzloss-dzloss*f_LOAD-d*z*f_LOAD_CEN

     ! Update passive scalars so that they don't change
     do ivar=6,nvar
        if(r%metal.and.ivar==r%imetal)cycle
        q(ivar)=gridp%uold(icellp,ivar)/max(gridp%uold(icellp,1),r%smallr)
        gridp%unew(icellp,ivar)=gridp%unew(icellp,ivar)+(dloss-dloss*f_LOAD-d*f_LOAD_CEN)*q(ivar)
     end do

     ! Update conservative variables in neighboring cells
     dm_ejecta = dloss*f_LOAD/dble(nSNnei)
     dm_load = dm_ejecta + d*f_LOAD_CEN/dble(nSNnei)

     ! Loop over solid angles
     do j=1,nSNnei

        ! Gather neighboring grid
        gridn => grid_nbor(j)%p
        icelln = icell_nbor(j)

        ! Neighboring cell properties
        vol_nei = dble(twotondim)**(hash_cell(0)-1-level_nbor(j))
        Z_nei = r%z_ave*0.02
        d_nei = max(gridn%uold(icelln,1),r%smallr)
        if(r%metal) Z_nei = gridn%uold(icelln,r%imetal)/d_nei

        ! Compute actual mass ratio
        f_w_cell = (dm_load+d_nei/8d0)/dm_ejecta - 1d0

        ! Compute critical mass ratio
        nH_nei = d_nei*scale_nH
        Zdepen = (max(0.01_dp,Z_nei/0.02))**(expZ_SN*2d0)
        f_w_crit = (A_SN/1d4)**2d0/(f_ESN*r%M_SNII)*num_sn**((expE_SN-1d0)*2d0)*nH_nei**(expN_SN*2d0)*Zdepen - 1d0
        f_w_crit = max(0d0,f_w_crit) ! safety

        ! Compute SN terminal momentum
        vload_rad = sqrt(2d0*f_ESN*e_SN*(1d0+f_w_crit))/f_CANCEL/(1d0+f_w_cell)/f_LOAD

        ! Determine in which phase is the blast wave
        if(f_w_cell.ge.f_w_crit)then ! radiative phase
           vload = vload_rad
        else ! adiabatic phase
           f_esn2 = 1d0-(1d0-f_ESN)*f_w_cell/f_w_crit
           vload = sqrt(2d0*f_esn2*e_SN/(1d0+f_w_cell))/f_LOAD
        endif
        if(vload > vload_rad) vload = vload_rad ! safety

        ! Compute radial momentun and kinetic energy
        p_solid = (1d0+f_w_cell)*dm_ejecta*vload
        ek_solid = p_solid*(vload*f_LOAD)/2d0

        ! Add mass, momentum and energy coming from central cell loading
        gridn%unew(icelln,1)=gridn%unew(icelln,1)+(dloss*f_LOAD+d*f_LOAD_CEN)/dble(nSNnei)/vol_nei
        gridn%unew(icelln,2)=gridn%unew(icelln,2)+(dloss*up*f_LOAD+d*u*f_LOAD_CEN)/dble(nSNnei)/vol_nei
        gridn%unew(icelln,3)=gridn%unew(icelln,3)+(dloss*vp*f_LOAD+d*v*f_LOAD_CEN)/dble(nSNnei)/vol_nei
        gridn%unew(icelln,4)=gridn%unew(icelln,4)+(dloss*wp*f_LOAD+d*w*f_LOAD_CEN)/dble(nSNnei)/vol_nei
        gridn%unew(icelln,5)=gridn%unew(icelln,5)+(ekloss*f_LOAD+(ekk+eth)*f_LOAD_CEN)/dble(nSNnei)/vol_nei

        ! Add metals
        if(r%metal)then
           gridn%unew(icelln,r%imetal)=gridn%unew(icelln,r%imetal)+(dzloss*f_LOAD+d*z*f_LOAD_CEN)/dble(nSNnei)/vol_nei
        endif

        ! Update passive scalars coming from central cell loading
        do ivar=6,nvar
           if(r%metal.and.ivar==r%imetal)cycle
           gridn%unew(icelln,ivar)=gridn%unew(icelln,ivar)+(dloss*f_LOAD+d*f_LOAD_CEN)*q(ivar)/dble(nSNnei)/vol_nei
        end do

        ! Add momentum and energy coming from the cold shell
        gridn%unew(icelln,2)=gridn%unew(icelln,2)+p_solid*vSNnei(1,j)/vol_nei
        gridn%unew(icelln,3)=gridn%unew(icelln,3)+p_solid*vSNnei(2,j)/vol_nei
        gridn%unew(icelln,4)=gridn%unew(icelln,4)+p_solid*vSNnei(3,j)/vol_nei
        gridn%unew(icelln,5)=gridn%unew(icelln,5)+ek_solid/vol_nei

     end do
     ! End loop over solid angle

     ! Unlock all octs
     call unlock_cache(s,gridp)
     do j=1,nSNnei
        gridn => grid_nbor(j)%p
        call unlock_cache(s,gridn)
     end do

     ! Update particle mass
     p%mp(ipart)=p%mp(ipart)-mejecta

     ! Update total mass of ejecta
     msn_loc=msn_loc+mejecta

  end do
  ! End loop over particles

  call close_cache(s,m%grid_dict)

end associate
#endif
#endif
end subroutine mechanical_feedback
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module feedback_module
