module sink_accretion_module

   type :: out_accretion_t
      real(kind=8)::mass
   end type out_accretion_t

contains
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################`
recursive subroutine r_sink_accretion(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(out_accretion_t)::output,next_output

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SINK_ACCRETION,pst%iUpper+1,input_size,output_size,ilevel)
     call r_sink_accretion(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
  else
     call sink_accretion(pst%s,pst%s%sink,ilevel,output%mass)
  endif

end subroutine r_sink_accretion
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine sink_accretion(s,p,ilevel,macc_loc)
  use constants
  use amr_parameters, only: ndim,twotondim,dp
  use amr_commons, only: nbor,oct
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t,cross
  use params_module
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
  real(kind=8)::macc_loc
  !==================================================================
  ! This is the RAMSES routine for sink (black hole) particle accretion.
  ! For now, it is focused on a simple mass-weighted Bondi-Hoyle-Lyttleton accretion scheme
  ! The routine modifies hydro variables unew, as well as sink particle properties.
  ! Written by Nicholas Choustikov (Feb 2025)
  !==================================================================
  ! Local variables
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,factG ! Units
  logical::ok
  real(dp)::rr,x,y,z,rrad
  integer::nBHnei,iBHnei
  real(dp),dimension(:,:),allocatable::xBHnei
  integer,dimension(:,:),allocatable::ckeynei
  real(dp),dimension(:),allocatable::vol
  real(dp)::dx_loc,vol_loc,vol_cell
  real(dp),dimension(1:ndim)::xcen,xnei,xrel
  integer,dimension(1:ndim)::ckey,ckey_nbor
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::i,j,k,ipart,icelln,ind,idim,nrad
  real(dp)::d,e,ethermal,r2_sink,v_bondi,cs_gas,cs,rho_gas,velocity
  real(dp)::weight
  real(dp),dimension(1:ndim)::vv,v_rel,x_acc,p_acc,l_acc,vel_gas
  type(oct),pointer::gridn
  real(dp)::dMBH_overdt,dMEd_overdt,m_acc,d_acc,m_gas,bondi_mass
  real(dp)::rho_inf,weighted_bondi!,dMdt_freefall,t_ff
  type(msg_large_realdp)::dummy_large_realdp

#ifdef HYDRO
#if NDIM==3

  associate(r=>s%r,g=>s%g,m=>s%m)

  if(r%verbose)write(*,*)'Entering accrete_sink...'

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Get all units and cell sizes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Gravitational constant
  factG=1
  if(r%cosmo)factG=3d0/4d0/twopi*g%omega_m*g%aexp

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Prepare for the B-spline interpolation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

  ! Compute number of cells within B-spline region
  nBHnei = int(r%sink_b_spline_order**ndim)
  allocate(xBHnei(1:ndim,1:nBHnei),ckeynei(1:ndim,1:nBHnei),vol(1:nBHnei))
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Open Cache
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! Open cache for array uold (fetch) and unew (flush)
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_large_realdp)/32,&
                pack=pack_fetch_refine,unpack=unpack_fetch_refine,&
                init=init_flush_godunov, flush=pack_flush_godunov,&
                combine=unpack_flush_godunov, bound=init_bound_refine)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Begin loop over sink particles
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! Loop over particles in Hilbert order
  hash_nbor(0) = ilevel+1
  xBHnei=0d0; ckeynei=0d0; vol=0d0
  macc_loc=0d0
  do ipart = p%headp(ilevel), p%tailp(ilevel)
     
     ! Black hole position
     xcen(1:ndim) = p%xp(ipart,1:ndim) / dx_loc

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Initialise B-spline interpolation
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! TODO: In principle these weights can be computed once for every sink
     if      (r%sink_b_spline_order==2)then
        call sink_B_spline_weights_CIC(s,xcen(1:ndim),xBHnei,ckeynei,vol)
     else if (r%sink_b_spline_order==3)then
        call sink_B_spline_weights_TSC(s,xcen(1:ndim),xBHnei,ckeynei,vol)
     else if (r%sink_b_spline_order==4)then
        call sink_B_spline_weights_PCS(s,xcen(1:ndim),xBHnei,ckeynei,vol)
     else
        write(*,*)'This is an unknown B-spline order'
        ckeynei = -1 ! To cause a seg-fault
     end if

     ! Initialise sink information at zero
     rho_gas=0d0; vel_gas=0d0; cs_gas=0d0; m_gas=0d0; weighted_bondi=0d0

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Collect local gas information
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

     ! Loop over all cells in the accretion region, collecting physical properties and hash IDs
     do j = 1,nBHnei

        ! Get neighbouring cell coordinates
        xnei(1:ndim) = xBHnei(1:ndim,j)

        ! Periodic boundary conditions
        do idim=1,ndim
           if(xnei(idim)<                0.0d0)xnei(idim)=xnei(idim)+m%ckey_max(ilevel+1)
           if(xnei(idim)>=m%ckey_max(ilevel+1))xnei(idim)=xnei(idim)-m%ckey_max(ilevel+1)
        end do

        ! Get neighboring cell at current level
        hash_nbor(1:ndim)  = ckeynei(1:ndim,j)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.)

        ! If missing then cycle
        if(.not.associated(gridn))cycle

        ! Get the B-spline weights for this cell (they should already be normalised)
        weight = vol(j)

        ! Get physical information
        d                = max(gridn%uold(icelln,1),r%smallr)
        vv(1)            =     gridn%uold(icelln,2)/d
        vv(2)            =     gridn%uold(icelln,3)/d
        vv(3)            =     gridn%uold(icelln,4)/d
        e                =     gridn%uold(icelln,5)
        ethermal         = (e - 0.5d0*d*sum(vv(:)**2)) / d
        cs               = max((r%gamma-1.0d0)*ethermal,r%smallc**2)*r%acc_sink_boost**(-2d0/3d0)

        ! Add to average (weighted) information
        rho_gas          = rho_gas         + d          * weight
        vel_gas(1:ndim)  = vel_gas(1:ndim) + vv(1:ndim) * weight
        cs_gas           = cs_gas          + cs         * weight
        m_gas            = m_gas           + d          * weight * vol_loc

        ! Compute local bondi rate
        !!! Compute Bondi boosts here (e.g. c_s boost due to unresolved fluctuations)
        ! Perhaps we should model these unresolved fluctuations with a log-normal distribution, with uncertainty given by sigma_v (as in star-formation routines)
        if(r%use_local_bondi_rate)then
           v_rel(1:ndim) = vv(1:ndim) - p%vp(ipart,1:ndim)
           if(r%bondi_use_vrel)then
              v_bondi    = sqrt(sum(v_rel(:)**2) + cs**2)
           else
              v_bondi    = cs
           end if
           if(r%bondi_use_gas_mass)then
              bondi_mass = p%mp(ipart) + d*vol_loc ! can add a check for m_gas > M_sink
           else
              bondi_mass = p%mp(ipart)
           end if
           r2_sink     = (factG * bondi_mass / v_bondi**2)**2
           xrel(1:ndim) = xnei(1:ndim) - xcen(1:ndim) 
           rho_inf = d / (bondi_alpha(sqrt(sum(xrel(:)**2))/(r2_sink+tiny(0.0_dp))**0.5d0))
           dMBH_overdt = 4.0d0 * pi * rho_inf * r2_sink * v_bondi
           weighted_bondi = weighted_bondi + dMBH_overdt*weight
        end if
     end do

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Compute overall accretion rate
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

     

     !!! Compute BHL accretion rate
     if(r%use_local_bondi_rate)then
        dMBH_overdt = weighted_bondi
     else
        ! Bondi velocity
        v_rel(1:ndim) = vel_gas(1:ndim) - p%vp(ipart,1:ndim)
        if(r%bondi_use_vrel)then
           velocity   = sum(v_rel(:)**2)
           v_bondi    = sqrt(velocity + cs**2)
        else
           v_bondi    = cs
        end if

        ! Bondi mass
        if(r%bondi_use_gas_mass)then
           bondi_mass = p%mp(ipart) + m_gas ! can add a check for m_gas > M_sink
        else
           bondi_mass = p%mp(ipart)
        end if

        ! Bondi radius
        r2_sink     = (factG * bondi_mass / v_bondi**2)**2

        ! Density at infinity (using extrapolation)
        rho_inf = rho_gas / (bondi_alpha(dble(r%sink_b_spline_order)*0.5d0*dx_loc/(r2_sink+tiny(0.0_dp))**0.5d0))

        ! Bondi-Hoyle-Lyttleton accretion rate
        dMBH_overdt = 4.0d0 * pi * rho_inf * r2_sink * v_bondi
     end if

     ! Eddington accretion rate, which introduces an optional cap
     dMEd_overdt = 4.0d0 * pi * factG_in_cgs * p%mp(ipart) * mH / (0.1d0 * sigma_T * c_cgs) * scale_t
     if(r%eddington_cap>0)dMBH_overdt = min(dMBH_overdt, dMEd_overdt*r%eddington_cap)

     if(r%verbose_sink)then
        write(*,*)'Sink properties:',rho_gas,cs_gas,v_bondi,r2_sink
     end if

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Accrete from local cells
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

     ! Loop over all cells in the accretion region, proceeding with accretion
     m_acc=0.0d0; x_acc=0.0d0; p_acc=0.0d0; l_acc=0.0d0
     do j = 1, nBHnei

        ! Compute neighbouring cell coordinates
        xnei(1:ndim) = xBHnei(1:ndim,j)
        xrel(1:ndim) = xnei(1:ndim) - xcen(1:ndim)

        ! Periodic boundary conditions
        do idim=1,ndim
           if(xnei(idim)<                0.0d0)xnei(idim)=xnei(idim) + m%ckey_max(ilevel+1)
           if(xnei(idim)>=m%ckey_max(ilevel+1))xnei(idim)=xnei(idim) - m%ckey_max(ilevel+1)
        end do
        ! Get neighboring cell at current level
        hash_nbor(1:ndim)  = ckeynei(1:ndim,j)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.)
        ! If missing cycle
        if(.not.associated(gridn))cycle

        ! Get physical information
        d          = max(gridn%uold(icelln,1),r%smallr)
        vv(1:ndim) =     gridn%uold(icelln,2:ndim+1)/d
        e          =     gridn%uold(icelln,5)/d
        ! We need to remove all non-thermal energies as they are not accreted

        ! Get the weight for this cell based on the B-spline interpolation
        weight = vol(j)

        ! Get accreted mass for this cell
        d_acc = dMBH_overdt * g%dtnew(p%levelp(ipart)) * weight / vol_loc

        ! Ensure that the accreted amount is positive
        !d_acc = min(d_acc, 0.25d0 * d)
        d_acc = max(d_acc, 0.0_dp)

        ! Accrete from the cell
        gridn%unew(icelln,1)         = gridn%unew(icelln,1)      - d_acc
        do idim=1,ndim
           gridn%unew(icelln,idim+1) = gridn%unew(icelln,idim+1) - d_acc * vel_gas(idim)
        end do
        gridn%unew(icelln,5)         = gridn%unew(icelln,5)      - d_acc * e
        ! TODO: Add passive scalar accretion here

        !!! Accretion onto the black hole
        ! Accreted mass
        m_acc = m_acc + d_acc * vol_loc
        ! Accreted relative center of mass
        ! This should be zero, as it equals dM * sum_i (x_i - x_p)*w_i which is zero by construction
        x_acc(1:ndim) = x_acc(1:ndim) + d_acc * xrel(1:ndim) * vol_loc * dx_loc
        ! Accreted relative momentum
        p_acc(1:ndim) = p_acc(1:ndim) + d_acc * vel_gas(1:ndim)  * vol_loc
        ! Accreted relative angular momentum
        l_acc(1:ndim) = l_acc(1:ndim) + d_acc * cross(xrel(1:ndim), vel_gas(1:ndim)) * vol_loc * dx_loc
     end do ! End loop over j

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Add accreted quantities to sink
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

     ! Add accreted properties to sink variables
     p%xp(ipart,1:ndim) = ( p%mp(ipart) * p%xp(ipart,1:ndim) + x_acc(1:ndim) ) / ( p%mp(ipart) + m_acc )
     p%vp(ipart,1:ndim) = ( p%mp(ipart) * p%vp(ipart,1:ndim) + p_acc(1:ndim) ) / ( p%mp(ipart) + m_acc )
     p%jp(ipart,1:ndim) = ( p%mp(ipart) * p%jp(ipart,1:ndim) + l_acc(1:ndim) ) / ( p%mp(ipart) + m_acc )
     p%mp(ipart)        =   p%mp(ipart) + m_acc

     ! Save accreted mass to total
     macc_loc = macc_loc + m_acc

     if(r%verbose_sink)then
        write(*,*)'Accretion:',ipart,dMBH_overdt,dMEd_overdt,p%mp(ipart)
        write(*,*)'Accreted properties:',ipart,m_acc/g%dtnew(p%levelp(ipart)),sqrt(sum(x_acc(:)**2)),sqrt(sum(p_acc(:)**2)),sqrt(sum(l_acc(:)**2))
     end if

  end do ! End loop over ipart

  if(r%verbose_sink)write(*,*)'Total accreted mass:',macc_loc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Close cache
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call close_cache(s,m%grid_dict)

  deallocate(xBHnei,ckeynei,vol)

  end associate

contains
   ! Routine to return alpha, defined as rho/rho_inf, for a critical
   ! Bondi accretion solution. The argument is x = r / r_Bondi.
   ! This is from Krumholz et al. 2004 (AJC)
   REAL(dp) function bondi_alpha(x)
     implicit none
     REAL(dp) x
     REAL(dp), PARAMETER :: XMIN=0.01d0, XMAX=2.0d0
     INTEGER, PARAMETER :: NTABLE=51
     REAL(dp) lambda_c, xtable, xtablep1, alpha_exp
     integer idx
     !     Table of alpha values. These correspond to x values that run from
     !     0.01 to 2.0 with uniform logarithmic spacing. The reason for
     !     this choice of range is that the asymptotic expressions are
     !     accurate to better than 2% outside this range.
     REAL(dp), PARAMETER, DIMENSION(NTABLE) :: alphatable = (/ &
          820.254, 701.882, 600.752, 514.341, 440.497, 377.381, 323.427, &
          277.295, 237.845, 204.1, 175.23, 150.524, 129.377, 111.27, 95.7613, &
          82.4745, 71.0869, 61.3237, 52.9498, 45.7644, 39.5963, 34.2989, &
          29.7471, 25.8338, 22.4676, 19.5705, 17.0755, 14.9254, 13.0714, &
          11.4717, 10.0903, 8.89675, 7.86467, 6.97159, 6.19825, 5.52812, &
          4.94699, 4.44279, 4.00497, 3.6246, 3.29395, 3.00637, 2.75612, &
          2.53827, 2.34854, 2.18322, 2.03912, 1.91344, 1.80378, 1.70804, &
          1.62439 /)
     !     Define a constant that appears in these formulae
     lambda_c    = exp(1.5d0) / 4
     !     Deal with the off-the-table cases
     if (x .le. XMIN) then
        bondi_alpha = lambda_c / sqrt(2d0 * x**ndim)
     else if (x .ge. XMAX) then
        bondi_alpha = exp(1d0/x)
     else
        !     We are on the table
        idx = floor ((NTABLE-1) * log(x/XMIN) / log(XMAX/XMIN))
        xtable = exp(log(XMIN) + idx*log(XMAX/XMIN)/(NTABLE-1))
        xtablep1 = exp(log(XMIN) + (idx+1)*log(XMAX/XMIN)/(NTABLE-1d0))
        alpha_exp = log(x/xtable) / log(xtablep1/xtable)
        !     Note the extra +1s below because of fortran 1 offset arrays
        bondi_alpha = alphatable(idx+1) * (alphatable(idx+2)/alphatable(idx+1))**alpha_exp
     end if
   end function bondi_alpha

#endif
#endif
end subroutine sink_accretion

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine sink_B_spline_weights_PCS(s,x,xnei,ckey,vol)
   use amr_parameters, only: ndim,dp,fourtondim
   use ramses_commons, only: ramses_t
   implicit none
   type(ramses_t)::s
   real(dp),dimension(1:ndim)::x
   real(dp),dimension(1:ndim,1:fourtondim)::xnei
   integer,dimension(1:ndim,1:fourtondim)::ckey
   real(dp),dimension(1:fourtondim)::vol
   !==================================================================
   ! Simple routine to compute B-spline cells and weights for a given sink
   ! Nicholas Choustikov
   !==================================================================
   integer::idim,j
   integer,dimension(1:ndim)::crr,cr,cl,cll
   real(dp)::xrr,xr,xl,xll
   real(dp),dimension(1:ndim)::wrr,wr,wl,wll

   associate(r=>s%r,m=>s%m)

   ! PCS at level nlevelmax; a particle contributes to 4 cells in each direction
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
      if(cll(idim)<0)cll(idim)=m%ckey_max(r%nlevelmax+1)-1
      if(cl (idim)<0)cl (idim)=m%ckey_max(r%nlevelmax+1)-1
      if(cr (idim)==m%ckey_max(r%nlevelmax+1))cr (idim)=0
      if(crr(idim)==m%ckey_max(r%nlevelmax+1))crr(idim)=0
   enddo

   !!! Compute volumes, positions and cartesian keys
   ! NOTE: Sinks always exist at ndim = 3

   ! Compute cloud volumes
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

   ! Compute cartesian keys
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

   ! Compute neighbour positions
   do j = 1,fourtondim
      do idim = 1,ndim
         xnei(idim,j) = dble(ckey(idim,j)) + 0.5d0
      end do
   end do

   end associate

end subroutine sink_B_spline_weights_PCS

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine sink_B_spline_weights_TSC(s,x,xnei,ckey,vol)
   use amr_parameters, only: ndim,dp,threetondim
   use ramses_commons, only: ramses_t
   implicit none
   type(ramses_t)::s
   real(dp),dimension(1:ndim)::x
   real(dp),dimension(1:ndim,1:threetondim)::xnei
   integer,dimension(1:ndim,1:threetondim)::ckey
   real(dp),dimension(1:threetondim)::vol
   !==================================================================
   ! Simple routine to compute B-spline cells and weights for a given sink
   ! Nicholas Choustikov
   !==================================================================
   integer::idim,j
   integer,dimension(1:ndim)::cr,cl,cc
   real(dp)::xr,xl,xc
   real(dp),dimension(1:ndim)::wr,wl,wc

   associate(r=>s%r,m=>s%m)

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
      if(cl(idim)<0)cl(idim)=m%ckey_max(r%nlevelmax+1)-1
      if(cr(idim)==m%ckey_max(r%nlevelmax+1))cr(idim)=0
   enddo

   !!! Compute volumes, positions and cartesian keys
   ! NOTE: Sinks always exist at ndim = 3

   ! Compute cloud volumes
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

   ! Compute cartesian keys
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

   ! Compute neighbour positions
   do j = 1,threetondim
      do idim = 1,ndim
         xnei(idim,j) = dble(ckey(idim,j)) + 0.5d0
      end do
   end do

   end associate

end subroutine sink_B_spline_weights_TSC

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

subroutine sink_B_spline_weights_CIC(s,x,xnei,ckey,vol)
   use amr_parameters, only: ndim,dp,twotondim
   use ramses_commons, only: ramses_t
   implicit none
   type(ramses_t)::s
   real(dp),dimension(1:ndim)::x
   real(dp),dimension(1:ndim,1:twotondim)::xnei
   integer,dimension(1:ndim,1:twotondim)::ckey
   real(dp),dimension(1:twotondim)::vol
   !==================================================================
   ! Simple routine to compute B-spline cells and weights for a given sink
   ! Nicholas Choustikov
   !==================================================================
   integer::idim,j
   integer,dimension(1:ndim)::ir,il
   real(dp),dimension(1:ndim)::dr,dl

   associate(r=>s%r,m=>s%m)

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
      if(il(idim)<0)il(idim)=m%ckey_max(r%nlevelmax+1)-1
      if(ir(idim)==m%ckey_max(r%nlevelmax+1))ir(idim)=0
   enddo

   !!! Compute volumes, positions and cartesian keys
   ! NOTE: Sinks always exist at ndim = 3

   ! Compute cloud volumes
   vol(1)=dl(1)*dl(2)*dl(3)
   vol(2)=dr(1)*dl(2)*dl(3)
   vol(3)=dl(1)*dr(2)*dl(3)
   vol(4)=dr(1)*dr(2)*dl(3)
   vol(5)=dl(1)*dl(2)*dr(3)
   vol(6)=dr(1)*dl(2)*dr(3)
   vol(7)=dl(1)*dr(2)*dr(3)
   vol(8)=dr(1)*dr(2)*dr(3)

   ! Compute cartesian keys
   ckey(1:3,1)=(/il(1),il(2),il(3)/)
   ckey(1:3,2)=(/ir(1),il(2),il(3)/)
   ckey(1:3,3)=(/il(1),ir(2),il(3)/)
   ckey(1:3,4)=(/ir(1),ir(2),il(3)/)
   ckey(1:3,5)=(/il(1),il(2),ir(3)/)
   ckey(1:3,6)=(/ir(1),il(2),ir(3)/)
   ckey(1:3,7)=(/il(1),ir(2),ir(3)/)
   ckey(1:3,8)=(/ir(1),ir(2),ir(3)/)

   ! Compute neighbour positions
   do j = 1,twotondim
      do idim = 1,ndim
         xnei(idim,j) = dble(ckey(idim,j)) + 0.5d0
      end do
   end do

   end associate

end subroutine sink_B_spline_weights_CIC

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

end module sink_accretion_module
