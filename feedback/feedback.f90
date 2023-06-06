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
     call thermal_feedback(pst%s,pst%s%s,ilevel,output%mass)
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
  ! =================================================================
  ! This is the RAMSES routine for supernovae feedback using
  ! a thermal energy dump. The routine is called at every time step
  ! and modify hydro variables unew.
  ! Written by Romain Teyssier (mini-ramses version in June 2023).
  !==================================================================
  ! Local variables
  integer,dimension(1:ndim)::ckey
  integer(kind=8),dimension(0:ndim)::hash_cell
  integer::i,ipart,icell,ind,idim
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx_loc,vol_loc,vol_cell
  real(dp)::mejecta,mloss,mzloss,zloss,ekinetic,ethermal
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
     call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.true.)

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
        call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.true.)
        if(.not.associated(gridp))ok_level=.false.

     end if

     if(.not. ok_level)then
        write(*,*)"Something went wrong in thermal_feedback"
        write(*,*)"Current level grid and coarser grid both dont exist..."
        stop
     endif

     ok_leaf = .not. gridp%refined(icell)
     if(.not. ok_leaf)then
        write(*,*)"Something went wrong in thermal_feedback"
        write(*,*)"Cell should be a leaf cell..."
        stop
     endif

     ! Compute supernova properties
     mejecta=r%eta_SNII*p%mp(ipart)
     mloss=mejecta/vol_cell
     ethermal=mloss*e_SN
     ekinetic=mloss*0.5d0*(p%vp(ipart,1)**2+p%vp(ipart,2)**2+p%vp(ipart,3)**2)
     zloss=r%yield_SNII+(1d0-r%yield_SNII)*p%zp(ipart)
     mzloss=mloss*zloss

     ! Update unew
     gridp%unew(icell,1)=gridp%unew(icell,1)+mloss
     gridp%unew(icell,2)=gridp%unew(icell,2)+mloss*p%vp(ipart,1)
     gridp%unew(icell,3)=gridp%unew(icell,3)+mloss*p%vp(ipart,2)
     gridp%unew(icell,4)=gridp%unew(icell,4)+mloss*p%vp(ipart,3)
     gridp%unew(icell,5)=gridp%unew(icell,5)+ekinetic+ethermal
     if(r%metal)gridp%unew(icell,r%imetal)=gridp%unew(icell,r%imetal)+mzloss

     ! If dual energy scheme is activated, update entropy
     if(r%entropy.and.r%dual_energy.GE.0)then
        dold = gridp%uold(icell,1)
        gridp%unew(icell,r%ientropy)=gridp%unew(icell,r%ientropy)+ethermal/dold**(r%gamma-1)*(r%gamma-1)
     endif

!     write(*,'("feedback ",6(1PE13.6,1X))')g%texp,dteff,birth_time,(dold+mloss)*scale_nH,ethermal/(dold+mloss)*scale_v**2/1.38d-16*1.66d-24,gridp%uold(icell,5)/dold*scale_v**2/1.38d-16*1.66d-24

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
  !
  integer::i
  type(out_feedback_t)::output_fbk

  mass_fbk=0d0
  do i=ilevel,pst%s%r%levelmin,-1
     if(i>pst%s%r%levelmin.and.pst%s%r%nsubcycle(i-1)==2)then
        if(pst%s%g%isubcycle(i)==1)then
           call r_mechanical_feedback(pst,i,1,output_fbk,2)
           mass_fbk=mass_fbk+output_fbk%mass
           exit
        else
           cycle
        endif
     else
        call r_mechanical_feedback(pst,i,1,output_fbk,2)
        mass_fbk=mass_fbk+output_fbk%mass
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
     call mechanical_feedback(pst%s,pst%s%s,ilevel,output%mass)
  endif

end subroutine r_mechanical_feedback
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine mechanical_feedback(s,p,ilevel,msn_loc)
  use amr_parameters, only: ndim,twotondim,dp
  use hydro_parameters, only: nvar
  use amr_commons, only: nbor,oct
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
  use nbors_utils
  use cache_commons
  use cache
  use marshal, only: pack_fetch_refine,unpack_fetch_refine
  use boundaries, only: init_bound_refine
  use godunov_fine_module, only: init_flush_godunov,pack_flush_godunov,unpack_flush_godunov
  use rho_fine_module, only: sort_hilbert
  use hilbert
  implicit none
  type(ramses_t)::s
  type(part_t)::p
  integer::ilevel
  real(kind=8)::msn_loc
  ! =================================================================
  ! This is the RAMSES routine for supernovae feedback using
  ! only mechanical energy, no thermal energy.
  ! The routine is called every other time step and modify hydro variables unew.
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
  ! Local variables
  type(part_t)::sn
  real(dp)::d,d_nei,dloss,ekloss,dzloss,dm_ejecta,dm_load,m_SN
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
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx_loc,vol_loc
  real(dp)::mejecta,mloss,mzloss,zloss,ekinetic,ethermal
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

  ! Supernovae specific energy in code units
  e_SN=r%E_SNII/(r%M_SNII*2d33)/scale_v**2

  ! Supernovae progenitor mass in code units
  m_SN=r%M_SNII*2d33/(scale_d*scale_l**3)

  ! First collect supernovae stars into the current grid CPU domain 
  call collect_sn(s,p,sn,ilevel,msn_loc)
!  if(sn%npart>0)write(*,*)g%myid,' found ',sn%npart,' supernovae'
  
  ! Sort supernovae stars according to their cell Hilbert key
  if(sn%npart>0)allocate(sn%sortp(1:sn%npart),sn%workp(1:sn%npart))
  do i=1,sn%npart
     sn%sortp(i)=i
  end do
  ix=0
  if(sn%npart>0)call sort_hilbert(r,g,sn,1,sn%npart,ix,0,1,ilevel)

  ! Group supernovae particles if they have the same parent cell
  if(sn%npart>0)then
     sn%zp=sn%zp*sn%mp
     do idim=1,ndim
        sn%vp(1:sn%npart,idim)=sn%vp(1:sn%npart,idim)*sn%mp
     end do
     ckey_ref=-1
     ipart_ref=-1
     do i=1,sn%npart
        ipart=sn%sortp(i)
        ckey(1:ndim)=int(sn%xp(ipart,1:ndim)/dx_loc)
        if(ckey(1)==ckey_ref(1).AND.ckey(2)==ckey_ref(2).AND.ckey(3)==ckey_ref(3))then
           sn%mp(ipart_ref)=sn%mp(ipart_ref)+sn%mp(ipart)
           sn%zp(ipart_ref)=sn%zp(ipart_ref)+sn%zp(ipart)
           sn%vp(ipart_ref,1:ndim)=sn%vp(ipart_ref,1:ndim)+sn%vp(ipart,1:ndim)
           sn%sortp(i)=0
        else
           ckey_ref=ckey
           ipart_ref=ipart
        endif
     end do
     sn%zp=sn%zp/sn%mp
     do idim=1,ndim
        sn%vp(1:sn%npart,idim)=sn%vp(1:sn%npart,idim)/sn%mp
     end do
  endif

  ! Deposit momentum and kinetic energy of supernovae cell into neighboring cells.

  ! Open cache for array uold (fetch) and unew (flush)
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_large_realdp)/32,&
                pack=pack_fetch_refine,unpack=unpack_fetch_refine,&
                init=init_flush_godunov, flush=pack_flush_godunov,&
                combine=unpack_flush_godunov, bound=init_bound_refine)

  ! Loop over SN particles
  do i=1,sn%npart

     ipart=sn%sortp(i)
     if(ipart==0)cycle ! duplicate

     ! Set pointers to null
     icellp=0; icelln=0
     nullify(gridp)
     nullify(gridn)
     do j=1,nSNnei
        nullify(grid_nbor(j)%p)
     end do

     ! Find parent cell at level ilevel
     ckey(1:ndim)=int(sn%xp(ipart,1:ndim)/dx_loc)
     xcen(1:ndim)=ckey(1:ndim)+0.5

     ! Get parent cell at level ilevel using cache
     hash_cell(0)=ilevel+1
     hash_cell(1:ndim)=ckey(1:ndim)
     call get_parent_cell(s,hash_cell,m%grid_dict,gridp,icellp,flush_cache=.true.,fetch_cache=.true.,lock=.true.)

     ok_level = associated(gridp)
     if(.not. ok_level)then
        write(*,*)"Something went wrong in mechanical_feedback"
        write(*,*)"Current level grid should exist..."
        stop
     endif

     ok_leaf = .not. gridp%refined(icellp)
     if(.not. ok_leaf)then
        write(*,*)"Something went wrong in mechanical_feedback"
        write(*,*)"Current level cell should be a leaf cell..."
        stop
     endif

     ! Compute supernova properties
     mejecta=sn%mp(ipart)
     up=sn%vp(ipart,1)
     vp=sn%vp(ipart,2)
     wp=sn%vp(ipart,3)
     dloss=mejecta/vol_loc
     ekloss=dloss*0.5d0*(up**2+vp**2+wp**2)
     zloss=r%yield_SNII+(1d0-r%yield_SNII)*sn%zp(ipart)
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
           if(xnei(idim)<                0.0d0)xnei(idim)=xnei(idim)+m%ckey_max(ilevel+1)
           if(xnei(idim)>=m%ckey_max(ilevel+1))xnei(idim)=xnei(idim)-m%ckey_max(ilevel+1)
        end do

        ! Get neighboring cell at ilevel
        ckey_nbor(1:ndim)=int(xnei(1:ndim))
        hash_nbor(0)=ilevel+1
        hash_nbor(1:ndim)=ckey_nbor(1:ndim)
        call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.,lock=.true.)

        ! If missing, get neighboring cell at ilevel-1
        if(.not.associated(gridn))then
           call unlock_cache(s,gridn)
           ckey_nbor(1:ndim)=int(xnei(1:ndim)/2.0)
           hash_nbor(0)=ilevel
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.,lock=.true.)

        ! If refined, get neighboring cell at ilevel+1
        else if (gridn%refined(icelln))then
           call unlock_cache(s,gridn)
           ckey_nbor(1:ndim)=int(xnei(1:ndim)*2.0)
           hash_nbor(0)=ilevel+2
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridn,icelln,flush_cache=.true.,fetch_cache=.true.,lock=.true.)
        endif

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

!     write(*,'(2(I3,1X),15(1PE14.7,1X))')g%myid,ilevel,xcen(1:3),num_sn,d,gridp%unew(icellp,1),d*scale_nH,dloss,f_LOAD
     
     ! Update conservative variables in neighboring cells
     dm_ejecta = dloss*f_LOAD/dble(nSNnei)
     dm_load = dm_ejecta + d*f_LOAD_CEN/dble(nSNnei)

     ! Loop over solid angles
     do j=1,nSNnei

        ! Gather neighboring grid
        gridn => grid_nbor(j)%p
        icelln = icell_nbor(j)

        ! Neighboring cell properties
        vol_nei = dble(twotondim)**(ilevel-level_nbor(j))
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

!        write(*,'(3(I3,1X),15(1PE14.7,1X))')g%myid,j,level_nbor(j),xcen(1:3),xcen(1:ndim)+xSNnei(1:ndim,j),dble(gridn%ckey(1:ndim)),dble(icelln),num_sn,nH_nei,f_w_cell,f_w_crit,p_solid
!        write(*,'(3(I3,1X),15(1PE14.7,1X))')g%myid,j,level_nbor(j),xcen(1:3),vSNnei(1:ndim,j),num_sn,nH_nei,Z_nei,f_w_cell,f_w_crit,p_solid
        
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

  end do
  ! End loop over SN cells

  call close_cache(s,m%grid_dict)

  ! Dealocate supernovae particles
  if(sn%npart>0)deallocate(sn%xp,sn%vp,sn%mp,sn%zp,sn%sortp,sn%workp)
  sn%npart=0

  end associate
#endif  
#endif
end subroutine mechanical_feedback
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine collect_sn(s,p,sn,ilevel,msn_loc)
  use amr_parameters, only: ndim,twotondim,dp,nhilbert
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use pm_commons, only: part_t
#ifndef WITHOUTMPI
  use mpi
#endif
  use hilbert
  implicit none
  type(ramses_t)::s
  type(part_t)::p,sn
  integer::ilevel
  real(kind=8)::msn_loc

  integer::n_loc,n_tot,ipart,i,icpu,idim,info
  real(dp)::dteff,birth_time,dx_loc,t_SN
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer(kind=8),dimension(1:ndim)::ix
  integer(kind=8),dimension(1:nhilbert)::hk
  real(dp),dimension(:,:),allocatable::x_loc,v_loc
  real(dp),dimension(:),allocatable::m_loc,z_loc,x_tmp
  integer,dimension(:),allocatable::cpu_loc,ind_loc
#ifndef WITHOUTMPI
  integer::nbuffer,countrecv,countsend,tag=101,istart
  integer,dimension(MPI_STATUS_SIZE,s%g%ncpu)::statuses
  integer,dimension(:),allocatable::send_cnt,send_oft
  integer,dimension(:),allocatable::recv_cnt,recv_oft
  integer,dimension(s%g%ncpu)::reqsend,reqrecv
#endif

#ifdef HYDRO
#if NDIM==3
  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Supernovae progenitors life time from Myr to proper time in code units
  t_SN=r%t_SNII*1d6*(365.*24.*3600.)/(scale_t/g%aexp**2)

  ! Total mass of ejecta
  msn_loc=0d0

  ! Mesh spacing in that level
  dx_loc=r%boxlen/2**ilevel 

  ! Last time step
  if (ilevel==r%levelmin)then
     dteff = g%dtnew(ilevel)
  else
     if(r%nsubcycle(ilevel-1)==1)then
        dteff = g%dtnew(ilevel)
     else
        dteff = g%dtnew(ilevel)+g%dtold(ilevel)
     endif
  endif

  ! Select only recently formed stars
  n_loc=0
  do ipart=p%headp(ilevel),p%tailp(ilevel)
     birth_time=p%tp(ipart) ! Proper time
     if(birth_time.lt.(g%texp-t_SN-dteff))cycle ! Already exploded
     if(birth_time.ge.(g%texp-t_SN))cycle ! Not old enough
     n_loc=n_loc+1
     p%sortp(n_loc)=ipart
  end do

  ! Allocate local supernovae particles
  if(n_loc>0)then
     allocate(x_loc(1:n_loc,1:ndim),v_loc(1:n_loc,1:ndim))
     allocate(m_loc(1:n_loc),z_loc(1:n_loc))
#ifndef WITHOUTMPI
     allocate(ind_loc(1:n_loc),cpu_loc(1:n_loc))
#endif
  endif

  ! Loop over supernovae explosions
  do i=1,n_loc
     ! Get star particle index
     ipart=p%sortp(i)
#ifndef WITHOUTMPI
     ! Get parent grid Hilbert key
     ix(1:ndim)=int(p%xp(ipart,1:ndim)/(2*dx_loc))
     hk(1:nhilbert)=hilbert_key(ix,ilevel-1)
     ! Determine parent grid processor
     cpu_loc(i) = m%domain_hilbert(ilevel)%get_rank(hk)
     ind_loc(i) = i
#endif
     ! Store progenitor properties
     x_loc(i,1:ndim)=p%xp(ipart,1:ndim)
     v_loc(i,1:ndim)=p%vp(ipart,1:ndim)
     m_loc(i)=r%eta_SNII*p%mp(ipart)
     z_loc(i)=p%zp(ipart)
     ! Reduce star particle mass
     p%mp(ipart)=p%mp(ipart)-m_loc(i)
     ! Update total mass of ejecta
     msn_loc=msn_loc+m_loc(i)
  end do

#ifndef WITHOUTMPI
  ! Sort supernovae explosions according to their CPU
  if(n_loc>0)call quick_sort_int_int(cpu_loc,ind_loc,n_loc)

  ! Swap supernovae properties in sorted arrays
  if(n_loc>0)then
     allocate(x_tmp(1:n_loc))
     do idim=1,ndim
        x_tmp=x_loc(1:n_loc,idim)
        do i=1,n_loc
           x_loc(i,idim)=x_tmp(ind_loc(i))
        end do
        x_tmp=v_loc(1:n_loc,idim)
        do i=1,n_loc
           v_loc(i,idim)=x_tmp(ind_loc(i))
        end do
     end do
     x_tmp=m_loc
     do i=1,n_loc
        m_loc(i)=x_tmp(ind_loc(i))
     end do
     x_tmp=z_loc
     do i=1,n_loc
        z_loc(i)=x_tmp(ind_loc(i))
     end do
     deallocate(x_tmp)
     deallocate(ind_loc)
  endif
  
  ! Compute number of particles to send and receive
  allocate(send_cnt(1:g%ncpu))
  allocate(recv_cnt(1:g%ncpu))
  allocate(send_oft(1:g%ncpu))
  allocate(recv_oft(1:g%ncpu))
  send_cnt=0; recv_cnt=0
  do i=1,n_loc
     icpu=cpu_loc(i)
     send_cnt(icpu)=send_cnt(icpu)+1
  end do
  if(n_loc>0)deallocate(cpu_loc)
  call MPI_ALLTOALL(send_cnt(1),1,MPI_INTEGER,recv_cnt(1),1,MPI_INTEGER,MPI_COMM_WORLD,info)
  n_tot=SUM(recv_cnt)
  recv_oft=0; send_oft=0
  do icpu=2,g%ncpu
     recv_oft(icpu)=recv_oft(icpu-1)+recv_cnt(icpu-1)
     send_oft(icpu)=send_oft(icpu-1)+send_cnt(icpu-1)
  end do  
#else
  n_tot=n_loc
#endif

  ! Allocate supernovae that belong to the local CPU grid-domain  
  sn%npart=n_tot
  if(n_tot>0)then
     allocate(sn%xp(1:n_tot,1:ndim),sn%vp(1:n_tot,1:ndim))
     allocate(sn%mp(1:n_tot),sn%zp(1:n_tot))
  endif

#ifndef WITHOUTMPI
  ! Swap positions
  do idim=1,ndim
     
     countrecv=0
     do icpu=1,g%ncpu
        nbuffer=recv_cnt(icpu)
        if(nbuffer>0)then
           countrecv=countrecv+1
           istart=recv_oft(icpu)+1
           call MPI_IRECV(sn%xp(istart,idim),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
        endif
     end do
     
     countsend=0
     do icpu=1,g%ncpu
        nbuffer=send_cnt(icpu)
        if(nbuffer>0) then
           countsend=countsend+1
           istart=send_oft(icpu)+1
           call MPI_ISEND(x_loc(istart,idim),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
        end if
     end do
     
     call MPI_WAITALL(countrecv,reqrecv,statuses,info)
     call MPI_WAITALL(countsend,reqsend,statuses,info)
     
  end do
  
  ! Swap velocities
  do idim=1,ndim
     
     countrecv=0
     do icpu=1,g%ncpu
        nbuffer=recv_cnt(icpu)
        if(nbuffer>0)then
           countrecv=countrecv+1
           istart=recv_oft(icpu)+1
           call MPI_IRECV(sn%vp(istart,idim),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
        endif
     end do
     
     countsend=0
     do icpu=1,g%ncpu
        nbuffer=send_cnt(icpu)
        if(nbuffer>0) then
           countsend=countsend+1
           istart=send_oft(icpu)+1
           call MPI_ISEND(v_loc(istart,idim),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
        end if
     end do
     
     call MPI_WAITALL(countrecv,reqrecv,statuses,info)
     call MPI_WAITALL(countsend,reqsend,statuses,info)
     
  end do

  ! Swap masses
  countrecv=0
  do icpu=1,g%ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
        call MPI_IRECV(sn%mp(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do
  
  countsend=0
  do icpu=1,g%ncpu
     nbuffer=send_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=send_oft(icpu)+1
        call MPI_ISEND(m_loc(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do
  
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)
  call MPI_WAITALL(countsend,reqsend,statuses,info)
  
  ! Swap metallicities
  countrecv=0
  do icpu=1,g%ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
        call MPI_IRECV(sn%zp(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do
  
  countsend=0
  do icpu=1,g%ncpu
     nbuffer=send_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=send_oft(icpu)+1
        call MPI_ISEND(z_loc(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do
  
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  deallocate(send_cnt,send_oft,recv_cnt,recv_oft)
#else
  if(n_tot>0)then
     sn%xp=x_loc
     sn%vp=v_loc
     sn%mp=m_loc
     sn%zp=z_tot
  endif
#endif

  ! Deallocate local arrays
  if(n_loc>0)then
     deallocate(x_loc,v_loc,m_loc,z_loc)
  endif

  end associate
#endif  
#endif
end subroutine collect_sn
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module feedback_module
