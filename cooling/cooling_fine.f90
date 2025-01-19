module cooling_fine_module
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_cooling_fine(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_COOLING_FINE,pst%iUpper+1,input_size,0,ilevel)
     call r_cooling_fine(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call cooling_fine(pst%s%r,pst%s%g,pst%s%m,pst%s%cool,pst%s%tables,ilevel)
  endif

end subroutine r_cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cooling_fine(r,g,m,c,tables,ilevel)
  use constants
  use amr_parameters, only:dp,ndim,nvector,twotondim
  use hydro_parameters, only: nener
  use rt_parameters, only: nions,nrtgroups,smallnp
  use amr_commons, only:run_t,global_t,mesh_t
  use cooling_module, only:cooling_t,solve_cooling,T2_min_fix,set_table
  use coolrates_module, only: neq_cooling_t
  use rt_cooling_module, only: rt_solve_cooling
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(cooling_t)::c
  type(neq_cooling_t)::tables
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for leaf cells at level ilevel
  !-------------------------------------------------------------------
  integer::i,ii,ind,igrid,idim,ngrid,nleaf
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nH_eos,nCOM,damp_factor,cooling_switch,t_blast
  integer,dimension(1:nvector)::ind_leaf
  real(kind=8),dimension(1:nvector)::nH,T2,delta_T2,ekk,err,emag
  real(kind=8),dimension(1:nvector)::T2min,Zsolar,boost
  logical,dimension(1:nvector)::cooling_on=.true.
  real(dp),dimension(nions, 1:nvector):: xion
#ifdef RT
  integer::ig,iNp
  real(dp),dimension(1:ndim)::Fpnew
  real(dp),dimension(nrtgroups, 1:nvector):: Np, dNpdt=0d0
  real(dp),dimension(ndim, nrtgroups, 1:nvector):: Fp, dFpdt=0
  real(dp),dimension(ndim, 1:nvector):: p_gas
  real(dp)::scale_Np,scale_Fp,Npnew
#endif
#if NENER>0
  integer::irad
#endif

#ifdef HYDRO
  if(r%verbose.and.g%myid==1)write(*,'("   Entering cooling_fine for level ",I2)')ilevel

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#ifdef RT
  call rt_units(scale_Np,scale_Fp)
#endif
  ! Density for isothermal/polytropic EOS in H/cc
  nH_eos = r%eos_nH

  ! This is used only for cosmological runs
  if(r%cosmo)then
     nCOM = 200.0*g%omega_b*rhoc*(g%h0/100)**2/g%aexp**3*c%X/mH
     nH_eos = MAX(nCOM,nH_eos)
  endif

  ! Loop over cells
  do ind=1,twotondim
     ! Loop over octs with vector sweeps
     do igrid=m%head(ilevel),m%tail(ilevel),nvector

        ! Collect a vector of leaf cells
        ngrid=MIN(nvector,m%tail(ilevel)-igrid+1)
        nleaf=0
        do i=1,ngrid
           if(.NOT. m%grid(igrid+i-1)%refined(ind))then
              nleaf=nleaf+1
              ind_leaf(nleaf)=igrid+i-1
           endif
        end do
        if(nleaf.eq.0)cycle

        ! Compute rho
        do i=1,nleaf
           nH(i)=MAX(m%grid(ind_leaf(i))%uold(ind,1),r%smallr)
        end do

        ! Compute metallicity in solar units
        if(r%metal)then
           do i=1,nleaf
              Zsolar(i)=m%grid(ind_leaf(i))%uold(ind,r%imetal)/nH(i)/0.02d0
           end do
        else
           do i=1,nleaf
              Zsolar(i)=r%z_ave
           end do
        endif

        ! Total energy
        do i=1,nleaf
           T2(i)=m%grid(ind_leaf(i))%uold(ind,5)
        end do

        ! Kinetic energy
        do i=1,nleaf
           ekk(i)=0.0d0
        end do
        do idim=1,3
           do i=1,nleaf
              ekk(i)=ekk(i)+0.5d0*m%grid(ind_leaf(i))%uold(ind,idim+1)**2/nH(i)
           end do
        end do

        ! Non-thermal energies
        do i=1,nleaf
           err(i)=0.0d0
        end do
#if NENER>0
        do irad=1,nener
           do i=1,nleaf
              err(i)=err(i)+m%grid(ind_leaf(i))%uold(ind,5+irad)
           end do
        end do
#endif
        ! Magnetic energy
        do i=1,nleaf
           emag(i)=0.0d0
        end do
#ifdef MHD
        do idim=1,3
           do i=1,nleaf
              emag(i)=emag(i)+0.125d0*(m%grid(ind_leaf(i))%bold(ind,idim) &
                   &                  +m%grid(ind_leaf(i))%bold(ind,idim+3))**2
           end do
        end do
#endif
        ! Gas thermal pressure
        do i=1,nleaf
           T2(i)=(r%gamma-1.0d0)*(T2(i)-ekk(i)-err(i)-emag(i))
        end do

        ! Compute T2=T/mu in Kelvin
        do i=1,nleaf
           T2(i)=T2(i)/nH(i)*scale_T2
        end do

        ! Compute nH in H/cc
        do i=1,nleaf
           nH(i)=nH(i)*scale_nH
        end do

        ! Compute radiation self-shielding factor
        if(r%self_shielding)then
           do i=1,nleaf
              boost(i)=MAX(exp(-nH(i)/0.01d0),1.0D-20)
           end do
        else
           do i=1,nleaf
              boost(i)=1
           end do
        endif

        ! Compute ionization fraction
        if(r%neq_chem) then
           do ii=0,nIons-1
              do i=1,nleaf
                 xion(1+ii,i) = m%grid(ind_leaf(i))%uold(ind,r%iIons+ii)/m%grid(ind_leaf(i))%uold(ind,1)
              end do
           end do
        endif

        ! Get photon densities and flux magnitudes
#ifdef RT
        do ig=1,nrtgroups
           iNp=1+(ig-1)*(ndim+1)
           do i=1,nleaf
              Np(ig,i)          = scale_Np * m%grid(ind_leaf(i))%rtuold(ind,iNp)
              Fp(1:ndim, ig, i) = scale_Fp * m%grid(ind_leaf(i))%rtuold(ind,iNp+1:iNp+ndim)
           enddo
        end do
#endif

        ! Compute gas momentum for radiation force
#ifdef RT
        do i=1,nleaf
           p_gas(1:ndim,i) = m%grid(ind_leaf(i))%uold(ind,2:1+ndim) * scale_d * scale_v
        end do
#endif
        !==========================================
        ! Compute minimum temperature  for cooling
        ! floor or for isothermal/polytropic EOS.
        ! You can put your own polytrope EOS here.
        !==========================================
        if(r%eos_type==1)then ! Strictly isothermal
           do i=1,nleaf
              T2min(i) = r%eos_T2
           end do
        else if(r%eos_type==2)then ! Power law
           do i=1,nleaf
              T2min(i) = r%eos_T2*(nH(i)/nH_eos)**(r%eos_index-1.0d0)
           end do
        else if(r%eos_type==3)then ! Smooth transition from isothermal to power law
           do i=1,nleaf
              T2min(i) = r%eos_T2*(1+(nH(i)/nH_eos)**(r%eos_index-1.0d0))
           end do
        else if(r%eos_type==4)then ! Sharp transition from isothermal to power law
           do i=1,nleaf
              if(nH(i)<nH_eos)then
                 T2min(i) = r%eos_T2
              else
                 T2min(i) = r%eos_T2*(nH(i)/nH_eos)**(r%eos_index-1.0d0)
              endif
           end do
        endif

        ! Compute thermal temperature by subtracting the polytrope
        if(r%cooling.or.r%cooling_ism.or.r%neq_chem)then
           do i=1,nleaf
              T2(i) = min(max(T2(i)-T2min(i),T2_min_fix),r%T2max)
           end do
        endif

        ! Compute cooling time step in second
        dtcool = g%dtnew(ilevel)*scale_t

        ! Smooth RT update
#ifdef RT
        if(r%rt_smooth) then
           do ig=1,nrtgroups
              iNp=1+(ig-1)*(ndim+1)
              do i=1,nleaf ! Calc addition per sec to Np, Fp for current dt
                 Npnew = scale_Np * m%grid(ind_leaf(i))%rtunew(ind,iNp)
                 Fpnew = scale_Fp * m%grid(ind_leaf(i))%rtunew(ind,iNp+1:iNp+ndim)
                 dNpdt(ig,i)   = (Npnew - Np(ig,i)) / dtcool
                 dFpdt(:,ig,i) = (Fpnew - Fp(:,ig,i)) / dtcool
              end do
           end do
        end if
#endif
        !==========================================
        ! Compute net cooling at constant nH
        !==========================================
        if(r%cooling)then
           ! Use classical ramses cooling
           call solve_cooling(c,nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)
        else if(r%cooling_ism)then
           ! Use cooling from cooling_module_frig described in Audit & Hennebelle 2005
           call solve_cooling_ism(nH,T2,dtcool,delta_T2,r%gamma,1.4d0,nleaf)
        else if(r%neq_chem)then
           call rt_solve_cooling(r, tables, T2, xion, &
#ifdef RT
                & Np, Fp, p_gas, dNpdt, dFpdt, &
#endif
                & nH, cooling_on, Zsolar, dtcool, nleaf)
        endif

        ! Compute rho in code units
        do i=1,nleaf
           nH(i) = nH(i)/scale_nH
        end do

        ! Update fluid momentum and kinetic energy due to radiation force
#ifdef RT
        do i=1,nleaf
           m%grid(ind_leaf(i))%uold(ind,2:1+ndim) = p_gas(1:ndim,i) / scale_d / scale_v
        end do
        do i=1,nleaf
           ekk(i)=0.0d0
        end do
        do idim=1,3
           do i=1,nleaf
              ekk(i)=ekk(i)+0.5d0*m%grid(ind_leaf(i))%uold(ind,idim+1)**2/nH(i)
           end do
        end do
#endif
        ! Compute radiated energy in code units
        if(r%cooling.or.r%cooling_ism)then
           ! Compute net energy sink
           do i=1,nleaf
              delta_T2(i) = delta_T2(i)*nH(i)/scale_T2/(r%gamma-1.0d0)
           end do
        endif

        ! Compute fluid internal energy in code units
        if(r%cooling.or.r%cooling_ism.or.r%neq_chem)then
           do i=1,nleaf
              T2(i) = T2(i)*nH(i)/scale_T2/(r%gamma-1.0d0)
           end do
        endif

        ! Compute polytrope internal energy in code units
        do i=1,nleaf
           T2min(i) = T2min(i)*nH(i)/scale_T2/(r%gamma-1.0d0)
        end do

        ! Update fluid internal energy
        if(r%cooling.or.r%cooling_ism)then
           do i=1,nleaf
              T2(i) = T2(i) + delta_T2(i)
           end do
        endif

        ! Update ionization fraction
        if(r%neq_chem) then
           do ii=0,nions-1
              do i=1,nleaf
                 m%grid(ind_leaf(i))%uold(ind,r%iIons+ii) = xion(1+ii,i)*nH(i)
              end do
           end do
        endif

        ! Update entropy if dual energy scheme is activated
        if(r%entropy.and.r%dual_energy.GE.0)then
           if(r%isothermal)then ! use only polytrope energy
              do i=1,nleaf
                 m%grid(ind_leaf(i))%uold(ind,r%ientropy) = (T2min(i))/nH(i)**(r%gamma-1)*(r%gamma-1)
              end do
           else if(r%cooling.or.r%cooling_ism.or.r%neq_chem)then
              do i=1,nleaf
                 m%grid(ind_leaf(i))%uold(ind,r%ientropy) = (T2(i) + T2min(i))/nH(i)**(r%gamma-1)*(r%gamma-1)
              end do
           endif
        endif

        ! Update total fluid energy
        if(r%isothermal)then ! use only polytrope energy
           do i=1,nleaf
              m%grid(ind_leaf(i))%uold(ind,5) = T2min(i) + ekk(i) + err(i) + emag(i)
           end do
        else if(r%cooling.or.r%cooling_ism.or.r%neq_chem)then ! add polytrope to thermal energy
           do i=1,nleaf
              m%grid(ind_leaf(i))%uold(ind,5) = T2(i) + T2min(i) + ekk(i) + err(i) + emag(i)
           end do
        endif

        ! Update photon densities and flux magnitudes
#ifdef RT
        do ig=1,nrtgroups
           iNp=1+(ig-1)*(ndim+1)
           do i=1,nleaf
              m%grid(ind_leaf(i))%rtuold(ind,iNp) = max(Np(ig,i)/scale_Np,smallNp)
              m%grid(ind_leaf(i))%rtuold(ind,iNp+1:iNp+ndim) = Fp(1:ndim,ig,i)/scale_Fp
           enddo
        end do
#endif
     end do
     ! End loop over grid
  end do
  ! End loop over cells

  ! Compute new cooling table
  if(r%cooling.and.ilevel==r%levelmin.and.r%cosmo)then
     if(g%myid==1)write(*,*)'Computing new cooling table'
     call set_table(c,dble(g%aexp))
  endif
#endif

end subroutine cooling_fine
end module cooling_fine_module
