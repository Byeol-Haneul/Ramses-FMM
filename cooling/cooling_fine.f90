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
     call cooling_fine(pst%s%r,pst%s%g,pst%s%m,pst%s%cool,ilevel)
  endif

end subroutine r_cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cooling_fine(r,g,m,c,ilevel)
  use constants
  use amr_parameters, only:dp,ndim,nvector,twotondim
  use amr_commons, only:run_t,global_t,mesh_t
  use cooling_module, only:cooling_t,solve_cooling,T2_min_fix,set_table
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(cooling_t)::c
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for leaf cells at level ilevel
  !-------------------------------------------------------------------
  integer::i,ind,igrid,idim,ngrid,nleaf
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nH_eos,nCOM,damp_factor,cooling_switch,t_blast
  integer,dimension(1:nvector)::ind_leaf
  real(kind=8),dimension(1:nvector)::nH,T2,delta_T2,ekk,err,emag
  real(kind=8),dimension(1:nvector)::T2min,Zsolar,boost
#ifdef SOLVERmhd
  integer::neul=5
#else
  integer::neul=ndim+2
#endif
#if NENER>0
  integer::irad
#endif

#ifdef HYDRO
  if(r%verbose.and.g%myid==1)write(*,'("   Entering cooling_fine for level ",I2)')ilevel

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

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

        ! Compute gas pressure (thermal+polytrope)
        do i=1,nleaf
           T2(i)=m%grid(ind_leaf(i))%uold(ind,neul)
        end do
        do i=1,nleaf
           ekk(i)=0.0d0
        end do
        do idim=2,neul-1
           do i=1,nleaf
              ekk(i)=ekk(i)+0.5d0*m%grid(ind_leaf(i))%uold(ind,idim)**2/nH(i)
           end do
        end do
        do i=1,nleaf
           err(i)=0.0d0
        end do
#if NENER>0
        do irad=0,nener-1
           do i=1,nleaf
              err(i)=err(i)+m%grid(ind_leaf(i))%uold(ind,inener+irad)
           end do
        end do
#endif
        do i=1,nleaf
           emag(i)=0.0d0
        end do
#ifdef SOLVERmhd
        do idim=1,3
           do i=1,nleaf
              emag(i)=emag(i)+0.125d0*(m%grid(ind_leaf(i))%uold(ind,idim+5) &
                   &                  +m%grid(ind_leaf(i))%uold(ind,idim+nvar))**2
           end do
        end do
#endif
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

        ! Compute radiation damping factor
        if(r%self_shielding)then
           do i=1,nleaf
              boost(i)=MAX(exp(-nH(i)/0.01d0),1.0D-20)
           end do
        else
           do i=1,nleaf
              boost(i)=1
           end do
        endif

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
        if(r%cooling)then
           do i=1,nleaf
              T2(i) = min(max(T2(i)-T2min(i),T2_min_fix),r%T2max)
           end do
        endif

        ! Compute cooling time step in second
        dtcool = g%dtnew(ilevel)*scale_t

        ! Compute net cooling at constant nH
        if(r%cooling)then
           if(r%cooling_ism) then
              ! Use cooling from cooling_module_frig described in Audit & Hennebelle 2005
              call solve_cooling_ism(nH,T2,dtcool,delta_T2,r%gamma,1.4d0,nleaf)
           else
              ! Use classical ramses cooling
              call solve_cooling(c,nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)
           endif
        endif

        ! Compute rho
        do i=1,nleaf
           nH(i) = nH(i)/scale_nH
        end do

        ! Deal with cooling
        if(r%cooling)then
           ! Compute net energy sink
           do i=1,nleaf
              delta_T2(i) = delta_T2(i)*nH(i)/scale_T2/(r%gamma-1.0d0)
           end do
           ! Compute initial fluid internal energy
           do i=1,nleaf
              T2(i) = T2(i)*nH(i)/scale_T2/(r%gamma-1.0d0)
           end do
        endif

        ! Compute polytrope internal energy
        do i=1,nleaf
           T2min(i) = T2min(i)*nH(i)/scale_T2/(r%gamma-1.0d0)
        end do

        ! Update fluid internal energy
        if(r%cooling)then
           do i=1,nleaf
              T2(i) = T2(i) + delta_T2(i)
           end do
        endif

        ! Update entropy if dual energy scheme is activated
        if(r%entropy.and.r%dual_energy.GE.0)then
           if(r%isothermal)then ! use only polytrope energy
              do i=1,nleaf
                 m%grid(ind_leaf(i))%uold(ind,r%ientropy) = (T2min(i))/nH(i)**(r%gamma-1)*(r%gamma-1)
              end do
           else if(r%cooling)then
              do i=1,nleaf
                 m%grid(ind_leaf(i))%uold(ind,r%ientropy) = (T2(i) + T2min(i))/nH(i)**(r%gamma-1)*(r%gamma-1)
              end do
           endif
        endif

        ! Update total fluid energy
        if(r%isothermal)then ! use only polytrope energy
           do i=1,nleaf
              m%grid(ind_leaf(i))%uold(ind,neul) = T2min(i) + ekk(i) + err(i) + emag(i)
           end do
        else if(r%cooling)then ! add polytrope to thermal energy
           do i=1,nleaf
              m%grid(ind_leaf(i))%uold(ind,neul) = T2(i) + T2min(i) + ekk(i) + err(i) + emag(i)
           end do
        endif

     end do
     ! End loop over grid
  end do
  ! End loop over cells

  ! Compute new cooling table
  if(r%cooling.and.ilevel==r%levelmin.and.r%cosmo.and..not.r%cooling_ism)then
     if(g%myid==1)write(*,*)'Computing new cooling table'
     call set_table(c,dble(g%aexp))
  endif
#endif

end subroutine cooling_fine
end module cooling_fine_module
