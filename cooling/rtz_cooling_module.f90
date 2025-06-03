! Non-equlibrium (in H2, HI, HII, HeI, HeII, HeIII)
! cooling module for radiation-hydrodynamics.
! For details, see Rosdahl et al. (2013), Rosdahl & Teyssier (2015),
! and Nickerson, Teyssier, & Rosdahl (2018).
! Joki Rosdahl, Sarah Nickerson, Andreas Bleuler, and Romain Teyssier.
! NOTE: T2=T/mu, Np = photon density, Fp = photon flux,
module rtz_cooling_module
  use amr_parameters, only: ndim, dp, nvector
  use amr_commons, only: run_t
  use hydro_parameters, only: nion
  use rt_parameters
  use coolrates_module
  use constants
  use rtz_module
  implicit none

  private   ! default
  public rtz_updateRTGroups_CoolConstants, initialize_elements &
       ,rtz_solve_cooling, n_elements

  real(kind=8),parameter::T2_min_fix=1d-2 ! Min temperature [K]

  real(kind=8),parameter::T_min=0.1, T_frac=0.1
  real(kind=8),parameter::x_min=1d-20, x_fm=1d-6, x_frac=0.1
  real(kind=8),parameter::Np_min=1d-13, Np_frac=0.2
  real(kind=8),parameter::Fp_frac=0.5

  ! IR group index
  integer,parameter::iIR=1

  ! RTZ STUFF
  integer, parameter :: n_elements = 27
  type(Element) :: elements(n_elements)

CONTAINS

SUBROUTINE initialize_elements()
   ! Initializes the atomic data we need for the RTZ module

   implicit none
   integer::i
   
   ! Initialize everything to zero
   do i=1,n_elements
      elements(i)%atomic_number = -1
      elements(i)%n_ions = -1
      elements(i)%atomic_mass = -1.0
      elements(i)%z_solar = -1.0
      elements(i)%G0_photo_rate = 0.0
      elements(i)%n_mol = 0
      elements(i)%depletion = 1.0
   enddo

   ! Element 1: Hydrogen
   elements(1)%atomic_number = 1
   elements(1)%atomic_mass = 1.008
   elements(1)%z_solar = 1.0
   elements(1)%G0_photo_rate = 0.0 ! No subionizing PI
   elements(1)%n_ions = 2
   elements(1)%n_mol = 0
   elements(1)%depletion = 1.0

   ! Element 2: Helium
   elements(2)%atomic_number = 2
   elements(2)%atomic_mass = 4.0026
   elements(2)%z_solar = 8.51E-02
   elements(2)%G0_photo_rate = 0. ! No subionizing PI
   elements(2)%n_ions = 3
   elements(2)%depletion = 1.0

   ! Element 6: Carbon
   elements(6)%atomic_number = 6
   elements(6)%atomic_mass = 12.0107
   elements(6)%z_solar = 2.69E-04
   elements(6)%G0_photo_rate = 3.39E-10
   elements(6)%n_ions = elements(6)%atomic_number + 1
   elements(6)%depletion = 0.5
    
   ! Element 7: Nitrogen
   elements(7)%atomic_number = 7
   elements(7)%atomic_mass = 14.0067
   elements(7)%z_solar = 6.76E-05
   elements(7)%G0_photo_rate = 0.0 ! No subionizing PI
   elements(7)%n_ions = elements(7)%atomic_number + 1
   elements(7)%depletion = 0.6

   ! Element 8: Oxygen
   elements(8)%atomic_number = 8
   elements(8)%atomic_mass = 15.9994
   elements(8)%z_solar = 4.90E-04
   elements(8)%G0_photo_rate = 0.0 ! No subionizing PI
   elements(8)%n_ions = elements(8)%atomic_number + 1
   elements(8)%depletion = 0.73

   ! Element 10: Neon
   elements(10)%atomic_number = 10
   elements(10)%atomic_mass = 20.1797
   elements(10)%z_solar = 8.51E-05
   elements(10)%G0_photo_rate = 0.0 ! No subionizing PI
   elements(10)%n_ions = elements(10)%atomic_number + 1
   elements(10)%depletion = 1.0

   ! Element 12: Magnesium
   elements(12)%atomic_number = 12
   elements(12)%atomic_mass = 24.305
   elements(12)%z_solar = 3.98E-05
   elements(12)%G0_photo_rate = 6.59E-11 
   elements(12)%n_ions = elements(12)%atomic_number + 1
   elements(12)%depletion = 0.16

   ! Element 14: Silicon
   elements(14)%atomic_number = 14
   elements(14)%atomic_mass = 28.0855
   elements(14)%z_solar = 3.24E-05
   elements(14)%G0_photo_rate = 4.47E-09
   elements(14)%n_ions = elements(14)%atomic_number + 1
   elements(14)%depletion = 0.1

   ! Element 16: Sulfur
   elements(16)%atomic_number = 16
   elements(16)%atomic_mass = 32.065
   elements(16)%z_solar = 1.32E-05
   elements(16)%G0_photo_rate = 1.13E-09
   elements(16)%n_ions = elements(16)%atomic_number + 1
   elements(16)%depletion = 1.0

   ! Element 26: Iron
   elements(26)%atomic_number = 26
   elements(26)%atomic_mass = 55.854
   elements(26)%z_solar = 3.16E-05
   elements(26)%G0_photo_rate = 4.71E-10
   elements(26)%n_ions = elements(26)%atomic_number + 1
   elements(26)%depletion = 0.01

END SUBROUTINE initialize_elements

!!$!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!!$SUBROUTINE update_UVrates(aexp)
!!$! Set the UV ionization and heating rates according to the given aexp.
!!$!-------------------------------------------------------------------------
!!$  use UV_module
!!$  use amr_parameters,only:haardt_madau
!!$  real(kind=8)::aexp
!!$!------------------------------------------------------------------------
!!$  UVrates=0.
!!$  if(.not. haardt_madau) RETURN
!!$
!!$  call inp_UV_rates_table(1./aexp - 1., UVrates, .true.)
!!$
!!$END SUBROUTINE update_UVrates

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE rtz_solve_cooling(r, tables, T2, xion, nElement, &
#ifdef RT
     & Np, Fp, p_gas, dNpdt, dFpdt, ilevel, &
#endif
     & dt, nCell)

#ifdef RT
  use neq_cooling_module, only: reduce_flux
#endif
  ! Semi-implicitly solve for new temperature, ionization states,
  ! photon density/flux, and gas velocity in a number of cells.
  ! Parameters:
  ! tables  => object containing all non equilibrium chemistry tables
  ! T2     <=> T/mu [K]
  ! xion   <=> 27x27 array with all of the mass fractions of each ionization state
  ! Np     <=> nrtgrp photon number densities [cm-3]
  ! Fp     <=> nrtgrp * ndim photon number fluxes [cm-2 s-1]
  ! p_gas  <=> ndim gas momentum densities [cm s-1 g cm-3]
  ! dNpdt   =>  Op split increment in photon densities during dt
  ! dFpdt   =>  Op split increment in photon flux magnitudes during dt
  ! nElement=>  Number density of each element [cm^-3] (length 27)
  ! c_switch=>  Cooling switch (1 for cool/heat, 0 for no cool/heat) (OFF)
  ! dt      =>  Timestep size             [s]
  ! nCell   =>  Number of cells (length of all the above vectors)
  ! ilevel  =>  Refinement levels
  !
  ! We use a slightly modified method of Anninos et al. (1997).
  !-------------------------------------------------------------------------
  implicit none
  type(run_t):: r
  type(neq_cooling_t):: tables
  real(kind=8),dimension(1:nvector):: T2
  real(kind=8),dimension(1:n_elements, 1:n_elements, 1:nvector):: xion
  real(kind=8),dimension(1:n_elements, 1:nvector):: nElement
  real(kind=8),dimension(1:nvector):: nH
#ifdef RT
  real(kind=8),dimension(1:ndim, 1:nvector):: p_gas
  real(kind=8),dimension(1:nrtgrp, 1:nvector):: Np, dNpdt
  real(kind=8),dimension(1:ndim, 1:nrtgrp, 1:nvector):: Fp, dFpdt
  integer::ilevel
#endif
!  logical,dimension(1:nvector):: c_switch
  real(kind=8)::dt
  integer::ncell
  !--------------------------------------------------------
  real(kind=8),dimension(1:nvector):: tLeft, ddt
  logical:: dt_ok
  real(kind=8):: dt_rec
  real(kind=8):: dT2
  real(kind=8),dimension(1:n_elements,1:n_elements):: dXion ! TODO(code): fix this
  integer::i, ia, nAct, nAct_next, loopcnt, code
  integer,dimension(1:nvector):: indAct              ! Active cell indexes
  real(kind=8):: one_over_x_FRAC, one_over_T_FRAC
#ifdef RT
  integer::ig
  real(kind=8):: one_over_rt_c_cgs, one_over_egy_IR_erg
  real(kind=8):: one_over_Np_FRAC, one_over_Fp_FRAC
  real(kind=8),dimension(1:ndim):: dp_gas
  real(kind=8),dimension(nrtgrp):: dNp
  real(kind=8),dimension(1:ndim, 1:nrtgrp):: dFp
  real(kind=8),dimension(1:nrtgrp):: group_egy_ratio, group_egy_erg
#endif
  integer*8,dimension(20)::loopCodes=0
  integer::iElement, ion_fracs, idx_max
  real(kind=8)::current_mass_frac
  integer:: i_interp, convergence_counter
  integer :: base_unit = 100
  integer :: element_unit, j
  character(len=50) :: element_filename

  ! Store some temporary variables reduce computations
  one_over_T_FRAC = 1d0 / T_FRAC
  one_over_x_FRAC = 1d0 / x_FRAC
#ifdef RT
  one_over_Np_FRAC = 1d0 / Np_FRAC
  one_over_Fp_FRAC = 1d0 / Fp_FRAC
  one_over_rt_c_cgs = 1d0 / tables%rt_c_cgs(ilevel)
  group_egy_erg(1:nrtgrp) = r%group_egy(1:nrtgrp) * eV2erg
  if(r%rt_isIR) then
     group_egy_ratio(1:nrtgrp) = r%group_egy(1:nrtgrp) / r%group_egy(iIR)
     one_over_egy_IR_erg = 1d0 / group_egy_erg(iIR)
  endif
#endif
  !-----------------------------------------------------------------------

  ! Check if we are running in equilibrium mode
  if (r%rtz_equilibrium_test) then

      ! Open files for all elements
      do i = 1, n_elements
         if (elements(i)%atomic_number .gt. 0) then
            element_unit = base_unit + i
            write(element_filename, '("element_", I0, "_ions.dat")') i
            open(unit=element_unit, file=element_filename, status='unknown')
         end if
      end do

      do i_interp = 1,300
         ! Initialize the convergence counter
         convergence_counter = 0

         ! Set the temperature
         r%neq_TConst = 10.d0**(((8.d0 - 2.d0) * (real(i_interp,kind=8) - 1.d0)/(300.d0-1.d0)) + 2.d0)

         tleft(1:ncell) = 1.d20             ! Set to an arbitrarily large number
         ddt(1:ncell) = 10000.d0 * 365.25d0 * 60.d0 * 60.d0 ! First guess at sub-timestep lengths

         do i=1,ncell
            indact(i) = i                   !      Set up indexes of active cells

            ! set nH -- TODO(code): maybe put bounds on nH here
            nH(i) = nElement(1,i)

            ! Ensure all state vars are legal:
            T2(i) = MAX(T2(i), T2_min_fix)

            ! Make sure the ionization fractions don't go below the min or max
            xion(1:n_elements,1:n_elements,i) = MIN(MAX(xion(1:n_elements,1:n_elements,i), x_MIN),1d0)

            ! Loop over each element and ensure that the ionization fractions
            ! sum to 1
            do iElement=1,n_elements
               ! Check if we are actually using that element
               if (elements(iElement)%atomic_number .gt. 0) then
                  ! Get the indices of the number of ions + number of molecules
                  ion_fracs = elements(iElement)%n_ions + elements(iElement)%n_mol

                  ! Sum the ionization mass fractions
                  current_mass_frac = sum(xion(iElement,1:ion_fracs,i))

                  ! Get the index of the maximum 
                  idx_max = MAXLOC(xion(iElement,1:ion_fracs,i), DIM=1)

                  ! Update the max ionization state so that it sums to 1
                  xion(iElement,idx_max,i) = 1.d0 - (current_mass_frac - xion(iElement,idx_max,i))
               end if
            end do
         end do

         ! Loop until all cells have tleft=0
         ! **********************************************
         nAct=nCell                                      ! Currently active cells
         loopcnt=0 !; n_cool_cells=n_cool_cells+nCell     !             Statistics
         do while (nAct .gt. 0)      ! Iterate while there are still active cells
            loopcnt=loopcnt+1 !  ;   tot_cool_loopcnt=tot_cool_loopcnt+nAct
            nAct_next=0                     ! Active cells for the next iteration
            do ia=1,nAct                             ! Loop over the active cells
               i = indAct(ia)                        !                 Cell index
               call rtz_cool_step(i)
               if(.not. dt_ok) then
                  convergence_counter = 0 ! Reset the convergence counter
                  ddt(i)=ddt(i)/2.                    ! Try again with smaller dt
                  nAct_next=nAct_next+1 ; indAct(nAct_next) = i
                  loopCodes(code) = loopCodes(code)+1
                  cycle
               endif
               convergence_counter = convergence_counter + 1
               ! Update the cell state (advance the time by ddt):
               T2(i) = T2(i) + dT2
               ! Check for convergence
               xion(:,:,i) = xion(:,:,i) + dXion(:,:)
               if (convergence_counter .gt. 1000) then
                  tleft(i) = 0.0 ! Finish the cell if we have reached convergence
               else
                  tleft(i)=tleft(i)-ddt(i)
                  ! Take at least 100 iterations
                  if (convergence_counter .lt. 1000) then
                     tleft(i)=max(tleft(i),ddt(i))
                  endif
               endif
               if(tleft(i) .gt. 0.) then           ! Not finished with this cell
                  nAct_next=nAct_next+1 ; indAct(nAct_next) = i
               else if(tleft(i) .lt. 0.) then        ! Overshot by abs(tleft(i))
                  print*,'In rtz_solve_cooling: tleft < 0  !!'
                  stop
               endif
               ddt(i)=min(dt_rec,tleft(i))    ! Use recommended dt from rtz_cool_step
            end do ! end loop over active cells
            nAct=nAct_next
         end do ! end iterative loop

         ! Write hydrogen for debugging
         write(*,*) r%neq_TConst, loopcnt, xion(1,1,1), xion(1,2,1)

         ! Write data to file
         do i = 1, n_elements
            if (elements(i)%atomic_number .gt. 0) then
               element_unit = base_unit + i
               write(element_unit,'(ES15.6, I8, *(ES15.6))') r%neq_TConst, loopcnt, &
                  (xion(i,j,1), j=1,elements(i)%n_ions)
            end if
         end do

      end do

      ! Close all element files
      do i = 1, n_elements
         if (elements(i)%atomic_number .gt. 0) then
            close(base_unit + i)
         end if
      end do

      write(*,*) '!************************************************!'
      stop "Program terminated due to equilibrium test"

  ! Otherwise perform the normal loop
  else
      tleft(1:ncell) = dt                !       Time left in dt for each cell
      ddt(1:ncell) = dt                  ! First guess at sub-timestep lengths

      do i=1,ncell
         indact(i) = i                   !      Set up indexes of active cells

         ! set nH -- TODO(code): maybe put bounds on nH here
         nH(i) = nElement(1,i)

         ! Ensure all state vars are legal:
         T2(i) = MAX(T2(i), T2_min_fix)

         ! Make sure the ionization fractions don't go below the min or max
         xion(1:n_elements,1:n_elements,i) = MIN(MAX(xion(1:n_elements,1:n_elements,i), x_MIN),1d0)

         ! Loop over each element and ensure that the ionization fractions
         ! sum to 1
         do iElement=1,n_elements
            ! Check if we are actually using that element
            if (elements(iElement)%atomic_number .gt. 0) then
               ! Get the indices of the number of ions + number of molecules
               ion_fracs = elements(iElement)%n_ions + elements(iElement)%n_mol

               ! Sum the ionization mass fractions
               current_mass_frac = sum(xion(iElement,1:ion_fracs,i))

               ! Get the index of the maximum 
               idx_max = MAXLOC(xion(iElement,1:ion_fracs,i), DIM=1)

               ! Update the max ioniztion state so that it sums to 1
               xion(iElement,idx_max,i) = 1.d0 - (current_mass_frac - xion(iElement,idx_max,i))
            end if
         end do
#ifdef RT
         do ig=1,nrtgrp
            Np(ig,i) = MAX(smallNp, Np(ig,i))
            call reduce_flux(Fp(:,ig,i),Np(ig,i)*tables%rt_c_cgs(ilevel))
         end do
#endif
      end do

      ! Loop until all cells have tleft=0
      ! **********************************************
      nAct=nCell                                      ! Currently active cells
      loopcnt=0 !; n_cool_cells=n_cool_cells+nCell     !             Statistics
      do while (nAct .gt. 0)      ! Iterate while there are still active cells
         loopcnt=loopcnt+1 !  ;   tot_cool_loopcnt=tot_cool_loopcnt+nAct
         nAct_next=0                     ! Active cells for the next iteration
         do ia=1,nAct                             ! Loop over the active cells
            i = indAct(ia)                        !                 Cell index
            call rtz_cool_step(i)
      !         if(loopcnt .gt. 100000) then
      !            call display_coolinfo(.true., loopcnt, i, dt-tleft(i), dt, ddt(i), nH(i), &
      ! #ifdef RT
      !                 &                Np(:,i), Fp(:,:,i), p_gas(:,i), dNp, dFp, dp_gas, ilevel, &
      ! #endif
      !                 &                T2(i), xion(:,:,i), dT2, dXion, code)
      !         endif
            if(.not. dt_ok) then
               ddt(i)=ddt(i)/2.                    ! Try again with smaller dt
               nAct_next=nAct_next+1 ; indAct(nAct_next) = i
               loopCodes(code) = loopCodes(code)+1
               cycle
            endif
            ! Update the cell state (advance the time by ddt):
            T2(i) = T2(i) + dT2
            xion(:,:,i) = xion(:,:,i) + dXion(:,:)
#ifdef RT
            Np(:,i) = Np(:,i) + dNp(:)
            Fp(:,:,i) = Fp(:,:,i) + dFp(:,:)
            p_gas(:,i) = p_gas(:,i) + dp_gas(:)
#endif
            tleft(i)=tleft(i)-ddt(i)
            if(tleft(i) .gt. 0.) then           ! Not finished with this cell
               nAct_next=nAct_next+1 ; indAct(nAct_next) = i
            else if(tleft(i) .lt. 0.) then        ! Overshot by abs(tleft(i))
               print*,'In rtz_solve_cooling: tleft < 0  !!'
               stop
            endif
            ddt(i)=min(dt_rec,tleft(i))    ! Use recommended dt from rtz_cool_step
         end do ! end loop over active cells
         nAct=nAct_next
      end do ! end iterative loop

  end if

contains

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  SUBROUTINE rtz_cool_step(icell)
    ! Compute change in cell state in timestep ddt(icell), or set in dt_rec
    ! a recommendation for new timestep if ddt(icell) proves too large.
    ! T2      => T/mu [K]                               -- dT2 is new value
    ! xion    => nion ionization fractions              --     dXion is new
    ! Np      => nrtgrp photon number densities [cm-3]  -- dNp is new value
    ! Fp      => nrtgrp * ndim photon fluxes [cm-2 s-1] -- dFp is new value
    ! p_gas   => ndim gas momenta [cm s-1 g cm-3]       --    dp_gas is new
    ! dNpdt   =>  Op split increment in photon densities during dt
    ! dFpdt   =>  Op split increment in photon flux magnitudes during dt
    ! nH      =>  Hydrogen number densities [cm-3]
    ! c_switch=>  Cooling switch (1 for cool/heat, 0 for no cool/heat) (OFF)
    ! dt      =>  Timestep size [s]
    ! dt_ok   <=  .f. if timestep constraints were broken, .t. otherwise
    ! dt_rec  <=  Recommended timesteps for next iteration
    ! code    <= Error code in cool step, if dt_ok=.f.
    !
    ! The original values, T2, xion etc, must stay unchanged, while dT2,
    ! dxion etc contain the new values (the difference at the end of the
    ! routine).
    !-----------------------------------------------------------------------
    use amr_commons
    use const
    use neq_cooling_module, only: getMu
    use collisional_ionization_module
    use recombination_module
    use charge_exchange_module
    use dust_recombination_module
    use photoionization_UVB_module
    use cosmic_ray_ionization_module
    implicit none
    integer, intent(in):: icell
    !-----------------------------------------------------------------------
    real(kind=8),dimension(nion):: alpha, beta, nN, nI
    real(kind=8):: dUU, fracMax, x_tot
    real(kind=8):: mu, TK, nHe, ne, neInit, Hrate
    real(kind=8):: xHI,dxHI, xH2=0d0,dXH2=0d0, xHeI,dxHeI
    real(kind=8):: Crate, dCdT2, X_nHkb, rate, dRate, cr, de=0d0
    real(kind=8):: photoRate, metal_tot, metal_prime, ss_factor, f_dust
#ifdef RT
    integer::igroup,idim
    real(kind=8),dimension(ndim):: dmom
    real(kind=8),dimension(nrtgrp):: recRad, phAbs, phSc, dustAbs
    real(kind=8),dimension(nrtgrp):: dustSc, kAbs_loc, kSc_loc
    real(kind=8),dimension(nrtgrp,nion)::signc
    real(kind=8):: rt_c_fraction, rt_c_cgs
    real(kind=8):: TR, one_over_C_v, E_rad, dE_T
    real(kind=8):: G0, eff_peh, cdex
    real(kind=8):: fluxMag, mom_fact
#endif
    real(kind=8):: rho, ncr
    logical:: newAtomicCons=.true.
    !-----------------------------------------------------------------------
    ! Variables specific to RTZ
    real(kind=8):: alpha_H2_loc, de_H2, xe
    real(kind=8):: dust_effective_number_density, dust_to_gas_mass_ratio_over_mw
    real(kind=8):: HI_number_density, HII_number_density
    real(kind=8):: paired_ion_number_density
    real(kind=8):: total_cosmic_ray_ionization_rate, H2_cosmic_ray_ionization_rate
    real(kind=8):: phi_s, cosmic_ray_scale_factor, primary_cosmic_ray_ionization_rate
    real(kind=8):: UV_background_G0
    integer:: atomic_number, n_ions, i_other_Element, i_other_Ion, iIon
    real(kind=8):: Zsolar
    !-----------------------------------------------------------------------

    ! RTZ variable initialization
    dust_to_gas_mass_ratio_over_mw = 0.d0
    total_cosmic_ray_ionization_rate = 1.d-16

    ! Scale solar metallicity based on oxygen abundance
    !Zsolar = 10.d0**((12.d0 + log10((nElement(8,icell)+1.d-20)/nElement(1,icell))) - 8.69d0)
    Zsolar = 1.d-40
    ! END RTZ variable initialization


#ifdef RT
    signc=tables%signc(:,:,ilevel)
    rt_c_fraction = r%rt_c_fraction(ilevel)
    rt_c_cgs = tables%rt_c_cgs(ilevel)
#endif
    dt_ok=.false.
    nHe=0.25*nH(icell)*r%Y_He/r%X_H       ! Helium number density
    ! U contains the original values, dU the updated ones
    dT2 = T2(icell) ; dXion(:,:) = xion(:,:,icell)
#ifdef RT
    dNp(:) = Np(:,icell) ; dFp(:,:) = Fp(:,:,icell)
    dp_gas(:) = p_gas(:,icell)
#endif

    mu = getMu(r, dXion, dT2)
    TK = dT2 * mu                                        !      Temperature
    if(r%neq_isTconst) TK=r%neq_Tconst                   ! Force constant T
    ne = getNe(dXion, nElement(:,icell))
    neInit = ne
    fracMax = 0d0 ! Max fractional update, to check if dt can be increased
    ss_factor = 1d0                  ! UV background self_shielding factor
    if(r%self_shielding) ss_factor = exp(-nH(icell)/1d-2)
    rho = nH(icell) / r%X_H * mH ! TODO(code): update this to the correct value

    ! RTZ -- initialize cosmiv ray variables
    ! Measure phi_s for secondary CR ionization
    phi_s = secondary_cr_rates(xe)
    ! Now calculate the total CR ionization rate
    primary_cosmic_ray_ionization_rate = total_cosmic_ray_ionization_rate / (1.d0 + phi_s)
    H2_cosmic_ray_ionization_rate = 2.d0 * primary_cosmic_ray_ionization_rate * (1.d0 + phi_s)
    ! END RTZ

#ifdef RT
    ! Set dust opacities--------------------------------------------------
    kAbs_loc = r%kappaAbs
    kSc_loc = r%kappaSc
    if(r%is_kIR_T) then                            ! k_IR depends on T_rad
          ! For the radiation temperature,  weigh the energy in each group
            ! by its opacity over IR opacity (derived from IR temperature)
       E_rad = group_egy_erg(iIR) * dNp(iIR)
       TR = max(0d0,(E_rad*rt_c_fraction/a_r)**0.25)      ! IR temperature
       kAbs_loc(iIR) = r%kappaAbs(iIR) * (TR/10d0)**2
       do iGroup=1,nrtgrp
          if(iGroup .ne. iIR)                                            &
               E_rad = E_rad + kAbs_loc(iGroup) / kAbs_loc(iIR)          &
               * group_egy_erg(iGroup) * dNp(iGroup)
       end do
       TR = max(0d0,(E_rad*rt_c_fraction/a_r)**0.25)  ! Rad. temperature
       if(r%rt_T_rad) then      ! Use radiation temperature for everything
          dT2 = TR/mu ;   TK = TR
       endif
                 ! Set the IR opacities according to the rad. temperature:
       kAbs_loc(iIR) = r%kappaAbs(iIR) * (TR/10d0)**2 * exp(-TR/1d3)
       kSc_loc(iIR)  = r%kappaSc(iIR)  * (TR/10d0)**2 * exp(-TR/1d3)
    endif ! if(is_kIR_T)
                         ! Set dust absorption and scattering rates [s-1]:
    ! TODO(code): update this to harley's dust model
    dustAbs(:) =kAbs_loc(:) *rho*Zsolar*f_dust*rt_c_cgs
    dustSc(iIR)=kSc_loc(iIR)*rho*Zsolar*f_dust*rt_c_cgs

    ! UPDATE PHOTON DENSITY AND FLUX *************************************
    if(r%rt_advect) then
       recRad(1:nrtgrp)=0. ; phAbs(1:nrtgrp)=0.
       ! Scattering rate; reduce the photon flux, but not photon density:
       phSc(1:nrtgrp)=0.

       ! HKnote: we do not allow OTSA with RTZ (for now)

       ! ABSORPTION/SCATTERING OF PHOTONS BY GAS
       do igroup=1,nrtgrp       ! ----------------Ionization absorbtion
          !TODO(code): need to loop over all elements and ionz here
          phAbs(igroup) = SUM(nN(:)*signc(igroup,:)*r%ssh2(igroup)) ! s-1
       end do
       ! IR, optical and UV depletion by dust absorption: ----------------
       ! IR scattering/abs on dust (abs after T update)
       if(r%rt_isIR) phSc(iIR)  = phSc(iIR) + dustSc(iIR)
       do igroup=1,nrtgrp      ! Deplete photons, since they go into IR
          if( .not. (r%rt_isIR .and. igroup.eq.iIR) ) & ! IR done elsewhere
               phAbs(igroup) = phAbs(igroup) + dustAbs(igroup)
       end do

       dmom(1:ndim)=0d0
       do igroup=1,nrtgrp     ! ----------------- Do the update of N and F
          dNp(igroup)= MAX(smallNp,                                      &
                        (ddt(icell)*(recRad(igroup)+dNpdt(igroup,icell)) &
                                    +dNp(igroup))                        &
                        / (1d0+ddt(icell)*phAbs(igroup)))

          dUU = ABS(dNp(igroup)-Np(igroup,icell))                        &
                /(Np(igroup,icell)+Np_MIN) * one_over_Np_FRAC
          if(dUU .gt. 1d0) then
             code=1 ;   RETURN                        ! ddt(icell) too big
          endif
          fracMax=MAX(fracMax,dUU)      ! To check if ddt can be increased

          do idim=1,ndim
             dFp(idim,igroup) = &
                  (ddt(icell)*dFpdt(idim,igroup,icell)+dFp(idim,igroup)) &
                  /(1d0+ddt(icell)*(phAbs(igroup)+phSc(igroup)))
          end do
          call reduce_flux(dFp(:,igroup),dNp(igroup)*rt_c_cgs)

          do idim=1,ndim
             dUU = ABS(dFp(idim,igroup)-Fp(idim,igroup,icell))           &
                  / (ABS(Fp(idim,igroup,icell))+Np_MIN*rt_c_cgs)         &
                  * one_over_Fp_FRAC
             if(dUU .gt. 1d0) then
                code=2 ;   RETURN                     ! ddt(icell) too big
             endif
             fracMax=MAX(fracMax,dUU)   ! To check if ddt can be increased
          end do

       end do

       do igroup=1,nrtgrp ! -------Momentum transfer from photons to gas:
          mom_fact = ddt(icell) * (phAbs(igroup) + phSc(igroup)) &
               * group_egy_erg(igroup) * one_over_clight

          if(r%rt_isoPress .and. .not. (r%rt_isIR .and. igroup==iIR)) then
             ! rt_isoPress: assume f=1, where f is reduced flux.
             fluxMag=sqrt(sum((dFp(:,igroup))**2))
             if(fluxMag .gt. 0d0) then
                mom_fact = mom_fact * dNp(igroup) / fluxMag
             else
                mom_fact = 0d0
             endif
          else
             mom_fact = mom_fact * one_over_rt_c_cgs
          end if

          do idim = 1, ndim
             dmom(idim) = dmom(idim) + dFp(idim,igroup) * mom_fact
          end do
       end do
       dp_gas = dp_gas + dmom * r%rt_pressBoost      ! update gas momentum

       ! Add absorbed UV/optical energy to IR:----------------------------
       if(r%rt_isIR) then
          do igroup=iIR+1,nrtgrp
             dNp(iIR) = dNp(iIR) + dustAbs(igroup) * ddt(icell)          &
                  * dNp(igroup) * group_egy_ratio(igroup)
          end do
       endif
       ! -----------------------------------------------------------------
    endif !if(rt)
#endif

    ! UPDATE TEMPERATURE *************************************************
    !if(c_switch(icell) .and. .not. rt_isTconst .and. .not. r%rt_T_rad) then
    if(.not. r%neq_isTconst .and. .not. r%rt_T_rad) then
       Hrate = 0.d0          !  Heating rate [erg cm-3 s-1]
       Crate = 0.d0          ! Cooling
       dCdT2 = dCdT2 * mu                            ! dC/dT2 = mu * dC/dT
       metal_tot=0d0 ; metal_prime=0d0                     ! Metal cooling
       !TODO(code) update X_nHkb to the correct value
      !  X_nHkb = r%X_H/(1.5 * nH(icell) * kB)        ! Multiplication factor
      !  rate  = X_nHkb*(Hrate - Crate - Zsolar*metal_tot)
      !  dRate = -X_nHkb*(dCdT2 + Zsolar*metal_prime)     ! dRate/dT2
      !                                                ! 1st order dt constr
      !  dUU   = ABS(MAX(T2_min_fix, T2(icell)+rate*ddt(icell))-T2(icell))
      !                                                       ! New T2 value
      !  dT2   = MAX(T2_min_fix &
      !             ,T2(icell)+rate*ddt(icell)/(1.-dRate*ddt(icell)))
      !  dUU   = MAX(dUU, ABS(dT2-T2(icell))) / (T2(icell)+T_MIN) &
      !                   *one_over_T_FRAC
      dT2 = 0.d0
      dUU = 0.d0
       if(dUU .gt. 1.) then                                     ! 10% rule
          code=3 ; RETURN
       endif
       fracMax=MAX(fracMax,dUU)
       TK=dT2*mu
    endif

#ifdef RT
    if(r%rt_isIR) then
       if(kAbs_loc(iIR) .gt. 0d0 .and. .not. r%rt_T_rad) then
          ! Evolve IR-Dust equilibrium temperature------------------------
          ! Delta (Cv T)= ( c_red/lambda E - c/lambda a T^4)
          !           / ( 1/Delta t + 4 c/lambda/C_v a T^3 + c_red/lambda)
          one_over_C_v = mH*mu*(r%gamma-1d0) / (rho*kB)
          E_rad = group_egy_erg(iIR) * dNp(iIR)
          dE_T = (rt_c_cgs * E_rad - c_cgs*a_r*TK**4)                    &
               /(1d0/(kAbs_loc(iIR) * Zsolar * rho * ddt(icell))  &
               +4d0*c_cgs * one_over_C_v *a_r*TK**3+rt_c_cgs)
          dT2 = dT2 + 1d0/mu * one_over_C_v * dE_T
          dNp(iIR) = dNp(iIR) - dE_T * one_over_egy_IR_erg

          dT2 = max(T2_min_fix,dT2)
          dNp(iIR) = max(dNp(iIR), smallNp)
          ! 10% rule for photon density:
          dUU = ABS(dNp(iIR)-Np(iIR,icell)) / (Np(iIR,icell)+Np_MIN)     &
                                            * one_over_Np_FRAC
          if(dUU .gt. 1.) then
             code=4 ;   RETURN
          endif
          fracMax=MAX(fracMax,dUU)

          dUU   = ABS(dT2-T2(icell)) / (T2(icell)+T_MIN) * one_over_T_FRAC
          if(dUU .gt. 1.) then                           ! 10% rule for T2
             code=5 ; RETURN
          endif
          fracMax=MAX(fracMax,dUU)
          TK=dT2*mu
          call reduce_flux(dFp(:,iIR),dNp(iIR)*rt_c_cgs)
       endif
    endif
#endif

    !/////////////////////////////////////////
    !//           UPDATE MOLECULES          //
    !/////////////////////////////////////////
    alpha_H2_loc = 0.d0

    !/////////////////////////////////////////
    !//       UPDATE IONIZATION STATES      //
    !/////////////////////////////////////////

    ! Initialize the change to zero
    dUU = 0.d0

    ! Get the effective dust number density
    dust_effective_number_density = nElement(1, icell) * dust_to_gas_mass_ratio_over_mw

    ! Loop over all elements
    do iElement = 1,n_elements
       if (elements(iElement)%atomic_number > 0) then
          ! Get the atomic number
          atomic_number = elements(iElement)%atomic_number

          ! Get the number of ions
          n_ions = elements(iElement)%n_ions

         !  write(*,*) iElement, atomic_number, n_ions, icell, nElement(iElement, icell)

          ! Loop over the number of ions
          do iIon = 1,n_ions

             !/////////////////////////
             !//       Creation      //
             !/////////////////////////
             cr = 0.d0

             ! Account for molecular hydrogen
            !  if (iElement.eq.1 .and. iIon.eq.0 .and. r%isH2) then
            !     ! Note: no factor of 2 needed since is 2*xH2
            !     cr = cr + de_H2 * dXion(1,2)
            !  end if

             ! Recombination of the more excited ionization state
             if (iIon.lt.n_ions) then  
                cr = cr + (recombination(TK, iIon+1, iElement) * ne * dXion(iElement,iIon+1))
             end if

             ! Collisional ionization of the less excited state
             if (iIon.gt.1) then 
                cr = cr + (collisional_ionization(TK, iIon-1, iElement) * ne * dXion(iElement,iIon-1))
             end if

             ! Photoionization of the less excited state
            !  if (iIon.gt.1) then 
            !     cr = cr + (HM12_UVB_z(iElement,iIon-1,1) * ss_factor * dXion(iElement,iIon-1))
            !  end if

             ! Photoionization by sub-ionizing ISRF --> only impacts lowest ionization states
            !  if (iIon.eq.2) then 
            !     cr = cr + (UV_background_G0 * elements(iElement)%G0_photo_rate * dXion(iElement,iIon-1))
            !  end if

             ! Cosmic ray ionization of the less excited state
            !  if (iIon > 1) then 
            !     cr = cr + (cosmic_ray_ionization_rates(iElement,iIon-1) * total_cosmic_ray_ionization_rate * dXion(iElement,iIon-1))
            !  end if

             ! Cosmic ray ionization of the less excited state from induced UV
            !  if (iIon.eq.2) then 
            !     cr = cr + (cosmic_ray_ionization_rates_induced_UV(iElement) * cosmic_ray_scale_factor * dXion(iElement,iIon-1))
            !  end if

             ! Recombination on dust from the more excited state
            !  if (iIon.lt.n_ions) then 
            !     cr = cr + (dust_recombination(iIon+1, i, TK, UV_background_G0, ne) * dust_effective_number_density * dXion(iElement,iIon+1))
            !  end if

             !TODO(code): add creation from local radiation

             !/////////////////////////
             !//     Destruction     //
             !/////////////////////////
             de = 0.d0

             ! Account for molecular hydrogen
            !  if (iElement.eq.1 .and. iIon.eq.1 .and. r%isH2) then
            !    de = de + alpha_H2_loc
            !  end if

             ! Collisional ionization 
             if (iIon .lt. n_ions) then 
                de = de + (collisional_ionization(TK, iIon, iElement) * ne)
             end if

             ! Photoionization
            !  if (iIon .lt. n_ions) then 
            !     de = de + (HM12_UVB_z(iElement,iIon,1) * ss_factor)
            !  end if

             ! Photoionization by sub-ionizing ISRF --> only impacts lowest ionization states
            !  if (iIon .eq. 1) then
            !     de = de + (UV_background_G0 * elements(iElement)%G0_photo_rate)
            !  end if

             ! Recombination
             if (iIon .gt. 1) then 
                de = de + (recombination(TK, iIon, iElement) * ne)
             end if 

             ! Cosmic ray ionization
            !  if (iIon .lt. n_ions) then
            !     de = de + (cosmic_ray_ionization_rates(iElement,iIon) * total_cosmic_ray_ionization_rate)
            !  end if

             ! Cosmic ray ionization from induced UV
            !  if (iIon .eq. 1) then 
            !     de = de + (cosmic_ray_ionization_rates_induced_UV(iElement) * cosmic_ray_scale_factor)
            !  end if

             ! Recombination on dust
            !  if (iIon .gt. 1) then 
            !     de = de + (dust_recombination(iIon, iElement, TK, UV_background_G0, ne) * dust_effective_number_density)
            !  end if

             !TODO(code): add destruction from local radiation

             !/////////////////////////
             !//   Charge Transfer   //
             !/////////////////////////
             ! Note, this was split off due to cross species
             !coupling. Saves us an extra double loop
             !Charge exchange
            !  if (iElement.eq.1) then !If element is hydrogen
            !     ! Loop over all other elements
            !     do i_other_Element = 2,n_elements
            !        if (elements(i_other_Element)%atomic_number > 0) then
            !           ! Loop over all other ionization states
            !           do i_other_Ion = 1,elements(i_other_Element)%n_ions
            !              ! Get the number density of the other ion
            !              paired_ion_number_density = nElement(i_other_Element, icell) * dXion(i_other_Element,i_other_Ion)

            !              if (iIon.eq.1) then !H
            !                 cr = cr + (charge_transfer_ionization(i_other_Ion,i_other_Element,TK) * dXion(iElement,iIon+1) * paired_ion_number_density) ! Example:  O + H+ => O+ + H
            !                 de = de + (charge_transfer_recombination(i_other_Ion,i_other_Element,TK) * paired_ion_number_density) ! Example:  O+ + H => O + H+
            !              else !H+
            !                 de = de + (charge_transfer_ionization(i_other_Ion,i_other_Element,TK) * paired_ion_number_density) ! Example:  O + H+ => O+ + H
            !                 cr = cr + (charge_transfer_recombination(i_other_Ion,i_other_Element,TK) * dXion(iElement,iIon-1) * paired_ion_number_density) ! Example:  O+ + H => O + H+
            !              end if

            !           end do !end loop over ionization states
            !        end if
            !     end do ! end loop over other elements
            !  else
            !     HI_number_density = dXion(1,1) * nElement(1, icell)
            !     HII_number_density = dXion(1,2) * nElement(1, icell)

            !     if (iIon.gt.1) then
            !        ! Ionization from less excited state
            !        cr = cr + (charge_transfer_ionization(iIon-1,iElement,TK) * dXion(iElement,iIon-1) * HII_number_density) ! Example:  O + H+ => O+ + H   

            !        ! Charge exchange recombination 
            !        de = de + (charge_transfer_recombination(iIon,iElement,TK) * HI_number_density) ! Example:  O+ + H => O + H+
            !     end if

            !     if (iIon.lt.n_ions) then
            !        ! Charge exchange ionization
            !        de = de + (charge_transfer_ionization(iIon,iElement,TK) * HII_number_density) ! Example:  O + H+ => O+ + H   

            !        ! Charge exchange recombination from the more excited state
            !        cr = cr + (charge_transfer_recombination(iIon+1,iElement,TK) * dXion(iElement,iIon+1) * HI_number_density) ! Example:  O+ + H => O + H+
            !     end if
            !  end if

             !/////////////////////////
             !//       Update        //
             !/////////////////////////
             dXion(iElement,iIon) = (cr*ddt(icell) + dXion(iElement,iIon))/(1.d0 + de*ddt(icell))
             dXion(iElement,iIon) = min(max(dXion(iElement,iIon),x_MIN),1.d0)

             ! Get the new electron fraction
             ne = getNe(dXion, nElement(:,icell))
             xe = ne / nElement(1, icell)
             phi_s = secondary_cr_rates(xe)
             primary_cosmic_ray_ionization_rate = total_cosmic_ray_ionization_rate / (1.d0 + phi_s)

             ! Check for convergence -- Fractional change in ion
             dUU = MAX(dUU,ABS((dXion(iElement,iIon)-xion(iElement,iIon,icell))/(xion(iElement,iIon,icell)+x_FM)))
             dUU = dUU * one_over_x_FRAC
             if(dUU .gt. 1.) then
                code=6 !TODO(code) update this code for each ion
                RETURN
             end if

             fracMax=MAX(fracMax,dUU)

             ! Check for convergence -- Fractional change in electrons
             dUU=ABS((ne-neInit)) / (neInit+x_FM) * one_over_x_FRAC
             if(dUU .gt. 1.) then
                code=8
                RETURN
             endif

             fracMax=MAX(fracMax,dUU)

          end do ! END ION LOOP

          ! REDUCE SO THAT IONIZATION FRACTIONS SUM TO 1
          ! Get the indices of the number of ions + number of molecules
          ion_fracs = elements(iElement)%n_ions + elements(iElement)%n_mol

          ! Sum the ionization mass fractions
          current_mass_frac = sum(dXion(iElement,1:ion_fracs))

          ! Get the index of the maximum 
          idx_max = MAXLOC(dXion(iElement,1:ion_fracs), DIM=1)

          ! Update the max ioniztion state so that it sums to 1
          dXion(iElement,idx_max) = 1.d0 - (current_mass_frac - dXion(iElement,idx_max))
       end if
    end do ! END ELEMENT LOOP

    ! UPDATE CONSTANT T WITH NEW MU **************************************
    if(r%neq_isTconst)then
       mu = getMu(r, dXion, dT2)
       dT2 = r%neq_Tconst/mu
    endif

    ! CLEAN UP AND RETURN ************************************************
    dT2 = dT2-T2(icell) ; dXion(:,:) = dXion(:,:)-xion(:,:,icell)
#ifdef RT
    dNp(:) = dNp(:)-Np(:,icell) ; dFp(:,:) = dFp(:,:)-Fp(:,:,icell)
    dp_gas(:)= dp_gas(:)-p_gas(:,icell)
#endif
    ! Now the dUs are really changes, not new values
    ! Check if we are safe to use a bigger timestep in next iteration:
    ! TODO(code) update with smarter timestep criteria
    if(fracMax .lt. 0.5) then
       dt_rec=ddt(icell)*2.
    else
       dt_rec=ddt(icell)
    endif
    dt_ok=.true.
    code=0

  END SUBROUTINE rtz_cool_step

END SUBROUTINE rtz_solve_cooling

!************************************************************************
SUBROUTINE rtz_updateRTGroups_CoolConstants(r,tables)
  ! TODO(code): update this for the metals
  ! Update photon group cooling and heating constants, to reflect an update
  ! in rt_c_cgs and in the cross-sections and energies in the groups.
  !------------------------------------------------------------------------
  implicit none
  type(run_t)::r
  type(neq_cooling_t)::tables
#ifdef RT
  !------------------------------------------------------------------------
  integer::iP, iI, i
  !------------------------------------------------------------------------
  do i=r%nlevelmax,r%levelmin,-1
    tables%signc(:,:,i) = r%group_csn*tables%rt_c_cgs(i)        ! [cm3 s-1]
    tables%sigec(:,:,i) = r%group_cse*tables%rt_c_cgs(i)        ! [cm3 s-1]
    do iP = 1,nrtgrp
      do iI = 1,nion               ! Photoheating rates for photons on ions
        tables%PHrate(iP,iI,i) =  eV2erg * &    ! See eq (19) in Aubert(08)
             (tables%sigec(iP,iI,i) * r%group_egy(iP)  &
             -tables%signc(iP,iI,i)*r%ionEvs(iI))
        tables%PHrate(iP,iI,i) = max(tables%PHrate(iP,iI,i),0d0)!Heating>0
      end do
    end do
  end do
#endif
END SUBROUTINE rtz_updateRTGroups_CoolConstants
!************************************************************************
FUNCTION getNe(xion,nion) result(ne)
  !Returns the electron number density by looping over all elements 
  !and summing their contributions
  implicit none

  real(KIND=8), intent(in)::xion(1:n_elements,1:n_elements)
  real(KIND=8), intent(in)::nion(1:n_elements)
  real(KIND=8)::ne
  integer::iIons, iElement, n_ions

  ne = 0.d0

  !Loop over all elements
  do iElement=1,n_elements
     if (elements(iElement)%atomic_number .gt. 0) then
        !Get the number of ions
        n_ions = elements(iElement)%n_ions

        !Loop over all ions 
        !Start loop at 2, no electrons in the ground state
        do iIons=2,n_ions
           ne = ne + (nion(iElement) * xion(iElement,iIons) * real(iIons - 1, kind=8)) 
        end do
     end if
  end do

END FUNCTION getNe

FUNCTION secondary_cr_rates(xe) result(phi_s)
  ! Secondary CR ionization rate
  implicit none

  real(KIND=8), intent(in):: xe
  real(KIND=8):: phi_s

  phi_s = (1.d0 - (xe / 1.2d0)) * (0.67d0 / (1.d0 + (xe / 0.05d0)))

END FUNCTION secondary_cr_rates

END MODULE rtz_cooling_module