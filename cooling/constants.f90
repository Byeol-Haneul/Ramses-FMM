module constants

  ! Numerical constants
  real(kind=8),parameter ::twopi        = 6.2831853d0
  real(kind=8),parameter ::pi           = twopi/2d0

  ! Physical constants
  ! Source:
  ! * SI - SI Brochure (2018)
  ! * PCAD - http://www.astro.wisc.edu/~dolan/constants.html
  ! * NIST - National Institute of Standards and Technology
  ! * IAU - Internatonal Astronomical Union resolution
  real(kind=8),parameter ::hplanck      = 6.6260702d-27 ! Planck const. [erg s]; SI
  real(kind=8),parameter ::eV2erg       = 1.6021766d-12 ! Electronvolt [erg]; SI
  real(kind=8),parameter ::kB           = 1.3806490d-16 ! Boltzmann const. [erg K-1]; SI
  real(kind=8),parameter ::clight       = 2.9979246d+10 ! Speed of light [cm s-1]; SI
  real(kind=8),parameter ::a_r          = 7.5657233d-15 ! Radiation density const. [erg cm-3 K-4]; SI (derived)
  real(kind=8),parameter ::mH           = 1.6605390d-24 ! H atom mass [g] = amu, i.e. atomic mass unit; NIST
  real(kind=8),parameter ::gravconst    = 6.6740800d-08 ! Gravitational const. [cm3 kg-1 s-2]; NIST
  real(kind=8),parameter ::sigma_T      = 6.6524587d-25 ! Thomson scattering cross-section [cm2]; NIST
  real(kind=8),parameter ::M_sun        = 1.9891000d+33 ! Solar Mass [g]; IAU
  real(kind=8),parameter ::L_sun        = 3.8280000d+33 ! Solar Lum [erg s-1]; IAU
  real(kind=8),parameter ::rhoc         = 1.8800000d-29 ! Crit. density [g cm-3]
  real(kind=8),parameter ::one_over_clight=3.335640484668562d-11 ! Save some computation
  ! Conversion factors - distance
  ! IAU 2012 convention:
  ! 1 pc = 648000 AU / pi
  ! 1 AU = 14 959 787 070 000 cm
  real(kind=8),parameter ::pc2cm        = 3.0856776d+18
  real(kind=8),parameter ::kpc2cm       = 3.0856776d+21
  real(kind=8),parameter ::Mpc2cm       = 3.0856776d+24
  real(kind=8),parameter ::Gpc2cm       = 3.0856776d+27

  ! Conversion factors - time
  ! Year definition follows IAU recommendation
  ! https://www.iau.org/publications/proceedings_rules/units/
  real(kind=8),parameter ::yr2sec       = 3.15576000d+07 ! Year [s]
  real(kind=8),parameter ::kyr2sec      = 3.15576000d+10 ! Kyr [s]
  real(kind=8),parameter ::Myr2sec      = 3.15576000d+13 ! Myr [s]
  real(kind=8),parameter ::Gyr2sec      = 3.15576000d+16 ! Gyr [s]
  real(kind=8),parameter ::sec2Gyr      = 3.16880878d-17 ! sec [Gyr]

  ! Non-equilibrium cooling constants
  real(kind=8),parameter::ionEv_HI    = 11.20d0
  real(kind=8),parameter::ionEv_HII   = 13.60d0
  real(kind=8),parameter::ionEv_HeII  = 24.59d0
  real(kind=8),parameter::ionEv_HeIII = 54.42d0
  
end module constants
