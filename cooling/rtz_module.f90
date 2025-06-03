module rtz_module
  implicit none

  type Element
      integer(KIND=4) :: atomic_number
      integer(KIND=4) :: n_ions
      integer(KIND=4) :: n_mol
      real(KIND=8)   :: atomic_mass
      real(KIND=8)   :: z_solar
      real(KIND=8)   :: G0_photo_rate
      real(KIND=8)   :: depletion
  end type Element

end module rtz_module