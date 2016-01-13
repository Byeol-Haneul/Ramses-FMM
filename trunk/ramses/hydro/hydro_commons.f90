module hydro_commons
  use amr_parameters
  use hydro_parameters
  real(dp)::mass_tot=0.0D0,mass_tot_0=0.0D0
  real(dp)::ana_xmi,ana_xma,ana_ymi,ana_yma,ana_zmi,ana_zma
  integer::nbins
  real(kind=8),parameter ::twopi   = 6.2831853d0
  real(kind=8),parameter ::hplanck = 6.6262000d-27
  real(kind=8),parameter ::eV      = 1.6022000d-12
  real(kind=8),parameter ::kB      = 1.3806200d-16
  real(kind=8),parameter ::clight  = 2.9979250d+10
  real(kind=8),parameter ::Gyr     = 3.1536000d+16
  real(kind=8)           ::X       = 0.76
  real(kind=8)           ::Y       = 0.24
  real(kind=8),parameter ::rhoc    = 1.8800000d-29
  real(kind=8),parameter ::mH      = 1.6600000d-24

  integer::iu1,iu2,ju1,ju2,ku1,ku2
  integer::if1,if2,jf1,jf2,kf1,kf2

  ! Work space for hydro kernel
  real(dp),dimension(:,:,:,:),allocatable::uloc
  real(dp),dimension(:,:,:,:),allocatable::gloc
  real(dp),dimension(:,:,:,:),allocatable::qloc
  real(dp),dimension(:,:,:),allocatable::cloc
  real(dp),dimension(:,:,:,:,:),allocatable::flux
  real(dp),dimension(:,:,:,:,:),allocatable::tmp
  real(dp),dimension(:,:,:,:,:),allocatable::dq
  real(dp),dimension(:,:,:,:,:),allocatable::qm
  real(dp),dimension(:,:,:,:,:),allocatable::qp
  real(dp),dimension(:,:,:,:),allocatable::fx
  real(dp),dimension(:,:,:,:),allocatable::tx
  real(dp),dimension(:,:,:),allocatable::divu
  logical ,dimension(:,:,:),allocatable::okoc

end module hydro_commons

module const
  use amr_parameters
  real(dp)::bigreal = 1.0e+30
  real(dp)::zero = 0.0
  real(dp)::one = 1.0
  real(dp)::two = 2.0
  real(dp)::three = 3.0
  real(dp)::four = 4.0
  real(dp)::two3rd = 0.6666666666666667
  real(dp)::half = 0.5
  real(dp)::third = 0.33333333333333333
  real(dp)::forth = 0.25
  real(dp)::sixth = 0.16666666666666667
end module const

