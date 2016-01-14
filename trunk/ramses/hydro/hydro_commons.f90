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

  ! Work space for hydro kernel
  integer::iu1,iu2,ju1,ju2,ku1,ku2 ! Hydro cell-grid size
  integer::if1,if2,jf1,jf2,kf1,kf2 ! Hydro face-grid size
  integer::io1,io2,jo1,jo2,ko1,ko2 ! Hydro oct-grid size
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
  logical ,dimension(:,:,:),allocatable::okloc
  integer ,dimension(:,:,:),allocatable::childloc
  integer ,dimension(:,:,:),allocatable::parentloc

  ! Work space for hydro kernel
  integer::iu1_2,iu2_2,ju1_2,ju2_2,ku1_2,ku2_2 ! Hydro cell-grid size
  integer::if1_2,if2_2,jf1_2,jf2_2,kf1_2,kf2_2 ! Hydro face-grid size
  integer::io1_2,io2_2,jo1_2,jo2_2,ko1_2,ko2_2 ! Hydro oct-grid size
  real(dp),dimension(:,:,:,:),allocatable::uloc_2
  real(dp),dimension(:,:,:,:),allocatable::gloc_2
  real(dp),dimension(:,:,:,:),allocatable::qloc_2
  real(dp),dimension(:,:,:),allocatable::cloc_2
  real(dp),dimension(:,:,:,:,:),allocatable::flux_2
  real(dp),dimension(:,:,:,:,:),allocatable::tmp_2
  real(dp),dimension(:,:,:,:,:),allocatable::dq_2
  real(dp),dimension(:,:,:,:,:),allocatable::qm_2
  real(dp),dimension(:,:,:,:,:),allocatable::qp_2
  real(dp),dimension(:,:,:,:),allocatable::fx_2
  real(dp),dimension(:,:,:,:),allocatable::tx_2
  real(dp),dimension(:,:,:),allocatable::divu_2
  logical ,dimension(:,:,:),allocatable::okloc_2
  integer ,dimension(:,:,:),allocatable::childloc_2
  integer ,dimension(:,:,:),allocatable::parentloc_2

  ! Work space for hydro kernel
  integer::iu1_4,iu2_4,ju1_4,ju2_4,ku1_4,ku2_4 ! Hydro cell-grid size
  integer::if1_4,if2_4,jf1_4,jf2_4,kf1_4,kf2_4 ! Hydro face-grid size
  integer::io1_4,io2_4,jo1_4,jo2_4,ko1_4,ko2_4 ! Hydro oct-grid size
  real(dp),dimension(:,:,:,:),allocatable::uloc_4
  real(dp),dimension(:,:,:,:),allocatable::gloc_4
  real(dp),dimension(:,:,:,:),allocatable::qloc_4
  real(dp),dimension(:,:,:),allocatable::cloc_4
  real(dp),dimension(:,:,:,:,:),allocatable::flux_4
  real(dp),dimension(:,:,:,:,:),allocatable::tmp_4
  real(dp),dimension(:,:,:,:,:),allocatable::dq_4
  real(dp),dimension(:,:,:,:,:),allocatable::qm_4
  real(dp),dimension(:,:,:,:,:),allocatable::qp_4
  real(dp),dimension(:,:,:,:),allocatable::fx_4
  real(dp),dimension(:,:,:,:),allocatable::tx_4
  real(dp),dimension(:,:,:),allocatable::divu_4
  logical ,dimension(:,:,:),allocatable::okloc_4
  integer ,dimension(:,:,:),allocatable::childloc_4
  integer ,dimension(:,:,:),allocatable::parentloc_4

  ! Work space for hydro kernel
  integer::iu1_8,iu2_8,ju1_8,ju2_8,ku1_8,ku2_8 ! Hydro cell-grid size
  integer::if1_8,if2_8,jf1_8,jf2_8,kf1_8,kf2_8 ! Hydro face-grid size
  integer::io1_8,io2_8,jo1_8,jo2_8,ko1_8,ko2_8 ! Hydro oct-grid size
  real(dp),dimension(:,:,:,:),allocatable::uloc_8
  real(dp),dimension(:,:,:,:),allocatable::gloc_8
  real(dp),dimension(:,:,:,:),allocatable::qloc_8
  real(dp),dimension(:,:,:),allocatable::cloc_8
  real(dp),dimension(:,:,:,:,:),allocatable::flux_8
  real(dp),dimension(:,:,:,:,:),allocatable::tmp_8
  real(dp),dimension(:,:,:,:,:),allocatable::dq_8
  real(dp),dimension(:,:,:,:,:),allocatable::qm_8
  real(dp),dimension(:,:,:,:,:),allocatable::qp_8
  real(dp),dimension(:,:,:,:),allocatable::fx_8
  real(dp),dimension(:,:,:,:),allocatable::tx_8
  real(dp),dimension(:,:,:),allocatable::divu_8
  logical ,dimension(:,:,:),allocatable::okloc_8
  integer ,dimension(:,:,:),allocatable::childloc_8
  integer ,dimension(:,:,:),allocatable::parentloc_8

  ! Work space for hydro kernel
  integer::iu1_16,iu2_16,ju1_16,ju2_16,ku1_16,ku2_16 ! Hydro cell-grid size
  integer::if1_16,if2_16,jf1_16,jf2_16,kf1_16,kf2_16 ! Hydro face-grid size
  integer::io1_16,io2_16,jo1_16,jo2_16,ko1_16,ko2_16 ! Hydro oct-grid size
  real(dp),dimension(:,:,:,:),allocatable::uloc_16
  real(dp),dimension(:,:,:,:),allocatable::gloc_16
  real(dp),dimension(:,:,:,:),allocatable::qloc_16
  real(dp),dimension(:,:,:),allocatable::cloc_16
  real(dp),dimension(:,:,:,:,:),allocatable::flux_16
  real(dp),dimension(:,:,:,:,:),allocatable::tmp_16
  real(dp),dimension(:,:,:,:,:),allocatable::dq_16
  real(dp),dimension(:,:,:,:,:),allocatable::qm_16
  real(dp),dimension(:,:,:,:,:),allocatable::qp_16
  real(dp),dimension(:,:,:,:),allocatable::fx_16
  real(dp),dimension(:,:,:,:),allocatable::tx_16
  real(dp),dimension(:,:,:),allocatable::divu_16
  logical ,dimension(:,:,:),allocatable::okloc_16
  integer ,dimension(:,:,:),allocatable::childloc_16
  integer ,dimension(:,:,:),allocatable::parentloc_16

  ! Work space for hydro kernel
  integer::iu1_32,iu2_32,ju1_32,ju2_32,ku1_32,ku2_32 ! Hydro cell-grid size
  integer::if1_32,if2_32,jf1_32,jf2_32,kf1_32,kf2_32 ! Hydro face-grid size
  integer::io1_32,io2_32,jo1_32,jo2_32,ko1_32,ko2_32 ! Hydro oct-grid size
  real(dp),dimension(:,:,:,:),allocatable::uloc_32
  real(dp),dimension(:,:,:,:),allocatable::gloc_32
  real(dp),dimension(:,:,:,:),allocatable::qloc_32
  real(dp),dimension(:,:,:),allocatable::cloc_32
  real(dp),dimension(:,:,:,:,:),allocatable::flux_32
  real(dp),dimension(:,:,:,:,:),allocatable::tmp_32
  real(dp),dimension(:,:,:,:,:),allocatable::dq_32
  real(dp),dimension(:,:,:,:,:),allocatable::qm_32
  real(dp),dimension(:,:,:,:,:),allocatable::qp_32
  real(dp),dimension(:,:,:,:),allocatable::fx_32
  real(dp),dimension(:,:,:,:),allocatable::tx_32
  real(dp),dimension(:,:,:),allocatable::divu_32
  logical ,dimension(:,:,:),allocatable::okloc_32
  integer ,dimension(:,:,:),allocatable::childloc_32
  integer ,dimension(:,:,:),allocatable::parentloc_32

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

