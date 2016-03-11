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
  integer,parameter::nx=2
  integer,parameter::iu1=-1,iu2=nx+2
  integer,parameter::ju1=(1-ndim/2)-1*(ndim/2),ju2=(1-ndim/2)+(nx+2)*(ndim/2)
  integer,parameter::ku1=(1-ndim/3)-1*(ndim/3),ku2=(1-ndim/3)+(nx+2)*(ndim/3)
  integer,parameter::if1=1,if2=nx+1
  integer,parameter::jf1=1,jf2=(1-ndim/2)+(nx+1)*(ndim/2)
  integer,parameter::kf1=1,kf2=(1-ndim/3)+(nx+1)*(ndim/3)
  integer,parameter::io1=0,io2=nx/2+1
  integer,parameter::jo1=(1-ndim/2),jo2=(1-ndim/2)+(nx/2+1)*(ndim/2)
  integer,parameter::ko1=(1-ndim/3),ko2=(1-ndim/3)+(nx/2+1)*(ndim/3)
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
  integer ,dimension(:,:,:,:),allocatable::nborloc

  ! Work space for hydro kernel
  integer,parameter::nx_2=4
  integer,parameter::iu1_2=-1,iu2_2=nx_2+2
  integer,parameter::ju1_2=(1-ndim/2)-1*(ndim/2),ju2_2=(1-ndim/2)+(nx_2+2)*(ndim/2)
  integer,parameter::ku1_2=(1-ndim/3)-1*(ndim/3),ku2_2=(1-ndim/3)+(nx_2+2)*(ndim/3)
  integer,parameter::if1_2=1,if2_2=nx_2+1
  integer,parameter::jf1_2=1,jf2_2=(1-ndim/2)+(nx_2+1)*(ndim/2)
  integer,parameter::kf1_2=1,kf2_2=(1-ndim/3)+(nx_2+1)*(ndim/3)
  integer,parameter::io1_2=0,io2_2=nx_2/2+1
  integer,parameter::jo1_2=(1-ndim/2),jo2_2=(1-ndim/2)+(nx_2/2+1)*(ndim/2)
  integer,parameter::ko1_2=(1-ndim/3),ko2_2=(1-ndim/3)+(nx_2/2+1)*(ndim/3)
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
  integer ,dimension(:,:,:,:),allocatable::nborloc_2

  ! Work space for hydro kernel
  integer,parameter::nx_4=8
  integer,parameter::iu1_4=-1,iu2_4=nx_4+2
  integer,parameter::ju1_4=(1-ndim/2)-1*(ndim/2),ju2_4=(1-ndim/2)+(nx_4+2)*(ndim/2)
  integer,parameter::ku1_4=(1-ndim/3)-1*(ndim/3),ku2_4=(1-ndim/3)+(nx_4+2)*(ndim/3)
  integer,parameter::if1_4=1,if2_4=nx_4+1
  integer,parameter::jf1_4=1,jf2_4=(1-ndim/2)+(nx_4+1)*(ndim/2)
  integer,parameter::kf1_4=1,kf2_4=(1-ndim/3)+(nx_4+1)*(ndim/3)
  integer,parameter::io1_4=0,io2_4=nx_4/2+1
  integer,parameter::jo1_4=(1-ndim/2),jo2_4=(1-ndim/2)+(nx_4/2+1)*(ndim/2)
  integer,parameter::ko1_4=(1-ndim/3),ko2_4=(1-ndim/3)+(nx_4/2+1)*(ndim/3)
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
  integer ,dimension(:,:,:,:),allocatable::nborloc_4

  ! Work space for hydro kernel
  integer,parameter::nx_8=16
  integer,parameter::iu1_8=-1,iu2_8=nx_8+2
  integer,parameter::ju1_8=(1-ndim/2)-1*(ndim/2),ju2_8=(1-ndim/2)+(nx_8+2)*(ndim/2)
  integer,parameter::ku1_8=(1-ndim/3)-1*(ndim/3),ku2_8=(1-ndim/3)+(nx_8+2)*(ndim/3)
  integer,parameter::if1_8=1,if2_8=nx_8+1
  integer,parameter::jf1_8=1,jf2_8=(1-ndim/2)+(nx_8+1)*(ndim/2)
  integer,parameter::kf1_8=1,kf2_8=(1-ndim/3)+(nx_8+1)*(ndim/3)
  integer,parameter::io1_8=0,io2_8=nx_8/2+1
  integer,parameter::jo1_8=(1-ndim/2),jo2_8=(1-ndim/2)+(nx_8/2+1)*(ndim/2)
  integer,parameter::ko1_8=(1-ndim/3),ko2_8=(1-ndim/3)+(nx_8/2+1)*(ndim/3)
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
  integer ,dimension(:,:,:,:),allocatable::nborloc_8

  ! Work space for hydro kernel
  integer,parameter::nx_16=32
  integer,parameter::iu1_16=-1,iu2_16=nx_16+2
  integer,parameter::ju1_16=(1-ndim/2)-1*(ndim/2),ju2_16=(1-ndim/2)+(nx_16+2)*(ndim/2)
  integer,parameter::ku1_16=(1-ndim/3)-1*(ndim/3),ku2_16=(1-ndim/3)+(nx_16+2)*(ndim/3)
  integer,parameter::if1_16=1,if2_16=nx_16+1
  integer,parameter::jf1_16=1,jf2_16=(1-ndim/2)+(nx_16+1)*(ndim/2)
  integer,parameter::kf1_16=1,kf2_16=(1-ndim/3)+(nx_16+1)*(ndim/3)
  integer,parameter::io1_16=0,io2_16=nx_16/2+1
  integer,parameter::jo1_16=(1-ndim/2),jo2_16=(1-ndim/2)+(nx_16/2+1)*(ndim/2)
  integer,parameter::ko1_16=(1-ndim/3),ko2_16=(1-ndim/3)+(nx_16/2+1)*(ndim/3)
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
  integer ,dimension(:,:,:,:),allocatable::nborloc_16

  ! Work space for hydro kernel
  integer,parameter::nx_32=64
  integer,parameter::iu1_32=-1,iu2_32=nx_32+2
  integer,parameter::ju1_32=(1-ndim/2)-1*(ndim/2),ju2_32=(1-ndim/2)+(nx_32+2)*(ndim/2)
  integer,parameter::ku1_32=(1-ndim/3)-1*(ndim/3),ku2_32=(1-ndim/3)+(nx_32+2)*(ndim/3)
  integer,parameter::if1_32=1,if2_32=nx_32+1
  integer,parameter::jf1_32=1,jf2_32=(1-ndim/2)+(nx_32+1)*(ndim/2)
  integer,parameter::kf1_32=1,kf2_32=(1-ndim/3)+(nx_32+1)*(ndim/3)
  integer,parameter::io1_32=0,io2_32=nx_32/2+1
  integer,parameter::jo1_32=(1-ndim/2),jo2_32=(1-ndim/2)+(nx_32/2+1)*(ndim/2)
  integer,parameter::ko1_32=(1-ndim/3),ko2_32=(1-ndim/3)+(nx_32/2+1)*(ndim/3)
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
  integer ,dimension(:,:,:,:),allocatable::nborloc_32

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

