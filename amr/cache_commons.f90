module cache_commons
  use amr_parameters, only: dp,ndim,twotondim
  use hydro_parameters, only: nvar
  use call_back
  
  ! Communication-related objects
  integer::mail_counter=0
  integer::request_id,flush_id
  integer,dimension(:),allocatable::reply_id

  integer::flush_tag=1000,msg_tag=100,request_tag=10

  ! Software cache parameters
  integer,parameter::operation_initflag=1,operation_upload=2,operation_godunov=3,operation_smooth=4
  integer,parameter::operation_hydro=5,operation_refine=6,operation_derefine=7,operation_loadbalance=8
  integer,parameter::operation_phi=9,operation_rho=10,operation_multipole=11,operation_cg=12
  integer,parameter::operation_build_mg=13,operation_restrict_mask=14,operation_mg=15
  integer,parameter::operation_restrict_res=16,operation_scan=17,operation_interpol=18
  integer,parameter::operation_split=19,operation_kick=20

  integer,parameter::domain_decompos_amr=1,domain_decompos_mg=2
  
  ! Message size
  integer,parameter::ntilemax=16  ! Fetch message buffer size
  integer,parameter::nflushmax=128  ! Flush message buffer size
  
  integer::size_msg_array ! Cache element size
  integer::size_request_array ! Comm. buffer sizes
  integer::size_flush_array
  integer::size_fetch_array

  ! Combiner rules
  integer,parameter::COMBINER_EXIST=1
  integer,parameter::COMBINER_CREATE=2
  
  integer::combiner_rule

  ! Message arrays
  integer(kind=4),dimension(:),allocatable::recv_request_array
  integer(kind=4),dimension(:),allocatable::send_request_array
  integer(kind=4),dimension(:),allocatable::recv_fetch_array
  integer(kind=4),dimension(:),allocatable::send_fetch_array
  integer(kind=4),dimension(:),allocatable::recv_flush_array
  integer(kind=4),dimension(:),allocatable::send_flush_array
  
  ! Message element type
  type msg_int4
     integer(kind=4),dimension(1:twotondim)::int4
  end type msg_int4
  type msg_realdp
     integer(kind=4),dimension(1:twotondim)::int4
     real(kind=dp),dimension(1:twotondim,1:nvar)::realdp
  end type msg_realdp
  type msg_small_realdp
     real(kind=dp),dimension(1:twotondim)::realdp
  end type msg_small_realdp
  type msg_twin_realdp
     real(kind=dp),dimension(1:twotondim)::realdp_phi
     real(kind=dp),dimension(1:twotondim)::realdp_dis
  end type msg_twin_realdp
  type msg_three_realdp
     real(kind=dp),dimension(1:twotondim)::realdp_phi
     real(kind=dp),dimension(1:twotondim)::realdp_phi_old
     real(kind=dp),dimension(1:twotondim)::realdp_dis
  end type msg_three_realdp
  type msg_large_realdp
     integer(kind=4),dimension(1:twotondim)::int4
     real(kind=dp),dimension(1:twotondim,1:nvar)::realdp_hydro
     real(kind=dp),dimension(1:twotondim,1:ndim+2)::realdp_poisson
  end type msg_large_realdp

  ! Cache call back functions
  type(cache_f)::pack_fetch,unpack_fetch
  type(cache_f)::pack_flush,unpack_flush
  type(cache_f)::init_flush
  
end module cache_commons
