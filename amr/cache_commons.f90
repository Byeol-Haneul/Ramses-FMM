module cache_commons
  use amr_parameters, only: dp,ndim,twotondim
  use hydro_parameters, only: nvar

  ! Communication-related objects
  integer::mail_counter=0
  integer::request_id,flush_id
  integer,allocatable,dimension(:)::reply_id
  logical::expect_flush=.false.

  ! Software cache parameters
  integer::cache_operation
  integer::operation_initflag=1,operation_upload=2,operation_godunov=3,operation_smooth=4
  integer::operation_hydro=5,operation_refine=6,operation_derefine=7,operation_loadbalance=8
  integer::operation_phi=9,operation_rho=10,operation_multipole=11,operation_cg=12
  integer::operation_build_mg=13,operation_restrict_mask=14,operation_mg=15
  integer::operation_restrict_res=16,operation_scan=17,operation_interpol=18
  integer::operation_split=19,operation_kick=20
  integer::cache_operation_type
  integer::operation_type_flag=1,operation_type_hydro=2,operation_type_poisson=3  
  integer::operation_type_refine=4,operation_type_mg=5,operation_type_interpol=6
  integer::domain_decompos_amr=1,domain_decompos_mg=2

  ! New MPI derived types
  integer::new_mpi_int4_msg,new_mpi_realdp_msg,new_mpi_request
  integer::new_mpi_small_realdp_msg,new_mpi_large_realdp_msg
  integer::new_mpi_twin_realdp_msg,new_mpi_three_realdp_msg
  integer::new_mpi_int4_flush,new_mpi_realdp_flush
  integer::new_mpi_small_realdp_flush,new_mpi_large_realdp_flush
  integer::new_mpi_twin_realdp_flush,new_mpi_three_realdp_flush
  integer::flush_tag=1000,msg_tag=100,request_tag=10

  ! Request message buffer
  type request
     integer(kind=4)::lev
     integer(kind=4),dimension(1:ndim)::ckey
  end type request

  ! Response message buffer
  integer,parameter::ntilemax=16
  type int4_msg
     sequence
     integer(kind=4)::type
     integer(kind=4)::ntile
     integer(kind=4),dimension(1:ntilemax)::lev
     integer(kind=4),dimension(1:ndim,1:ntilemax)::ckey
     integer(kind=4),dimension(1:twotondim,1:ntilemax)::int4
  end type int4_msg
  type realdp_msg
     sequence
     integer(kind=4)::type
     integer(kind=4)::ntile
     integer(kind=4),dimension(1:ntilemax)::lev
     integer(kind=4),dimension(1:ndim,1:ntilemax)::ckey
     integer(kind=4),dimension(1:twotondim,1:ntilemax)::int4
     real(kind=dp),dimension(1:twotondim,1:nvar,1:ntilemax)::realdp
  end type realdp_msg
  type small_realdp_msg
     sequence
     integer(kind=4)::type
     integer(kind=4)::ntile
     integer(kind=4),dimension(1:ntilemax)::lev
     integer(kind=4),dimension(1:ndim,1:ntilemax)::ckey
     real(kind=dp),dimension(1:twotondim,1:ntilemax)::realdp
  end type small_realdp_msg
  type twin_realdp_msg
     sequence
     integer(kind=4)::type
     integer(kind=4)::ntile
     integer(kind=4),dimension(1:ntilemax)::lev
     integer(kind=4),dimension(1:ndim,1:ntilemax)::ckey
     real(kind=dp),dimension(1:twotondim,1:ntilemax)::realdp_phi
     real(kind=dp),dimension(1:twotondim,1:ntilemax)::realdp_dis
  end type twin_realdp_msg
  type three_realdp_msg
     sequence
     integer(kind=4)::type
     integer(kind=4)::ntile
     integer(kind=4),dimension(1:ntilemax)::lev
     integer(kind=4),dimension(1:ndim,1:ntilemax)::ckey
     real(kind=dp),dimension(1:twotondim,1:ntilemax)::realdp_phi
     real(kind=dp),dimension(1:twotondim,1:ntilemax)::realdp_phi_old
     real(kind=dp),dimension(1:twotondim,1:ntilemax)::realdp_dis
  end type three_realdp_msg
  type large_realdp_msg
     sequence
     integer(kind=4)::type
     integer(kind=4)::ntile
     integer(kind=4),dimension(1:ntilemax)::lev
     integer(kind=4),dimension(1:ndim,1:ntilemax)::ckey
     integer(kind=4),dimension(1:twotondim,1:ntilemax)::int4
#ifdef HYDRO
     real(kind=dp),dimension(1:twotondim,1:nvar,1:ntilemax)::realdp_hydro
#endif
#ifdef GRAV
     real(kind=dp),dimension(1:twotondim,1:ndim+2,1:ntilemax)::realdp_poisson
#endif
  end type large_realdp_msg

  ! Fush message buffer
  integer,parameter::nflushmax=128
  type int4_flush
     sequence
     integer(kind=4)::nflush
     integer(kind=4),dimension(1:nflushmax)::lev
     integer(kind=4),dimension(1:ndim,1:nflushmax)::ckey
     integer(kind=4),dimension(1:twotondim,1:nflushmax)::int4
  end type int4_flush
  type realdp_flush
     sequence
     integer(kind=4)::nflush
     integer(kind=4),dimension(1:nflushmax)::lev
     integer(kind=4),dimension(1:ndim,1:nflushmax)::ckey
     integer(kind=4),dimension(1:twotondim,1:nflushmax)::int4
     real(kind=dp),dimension(1:twotondim,1:nvar,1:nflushmax)::realdp
  end type realdp_flush
  type small_realdp_flush
     sequence
     integer(kind=4)::nflush
     integer(kind=4),dimension(1:nflushmax)::lev
     integer(kind=4),dimension(1:ndim,1:nflushmax)::ckey
     real(kind=dp),dimension(1:twotondim,1:nflushmax)::realdp
  end type small_realdp_flush
  type twin_realdp_flush
     sequence
     integer(kind=4)::nflush
     integer(kind=4),dimension(1:nflushmax)::lev
     integer(kind=4),dimension(1:ndim,1:nflushmax)::ckey
     real(kind=dp),dimension(1:twotondim,1:nflushmax)::realdp_phi
     real(kind=dp),dimension(1:twotondim,1:nflushmax)::realdp_dis
  end type twin_realdp_flush
  type three_realdp_flush
     sequence
     integer(kind=4)::nflush
     integer(kind=4),dimension(1:nflushmax)::lev
     integer(kind=4),dimension(1:ndim,1:nflushmax)::ckey
     real(kind=dp),dimension(1:twotondim,1:nflushmax)::realdp_phi
     real(kind=dp),dimension(1:twotondim,1:nflushmax)::realdp_phi_old
     real(kind=dp),dimension(1:twotondim,1:nflushmax)::realdp_dis
  end type three_realdp_flush
  type large_realdp_flush
     sequence
     integer(kind=4)::nflush
     integer(kind=4),dimension(1:nflushmax)::lev
     integer(kind=4),dimension(1:ndim,1:nflushmax)::ckey
     integer(kind=4),dimension(1:twotondim,1:nflushmax)::int4
#ifdef HYDRO
     real(kind=dp),dimension(1:twotondim,1:nvar,1:nflushmax)::realdp_hydro
#endif
#ifdef GRAV
     real(kind=dp),dimension(1:twotondim,1:ndim+2,1:nflushmax)::realdp_poisson
#endif
  end type large_realdp_flush

  ! Communication buffers
  type(request)::recv_request
  type(int4_msg),allocatable,dimension(:)::reply_flag
  type(realdp_msg),allocatable,dimension(:)::reply_hydro
  type(twin_realdp_msg),allocatable,dimension(:)::reply_mg
  type(three_realdp_msg),allocatable,dimension(:)::reply_interpol
  type(small_realdp_msg),allocatable,dimension(:)::reply_poisson
  type(large_realdp_msg),allocatable,dimension(:)::reply_refine
  type(int4_flush)::recv_flush_flag
  type(realdp_flush)::recv_flush_hydro
  type(twin_realdp_flush)::recv_flush_mg
  type(three_realdp_flush)::recv_flush_interpol
  type(small_realdp_flush)::recv_flush_poisson
  type(large_realdp_flush)::recv_flush_refine
  type(int4_flush),allocatable,dimension(:)::send_flush_flag
  type(realdp_flush),allocatable,dimension(:)::send_flush_hydro
  type(twin_realdp_flush),allocatable,dimension(:)::send_flush_mg
  type(three_realdp_flush),allocatable,dimension(:)::send_flush_interpol
  type(small_realdp_flush),allocatable,dimension(:)::send_flush_poisson
  type(large_realdp_flush),allocatable,dimension(:)::send_flush_refine

end module cache_commons
