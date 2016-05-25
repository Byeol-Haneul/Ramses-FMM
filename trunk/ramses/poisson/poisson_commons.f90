module poisson_commons 
  use amr_commons
  use hash
  use poisson_parameters

  ! Minimum MG level
  integer :: levelmin_mg

  ! MG hash table
  type(hash_table)::mg_dict

  ! MG grid
  integer,allocatable,dimension(:)::head_mg,tail_mg,noct_mg,noct_tot_mg
  integer::ifree_mg

  ! Multigrid safety switch
  logical, allocatable, dimension(:) :: safe_mode

  ! Multipole coefficients
  real(dp),dimension(1:ndim+1)::multipole

  ! Conjugate Gradient communicator
  integer::send_tot,recv_tot
  integer,dimension(:),allocatable::send_cnt,recv_cnt,send_oft,recv_oft
  integer,dimension(:),allocatable::grid_recv_buf
  integer,dimension(:,:),allocatable::nbor_indx
  real(dp),dimension(:),allocatable::phi_send_buf,phi_recv_buf
  real(dp),dimension(:,:),allocatable::phi_remote

  ! Temporary workspace
  integer,dimension(:),allocatable::nremote
  integer,dimension(:,:),allocatable::nalltoall,nalltoall_tot
  integer,dimension(:),allocatable::x_send_buf,y_send_buf,z_send_buf
  integer,dimension(:),allocatable::x_recv_buf,y_recv_buf,z_recv_buf

  ! Multigrid communicator
  type comm_mg
     integer::send_tot
     integer::recv_tot
     integer ,dimension(:,:),pointer::nbor_indx
     real(dp),dimension(:,:),pointer::phi_remote
     real(dp),dimension(:,:),pointer::dis_remote
     real(dp),dimension(:)  ,pointer::phi_send_buf
     real(dp),dimension(:)  ,pointer::phi_recv_buf
     integer ,dimension(:)  ,pointer::grid_recv_buf
     integer ,dimension(:)  ,pointer::send_cnt
     integer ,dimension(:)  ,pointer::send_oft
     integer ,dimension(:)  ,pointer::recv_cnt
     integer ,dimension(:)  ,pointer::recv_oft
  end type comm_mg
  type(comm_mg),allocatable,dimension(:)::buffer_mg

end module poisson_commons
