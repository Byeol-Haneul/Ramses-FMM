subroutine init_amr
  use amr_commons
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif
  integer::ilevel,icpu,info
  integer::intex,realdpex
  integer(kind=8)::max_key
  integer,dimension(1:10)::new_type_disp,new_type_type,new_type_length
  real(kind=4)::real_mem,real_mem_tot

  if(verbose.and.myid==1)write(*,*)'Entering init_amr'

  ! Initial time step for each level
  dtold=0.0D0
  dtnew=0.0D0

  ! Set up cache size
  ncachemax=100000 !0.1*ngridmax

  ! Allocate main oct array
  allocate(grid(1:ngridmax+ncachemax))

  ! Allocate cache-related arrays
  allocate(dirty(1:ncachemax))
  allocate(locked(1:ncachemax))
  allocate(occupied(1:ncachemax))
  allocate(parent_cpu(1:ncachemax))
  dirty=.false.
  occupied=.false.
  locked=.false.
  free_cache=1; ncache=0
  allocate(lev_null(1:ncachemax))
  allocate(ckey_null(1:ndim,1:ncachemax))
  allocate(occupied_null(1:ncachemax))
  occupied_null=.false.
  free_null=1; nnull=0

  ! Allocate oct hash table
  if(verbose.and.myid==1)write(*,*)'Initialize empty hash'
  call init_empty_hash(grid_dict,2*(ngridmax+ncachemax))

  ! Set initial cpu boundaries
  ! Set maximum Cartesian key per level
  if(verbose.and.myid==1)write(*,*)'Initialize level cpu boundaries'
  allocate(bound_key_level(0:ncpu,levelmin:nlevelmax))
  allocate(ckey_max(levelmin:nlevelmax))
  do ilevel=levelmin,nlevelmax
     ckey_max(ilevel)=2**(ilevel-1)
     max_key=2**((ilevel-1)*ndim)
     do icpu=1,ncpu-1
        bound_key_level(icpu,ilevel) = (icpu*max_key)/ncpu
     end do
     bound_key_level(0,ilevel) = 0
     bound_key_level(ncpu,ilevel) = max_key
  end do

  ! Allocate head, tail and numbers for each level
  if(verbose.and.myid==1)write(*,*)'Initialize oct decomposition'
  allocate(head(levelmin:nlevelmax))
  allocate(tail(levelmin:nlevelmax))
  allocate(noct(levelmin:nlevelmax))
  allocate(noct_min(levelmin:nlevelmax))
  allocate(noct_max(levelmin:nlevelmax))
  allocate(noct_tot(levelmin:nlevelmax))
  head=1       ! Head oct in the level
  tail=0       ! Tail oct in the level
  noct=0       ! Number of oct in the level and in the cpu
  noct_tot=0   ! Total number of oct in the level (all cpus)
  noct_min=0   ! Minimum number of oct across all cpus
  noct_max=0   ! Maximum number of oct across all cpus
  noct_used=0  ! Number of oct used across all levels
  noct_used_tot=0  ! Total number of oct used (all cpus)

#ifndef WITHOUTMPI  
  ! Allocate all communication and cache-related variables
  allocate(reply_id(1:ncpu))
  allocate(reply_flag(1:ncpu))
  allocate(reply_hydro(1:ncpu))
  allocate(send_flush_flag(1:ncpu))
  allocate(send_flush_hydro(1:ncpu))
  do icpu=1,ncpu
     send_flush_flag(icpu)%nflush=0
     send_flush_hydro(icpu)%nflush=0
  end do

  ! Create and commit MPI derived types
  call MPI_TYPE_EXTENT(MPI_INTEGER,intex,info)
  call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,realdpex,info)

  ! New type for int4_msg
  new_type_disp(1)=0
  new_type_disp(2)=2*intex
  new_type_disp(3)=(2+ntilemax*(1+ndim))*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_INTEGER
  new_type_length(1)=2
  new_type_length(2)=ntilemax*(1+ndim)
  new_type_length(3)=ntilemax*(twotondim)
  call MPI_TYPE_STRUCT(3,new_type_length,new_type_disp,new_type_type,new_mpi_int4_msg,info)
  call MPI_TYPE_COMMIT(new_mpi_int4_msg,info)

  ! New type for int4_flush
  new_type_disp(1)=0
  new_type_disp(2)=intex
  new_type_disp(3)=(1+nflushmax*(1+ndim))*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_INTEGER
  new_type_length(1)=1
  new_type_length(2)=nflushmax*(1+ndim)
  new_type_length(3)=nflushmax*(twotondim)
  call MPI_TYPE_STRUCT(3,new_type_length,new_type_disp,new_type_type,new_mpi_int4_flush,info)
  call MPI_TYPE_COMMIT(new_mpi_int4_flush,info)

  ! New type for realdp_msg
  new_type_disp(1)=0
  new_type_disp(2)=2*intex
  new_type_disp(3)=(2+ntilemax*(1+ndim))*intex
  new_type_disp(4)=(2+ntilemax*(1+ndim)+ntilemax*(twotondim))*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_INTEGER
  new_type_type(4)=MPI_DOUBLE_PRECISION
  new_type_length(1)=2
  new_type_length(2)=ntilemax*(1+ndim)
  new_type_length(3)=ntilemax*(twotondim)
  new_type_length(4)=ntilemax*(twotondim*nvar)
  call MPI_TYPE_STRUCT(4,new_type_length,new_type_disp,new_type_type,new_mpi_realdp_msg,info)
  call MPI_TYPE_COMMIT(new_mpi_realdp_msg,info)

  ! New type for realdp_flush
  new_type_disp(1)=0
  new_type_disp(2)=intex
  new_type_disp(3)=(1+nflushmax*(1+ndim))*intex
  new_type_disp(4)=(1+nflushmax*(1+ndim)+nflushmax*(twotondim))*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_INTEGER
  new_type_type(4)=MPI_DOUBLE_PRECISION
  new_type_length(1)=1
  new_type_length(2)=nflushmax*(1+ndim)
  new_type_length(3)=nflushmax*(twotondim)
  new_type_length(4)=nflushmax*(twotondim*nvar)
  call MPI_TYPE_STRUCT(4,new_type_length,new_type_disp,new_type_type,new_mpi_realdp_flush,info)
  call MPI_TYPE_COMMIT(new_mpi_realdp_flush,info)

  ! New type for request
  new_type_disp(1)=0
  new_type_disp(2)=intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_length(1)=1
  new_type_length(2)=ndim
  call MPI_TYPE_STRUCT(2,new_type_length,new_type_disp,new_type_type,new_mpi_request,info)
  call MPI_TYPE_COMMIT(new_mpi_request,info)

#endif

end subroutine init_amr



