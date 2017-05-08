subroutine init_amr(r,g,m)
  use amr_commons
  use poisson_commons
  use hash
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif

  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m

  integer::ilevel,icpu,info,igrid,i
  integer::intex,realdpex,msg_size
  integer(kind=8)::max_key
  integer,dimension(1:10)::new_type_disp,new_type_type,new_type_length,new_type_address
  real(kind=4)::real_mem,real_mem_tot
  character(len=5)::nchar
  character(len=80)::file_params
  integer::ncpu_file,levelmin_file,nlevelmax_file

  if(verbose.and.myid==1)write(*,*)'Entering init_amr'

  ! Initial time step for each level
  dtold=0.0D0
  dtnew=0.0D0
  g%dtold=0.0D0
  g%dtnew=0.0D0

  ! Allocate main oct array
  allocate(grid(1:ngridmax+ncachemax))
  do igrid=1,ngridmax+ncachemax
     grid(igrid)%lev=0
  end do
  allocate(m%grid(1:r%ngridmax+r%ncachemax))
  do igrid=1,r%ngridmax+r%ncachemax
     m%grid(igrid)%lev=0
  end do

  ! Allocate cache-related arrays
  allocate(dirty(1:ncachemax))
  allocate(locked(1:ncachemax))
  allocate(occupied(1:ncachemax))
  allocate(parent_cpu(1:ncachemax))
  dirty=.false.
  locked=.false.
  occupied=.false.
  free_cache=1; ncache=0
  allocate(m%dirty(1:r%ncachemax))
  allocate(m%locked(1:r%ncachemax))
  allocate(m%occupied(1:r%ncachemax))
  allocate(m%parent_cpu(1:r%ncachemax))
  m%dirty=.false.
  m%locked=.false.
  m%occupied=.false.
  m%free_cache=1; m%ncache=0

  allocate(lev_null(1:ncachemax))
  allocate(ckey_null(1:ndim,1:ncachemax))
  allocate(occupied_null(1:ncachemax))
  occupied_null=.false.
  free_null=1; nnull=0
  allocate(m%lev_null(1:r%ncachemax))
  allocate(m%ckey_null(1:ndim,1:r%ncachemax))
  allocate(m%occupied_null(1:r%ncachemax))
  m%occupied_null=.false.
  m%free_null=1; m%nnull=0

  ! Allocate hash table for AMR data
  if(verbose.and.myid==1)write(*,*)'Initialize empty hash'
  call init_empty_hash(grid_dict,2*(ngridmax+ncachemax),'simple')
  call init_empty_hash(m%grid_dict,2*(r%ngridmax+r%ncachemax),'simple')

  ! Allocate another smaller hash table for multigrid data
  if(poisson)then
     call init_empty_hash(mg_dict,2*(ngridmax+ncachemax)/7,'simple')
     call init_empty_hash(m%mg_dict,2*(r%ngridmax+r%ncachemax)/7,'simple')
  endif

  ! Set initial cpu boundaries
  ! Set maximum Cartesian key per level
  if(verbose.and.myid==1)write(*,*)'Initialize level cpu boundaries'

  allocate(ckey_max(1:nlevelmax+1))
  allocate(hkey_max(1:nhilbert,1:nlevelmax+1))
  hkey_max=0
  allocate(m%ckey_max(1:r%nlevelmax+1))
  allocate(m%hkey_max(1:nhilbert,1:r%nlevelmax+1))
  m%hkey_max=0

  allocate(domain(1:nlevelmax+1))
  do ilevel=1,nlevelmax+1
     call domain(ilevel)%create(ncpu*overload)
  end do
  allocate(m%domain(1:r%nlevelmax+1))
  do ilevel=1,r%nlevelmax+1
     call m%domain(ilevel)%create(g%ncpu*r%overload)
  end do

  ! Make sure that the coarsest level uses only one Hilbert integer
  if(levelmin>(levels_per_key(ndim)+1))then
     write(*,*)'levelmin is way too large !'
     write(*,*)'are you crazy ?'
     stop
  endif

  ! Set bounds for Hilbert keys for coarse levels
  do ilevel=1,levelmin
     ckey_max(ilevel)=2**(ilevel-1)
     max_key=2**((ilevel-1)*ndim)
     hkey_max(1,ilevel)=max_key
     do icpu=1,ncpu-1
        domain(ilevel)%b(1,icpu) = (icpu*max_key)/ncpu
     end do
     domain(ilevel)%b(1,0) = 0
     domain(ilevel)%b(1,ncpu) = max_key
  end do
  do ilevel=1,r%levelmin
     m%ckey_max(ilevel)=2**(ilevel-1)
     max_key=2**((ilevel-1)*ndim)
     m%hkey_max(1,ilevel)=max_key
     do icpu=1,g%ncpu-1
        m%domain(ilevel)%b(1,icpu) = (icpu*max_key)/g%ncpu
     end do
     m%domain(ilevel)%b(1,0) = 0
     m%domain(ilevel)%b(1,ncpu) = max_key
  end do

  ! Set bounds for Hilbert keys for fine levels
  do ilevel=levelmin+1,nlevelmax+1
     ckey_max(ilevel) = 2**(ilevel-1)
     hkey_max(1:nhilbert,ilevel) = refine_key(hkey_max(1:nhilbert,ilevel-1),ilevel-2)
     ! Multiply the bounds by twotondim
     do icpu=0,ncpu
        domain(ilevel)%b(1:nhilbert,icpu) = refine_key(domain(ilevel-1)%b(1:nhilbert,icpu),ilevel-2)
     end do
  end do
  do ilevel=r%levelmin+1,r%nlevelmax+1
     m%ckey_max(ilevel) = 2**(ilevel-1)
     m%hkey_max(1:nhilbert,ilevel) = refine_key(m%hkey_max(1:nhilbert,ilevel-1),ilevel-2)
     ! Multiply the bounds by twotondim
     do icpu=0,g%ncpu
        m%domain(ilevel)%b(1:nhilbert,icpu) = refine_key(m%domain(ilevel-1)%b(1:nhilbert,icpu),ilevel-2)
     end do
  end do

  ! Allocate head, tail and numbers for each level
  if(verbose.and.myid==1)write(*,*)'Initialize oct decomposition'
  allocate(head(levelmin:nlevelmax))
  allocate(tail(levelmin:nlevelmax))
  allocate(head_cache(1:nlevelmax))
  allocate(tail_cache(1:nlevelmax))
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
  allocate(m%head(r%levelmin:r%nlevelmax))
  allocate(m%tail(r%levelmin:r%nlevelmax))
  allocate(m%head_cache(1:r%nlevelmax))
  allocate(m%tail_cache(1:r%nlevelmax))
  allocate(m%noct(r%levelmin:r%nlevelmax))
  allocate(m%noct_min(r%levelmin:r%nlevelmax))
  allocate(m%noct_max(r%levelmin:r%nlevelmax))
  allocate(m%noct_tot(r%levelmin:r%nlevelmax))
  m%head=1       ! Head oct in the level
  m%tail=0       ! Tail oct in the level
  m%noct=0       ! Number of oct in the level and in the cpu
  m%noct_tot=0   ! Total number of oct in the level (all cpus)
  m%noct_min=0   ! Minimum number of oct across all cpus
  m%noct_max=0   ! Maximum number of oct across all cpus
  m%noct_used=0  ! Number of oct used across all levels
  m%noct_used_tot=0  ! Total number of oct used (all cpus)

  if(nrestart>0)then
     ! Read parameters from restart file
     call title(nrestart,nchar)
     file_params='output_'//TRIM(nchar)//'/params.out'
     call input_params(file_params,ncpu_file,levelmin_file,nlevelmax_file)
     if(myid==1)write(*,'(" Restarting from output number ",I8)')nrestart
     if(myid==1)write(*,'(" Restart snapshot has ",I8," files")')ncpu_file
  endif

#ifndef WITHOUTMPI  
  ! Allocate all communication and cache-related variables
  allocate(reply_id(1:ncpu))
  allocate(reply_interpol(1:ncpu))
  allocate(reply_mg(1:ncpu))
  allocate(reply_flag(1:ncpu))
  allocate(reply_hydro(1:ncpu))
  allocate(reply_poisson(1:ncpu))
  allocate(reply_refine(1:ncpu))
  allocate(send_flush_interpol(1:ncpu))
  allocate(send_flush_mg(1:ncpu))
  allocate(send_flush_flag(1:ncpu))
  allocate(send_flush_hydro(1:ncpu))
  allocate(send_flush_poisson(1:ncpu))
  allocate(send_flush_refine(1:ncpu))
  do icpu=1,ncpu
     send_flush_interpol(icpu)%nflush=0
     send_flush_mg(icpu)%nflush=0
     send_flush_flag(icpu)%nflush=0
     send_flush_hydro(icpu)%nflush=0
     send_flush_poisson(icpu)%nflush=0
     send_flush_refine(icpu)%nflush=0
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

  ! New type for small_realdp_msg
  new_type_disp(1)=0
  new_type_disp(2)=2*intex
  new_type_disp(3)=(2+ntilemax*(1+ndim))*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_DOUBLE_PRECISION
  new_type_length(1)=2
  new_type_length(2)=ntilemax*(1+ndim)
  new_type_length(3)=ntilemax*(twotondim)
  call MPI_TYPE_STRUCT(3,new_type_length,new_type_disp,new_type_type,new_mpi_small_realdp_msg,info)
  call MPI_TYPE_COMMIT(new_mpi_small_realdp_msg,info)

  ! New type for twin_realdp_msg
  new_type_disp(1)=0
  new_type_disp(2)=2*intex
  new_type_disp(3)=(2+ntilemax*(1+ndim))*intex
  new_type_disp(4)=(2+ntilemax*(1+ndim))*intex+ntilemax*twotondim*realdpex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_DOUBLE_PRECISION
  new_type_type(4)=MPI_DOUBLE_PRECISION
  new_type_length(1)=2
  new_type_length(2)=ntilemax*(1+ndim)
  new_type_length(3)=ntilemax*(twotondim)
  new_type_length(4)=ntilemax*(twotondim)
  call MPI_TYPE_STRUCT(4,new_type_length,new_type_disp,new_type_type,new_mpi_twin_realdp_msg,info)
  call MPI_TYPE_COMMIT(new_mpi_twin_realdp_msg,info)

  ! New type for three_realdp_msg
  new_type_disp(1)=0
  new_type_disp(2)=2*intex
  new_type_disp(3)=(2+ntilemax*(1+ndim))*intex
  new_type_disp(4)=(2+ntilemax*(1+ndim))*intex+ntilemax*twotondim*realdpex
  new_type_disp(5)=(2+ntilemax*(1+ndim))*intex+2*ntilemax*twotondim*realdpex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_DOUBLE_PRECISION
  new_type_type(4)=MPI_DOUBLE_PRECISION
  new_type_type(5)=MPI_DOUBLE_PRECISION
  new_type_length(1)=2
  new_type_length(2)=ntilemax*(1+ndim)
  new_type_length(3)=ntilemax*(twotondim)
  new_type_length(4)=ntilemax*(twotondim)
  new_type_length(5)=ntilemax*(twotondim)
  call MPI_TYPE_STRUCT(5,new_type_length,new_type_disp,new_type_type,new_mpi_three_realdp_msg,info)
  call MPI_TYPE_COMMIT(new_mpi_three_realdp_msg,info)

  ! New type for large_realdp_msg
  new_type_disp(1)=0
  new_type_disp(2)=2*intex
  new_type_disp(3)=new_type_disp(2)+ntilemax*(1+ndim)*intex
  new_type_disp(4)=new_type_disp(3)+ntilemax*twotondim*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_INTEGER
  new_type_length(1)=2
  new_type_length(2)=ntilemax*(1+ndim)
  new_type_length(3)=ntilemax*twotondim
  msg_size=3

#ifdef HYDRO
  msg_size=msg_size+1
  new_type_length(msg_size)=ntilemax*twotondim*nvar
  new_type_type(msg_size)=MPI_DOUBLE_PRECISION
  new_type_disp(msg_size+1)=new_type_disp(msg_size)+ntilemax*twotondim*nvar*realdpex
#endif

#ifdef GRAV
  msg_size=msg_size+1
  new_type_length(msg_size)=ntilemax*twotondim*(ndim+2)
  new_type_type(msg_size)=MPI_DOUBLE_PRECISION
  new_type_disp(msg_size+1)=new_type_disp(msg_size)+ntilemax*twotondim*(ndim+2)*realdpex
#endif

  call MPI_TYPE_STRUCT(msg_size,new_type_length,new_type_disp,new_type_type,new_mpi_large_realdp_msg,info)
  call MPI_TYPE_COMMIT(new_mpi_large_realdp_msg,info)

  ! New type for small_realdp_flush
  new_type_disp(1)=0
  new_type_disp(2)=intex
  new_type_disp(3)=(1+nflushmax*(1+ndim))*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_DOUBLE_PRECISION
  new_type_length(1)=1
  new_type_length(2)=nflushmax*(1+ndim)
  new_type_length(3)=nflushmax*twotondim
  call MPI_TYPE_STRUCT(3,new_type_length,new_type_disp,new_type_type,new_mpi_small_realdp_flush,info)
  call MPI_TYPE_COMMIT(new_mpi_small_realdp_flush,info)

  ! New type for twin_realdp_flush
  new_type_disp(1)=0
  new_type_disp(2)=intex
  new_type_disp(3)=(1+nflushmax*(1+ndim))*intex
  new_type_disp(4)=(1+nflushmax*(1+ndim))*intex+nflushmax*twotondim*realdpex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_DOUBLE_PRECISION
  new_type_type(4)=MPI_DOUBLE_PRECISION
  new_type_length(1)=1
  new_type_length(2)=nflushmax*(1+ndim)
  new_type_length(3)=nflushmax*twotondim
  new_type_length(4)=nflushmax*twotondim
  call MPI_TYPE_STRUCT(4,new_type_length,new_type_disp,new_type_type,new_mpi_twin_realdp_flush,info)
  call MPI_TYPE_COMMIT(new_mpi_twin_realdp_flush,info)

  ! New type for three_realdp_flush
  new_type_disp(1)=0
  new_type_disp(2)=intex
  new_type_disp(3)=(1+nflushmax*(1+ndim))*intex
  new_type_disp(4)=(1+nflushmax*(1+ndim))*intex+nflushmax*twotondim*realdpex
  new_type_disp(5)=(1+nflushmax*(1+ndim))*intex+2*nflushmax*twotondim*realdpex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_DOUBLE_PRECISION
  new_type_type(4)=MPI_DOUBLE_PRECISION
  new_type_type(5)=MPI_DOUBLE_PRECISION
  new_type_length(1)=1
  new_type_length(2)=nflushmax*(1+ndim)
  new_type_length(3)=nflushmax*twotondim
  new_type_length(4)=nflushmax*twotondim
  new_type_length(5)=nflushmax*twotondim
  call MPI_TYPE_STRUCT(5,new_type_length,new_type_disp,new_type_type,new_mpi_three_realdp_flush,info)
  call MPI_TYPE_COMMIT(new_mpi_three_realdp_flush,info)

  ! New type for large_realdp_flush
  new_type_disp(1)=0
  new_type_disp(2)=intex
  new_type_disp(3)=new_type_disp(2)+nflushmax*(1+ndim)*intex
  new_type_disp(4)=new_type_disp(3)+nflushmax*twotondim*intex
  new_type_type(1)=MPI_INTEGER
  new_type_type(2)=MPI_INTEGER
  new_type_type(3)=MPI_INTEGER
  new_type_length(1)=1
  new_type_length(2)=nflushmax*(1+ndim)
  new_type_length(3)=nflushmax*twotondim
  msg_size=3

#ifdef HYDRO
  msg_size=msg_size+1
  new_type_type(msg_size)=MPI_DOUBLE_PRECISION
  new_type_length(msg_size)=nflushmax*twotondim*nvar
  new_type_disp(msg_size+1)=new_type_disp(msg_size)+nflushmax*twotondim*nvar*realdpex
#endif

#ifdef GRAV
  msg_size=msg_size+1
  new_type_type(msg_size)=MPI_DOUBLE_PRECISION
  new_type_length(msg_size)=nflushmax*twotondim*(ndim+2)
  new_type_disp(msg_size+1)=new_type_disp(msg_size)+nflushmax*twotondim*(ndim+2)*realdpex
#endif

  call MPI_TYPE_STRUCT(msg_size,new_type_length,new_type_disp,new_type_type,new_mpi_large_realdp_flush,info)
  call MPI_TYPE_COMMIT(new_mpi_large_realdp_flush,info)

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



