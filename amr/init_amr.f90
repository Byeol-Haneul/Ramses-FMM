!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_set_add(pst,input_size,output_size,iUpper)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer::iUpper

  integer::n,iLower,iMiddle

  associate(mdl=>pst%s%mdl)

 iLower = mdl%myid-1
 n = iUpper - iLower
 iMiddle = (iUpper + iLower) / 2
 
 if(n>1)then
    pst%iUpper = iMiddle
    pst%nLower = iMiddle - iLower
    pst%nUpper = iUpper - iMiddle
    allocate(pst%pLower)
    pst%pLower%s => pst%s
    call mdl_send_request(mdl,MDL_SET_ADD,pst%iUpper+1,input_size,output_size,iUpper)
    call r_set_add(pst%pLower,input_size,output_size,iMiddle)
    call mdl_get_reply(mdl,pst%iUpper+1,output_size)
 end if
 
 end associate

end subroutine r_set_add
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_init_amr(pst,input_size,output_size)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_INIT_AMR,pst%iUpper+1,input_size,output_size)
     call r_init_amr(pst%pLower,input_size,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     call init_amr(pst%s%r,pst%s%g,pst%s%m)
  endif

end subroutine r_init_amr
!###############################################
!###############################################
!###############################################
!###############################################
subroutine init_amr(r,g,m)
  use amr_parameters, ONLY: nhilbert
  use amr_commons, ONLY: run_t, global_t, mesh_t
  use hash
  use hilbert
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m

  ! Local variables
  integer::ilevel,icpu,igrid
  integer(kind=8)::max_key
  character(len=5)::nchar
  character(len=80)::file_params
  integer::ncpu_file,levelmin_file,nlevelmax_file

  ! Initial time step for each level
  g%dtold=0.0D0
  g%dtnew=0.0D0

  ! Allocate main oct array
  allocate(m%grid(1:r%ngridmax+r%ncachemax))
  do igrid=1,r%ngridmax+r%ncachemax
     m%grid(igrid)%lev=0
  end do

  ! Allocate cache-related arrays
  allocate(m%dirty(1:r%ncachemax))
  allocate(m%locked(1:r%ncachemax))
  allocate(m%occupied(1:r%ncachemax))
  allocate(m%parent_cpu(1:r%ncachemax))
  m%dirty=.false.
  m%locked=.false.
  m%occupied=.false.
  m%free_cache=1; m%ncache=0

  allocate(m%lev_null(1:r%ncachemax))
  allocate(m%ckey_null(1:ndim,1:r%ncachemax))
  allocate(m%occupied_null(1:r%ncachemax))
  m%occupied_null=.false.
  m%free_null=1; m%nnull=0

  ! Allocate hash table for AMR data
  if(r%verbose.and.g%myid==1)write(*,*)'Initialize empty hash'
  call init_empty_hash(m%grid_dict,2*(r%ngridmax+r%ncachemax),'simple')

  ! Allocate another smaller hash table for multigrid data
  if(r%poisson)then
     call init_empty_hash(m%mg_dict,2*(r%ngridmax+r%ncachemax)/7,'simple')
  endif

  ! Set initial cpu boundaries
  ! Set maximum Cartesian key per level
  if(r%verbose.and.g%myid==1)write(*,*)'Initialize level cpu boundaries'

  allocate(m%ckey_max(1:r%nlevelmax+1))
  allocate(m%hkey_max(1:nhilbert,1:r%nlevelmax+1))
  m%hkey_max=0

  allocate(m%domain(1:r%nlevelmax+1))
  do ilevel=1,r%nlevelmax+1
     call m%domain(ilevel)%create(g%myid,g%ncpu,g%ncpu)
  end do

  ! Make sure that the coarsest level uses only one Hilbert integer
  if(r%levelmin>(levels_per_key(ndim)+1))then
     write(*,*)'levelmin is way too large !'
     write(*,*)'are you crazy ?'
     stop
  endif

  ! Set bounds for Hilbert keys for coarse levels
  do ilevel=1,r%levelmin
     m%ckey_max(ilevel)=2**(ilevel-1)
     max_key=2**((ilevel-1)*ndim)
     m%hkey_max(1,ilevel)=max_key
     do icpu=1,g%ncpu-1
        m%domain(ilevel)%b(1,icpu) = (icpu*max_key)/g%ncpu
     end do
     m%domain(ilevel)%b(1,0) = 0
     m%domain(ilevel)%b(1,g%ncpu) = max_key
  end do

  ! Set bounds for Hilbert keys for fine levels
  do ilevel=r%levelmin+1,r%nlevelmax+1
     m%ckey_max(ilevel) = 2**(ilevel-1)
     m%hkey_max(1:nhilbert,ilevel) = refine_key(m%hkey_max(1:nhilbert,ilevel-1),ilevel-2)
     ! Multiply the bounds by twotondim
     do icpu=0,g%ncpu
        m%domain(ilevel)%b(1:nhilbert,icpu) = refine_key(m%domain(ilevel-1)%b(1:nhilbert,icpu),ilevel-2)
     end do
  end do

  ! Allocate head, tail and numbers for each level
  if(r%verbose.and.g%myid==1)write(*,*)'Initialize oct decomposition'
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

  if(r%nrestart>0)then
     ! Read parameters from restart file
     call title(r%nrestart,nchar)
     file_params='output_'//TRIM(nchar)//'/params.out'
     call input_params(r,g,file_params,ncpu_file,levelmin_file,nlevelmax_file)
     if(g%myid==1)write(*,'(" Restarting from output number ",I8)')r%nrestart
     if(g%myid==1)write(*,'(" Restart snapshot has ",I8," files")')ncpu_file
  endif

end subroutine init_amr
!###############################################
!###############################################
!###############################################
!###############################################



