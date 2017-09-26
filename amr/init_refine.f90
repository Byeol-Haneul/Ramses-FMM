!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_basegrid_2(r,g,m,p)
  use amr_parameters, only: nvector,nhilbert,ndim,twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use hilbert
  use hash, only: hash_set
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  !------------------------------------------
  ! This routine builds the coarse level grid
  !------------------------------------------
  integer::i,igrid,ioct,ilev
  integer(kind=8)::ikey
  integer(kind=8),dimension(1:nvector,1:nhilbert)::hk=0
  integer(kind=8),dimension(1:nvector,1:ndim)::ix=0
  integer(kind=8),dimension(0:ndim)::hash_key
  integer(kind=8),dimension(1:nhilbert,1:r%nlevelmax)::key_ref
  integer(kind=8),dimension(1:nhilbert)::coarse_key
  integer,dimension(1:r%nlevelmax)::n_same,npatch

  if(g%myid==1)write(*,*)'Building initial base grid'
  g%init=.true.

  ! Loop over the base level grid
  igrid=0
  do ikey=m%domain(r%levelmin)%b(1,g%myid-1), m%domain(r%levelmin)%b(1,g%myid)-1
     hk(1,1)=ikey
     ! Compute Cartesian index from Hilbert index
     call hilbert_reverse(ix,hk,r%levelmin-1,1)
     ! Insert new grid in main array
     igrid=igrid+1
     if(igrid==1)m%head(r%levelmin)=1
     m%tail(r%levelmin)=igrid
     m%noct(r%levelmin)=m%noct(r%levelmin)+1
     m%noct_used=m%noct_used+1
     m%grid(igrid)%lev=r%levelmin
     m%grid(igrid)%ckey(1:ndim)=int(ix(1,1:ndim),kind=4)
     m%grid(igrid)%hkey(1:nhilbert)=hk(1,1:nhilbert)
     m%grid(igrid)%refined(1:twotondim)=.false.
     ! Insert new grid in hash table
     hash_key(0)=r%levelmin
     hash_key(1:ndim)=ix(1,1:ndim)
     call hash_set(m%grid_dict,hash_key,igrid)
  end do

  !-----------
  ! Super-octs
  !-----------
  do i=1,r%levelmin
     npatch(i)=twotondim**i
  end do
  ilev=r%levelmin
  n_same=0
  key_ref=0
  key_ref(1,1:r%nlevelmax)=-1
  do ioct=m%head(ilev),m%tail(ilev)
     m%grid(ioct)%superoct=1
     coarse_key(1:nhilbert)=m%grid(ioct)%hkey(1:nhilbert)
     do i=1,MIN(ilev-1,r%nsuperoct)
        coarse_key(1:nhilbert)=coarsen_key(coarse_key(1:nhilbert),ilev-1) ! ilev-1 used to speed up only
        if(eq_keys(coarse_key(1:nhilbert),key_ref(1:nhilbert,i)))then
           n_same(i)=n_same(i)+1
        else
           n_same(i)=1
           key_ref(1:nhilbert,i)=coarse_key(1:nhilbert)
        endif
        if(n_same(i).EQ.npatch(i))then
           m%grid(ioct-npatch(i)+1:ioct)%superoct=npatch(i)
        endif
     end do
  end do

  !---------------------
  ! Total number of octs
  !---------------------
  m%noct_tot(r%levelmin)=m%noct(r%levelmin)
  m%noct_min(r%levelmin)=m%noct(r%levelmin)
  m%noct_max(r%levelmin)=m%noct(r%levelmin)
  m%noct_used_max=m%noct_used
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(m%noct(r%levelmin),m%noct_tot(r%levelmin),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(m%noct(r%levelmin),m%noct_min(r%levelmin),1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(m%noct(r%levelmin),m%noct_max(r%levelmin),1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(m%noct_used,m%noct_used_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif

  if(r%hydro)call init_flow_fine_2(r,g,m,r%levelmin)
  if(r%pic)call rho_fine_2(r,g,m,p,r%levelmin)
  
end subroutine init_refine_basegrid_2
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_adaptive_2(r,g,m,p)
  !--------------------------------------------------------------
  ! This routine builds additional refinements to the
  ! the initial AMR grid for filetype ne 'grafic'
  !--------------------------------------------------------------
  use amr_parameters, only: dp
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p

  ! Local variables
  integer::ilevel,i

  if(g%myid==1)write(*,*)'Building initial adaptive grid thread safe'

  do i=r%levelmin,r%nlevelmax+1

     call refine_fine_2(r,g,m,r%levelmin)
#ifndef WITHOUTMPI
     call load_balance_2(r,g,m,r%levelmin)
#endif

     do ilevel=r%nlevelmax,r%levelmin,-1
        if(r%hydro)then
           call init_flow_fine_2(r,g,m,ilevel)
           call upload_fine_2(r,g,m,ilevel)
        endif
     end do

     call rho_fine_2(r,g,m,p,r%levelmin)

     do ilevel=r%nlevelmax,r%levelmin,-1
        call flag_fine_2(r,g,m,ilevel,2)
     end do

  end do

  do ilevel=r%levelmin,r%nlevelmax
     call write_screen_2(m,ilevel)
  end do

  g%init=.false.
  
end subroutine init_refine_adaptive_2
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_restart_2(r,g,m)
  !--------------------------------------------------------------
  ! This routine builds from a RAMSES restart file
  ! the initial AMR grid.
  !--------------------------------------------------------------
  use amr_parameters, only: dp,nhilbert,ndim,twotondim,nvector
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  use hash
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m

  ! Local variables
  integer::ilevel,ncpu_file,levelmin_file,nlevelmax_file
  integer::icpu,iskip_amr=0,iskip_hydro=0,iskip_grav=0,ilun
  integer::i,ind,istart,iend,noct_tmp,ilev,ioct
  integer::igrid,igrid_level,nleft,nright,ileft,iright
  integer::levelmin_max,nlevelmax_min
  character(LEN=80)::file_params,file_amr,file_hydro,file_grav
  character(LEN=5)::nchar,ncharcpu

  integer,dimension(:),allocatable::noct_file,noct_skip,ntarget_cum,noct_cum
  integer(kind=8),allocatable,dimension(:,:)::bound_key_target
  integer(kind=8),allocatable,dimension(:,:)::bound_key_target_tot

  integer(kind=4), dimension(1:nvector),save::dummy_state
  integer(kind=8), dimension(1:nvector,1:nhilbert),save::hk
  integer(kind=8), dimension(1:nvector,1:ndim),save::ix
  integer(kind=8), dimension(0:ndim)::hash_key
  integer(kind=8),dimension(1:nhilbert,1:r%nlevelmax)::key_ref
  integer(kind=8),dimension(1:nhilbert)::coarse_key,one_key,zero_key
  integer,dimension(1:r%nlevelmax)::n_same,npatch
  integer(kind=8)::ipos

  integer,dimension(1:ndim)::ckey
  logical,dimension(1:twotondim)::refined
  real(dp),dimension(1:twotondim,1:nvar)::uold
  real(dp),dimension(1:twotondim,1:ndim)::f
  real(dp),dimension(1:twotondim)::phi

  if(r%verbose)write(*,*)'Building adaptive grid from restart file',r%nrestart

  ! Read parameters from restart file
  call title(r%nrestart,nchar)
  file_params='output_'//TRIM(nchar)//'/params.out'
  call input_params_2(r,g,file_params,ncpu_file,levelmin_file,nlevelmax_file)
  if(g%myid==1)write(*,'(" Restart snapshot has levelmin=",I4)')levelmin_file
  if(g%myid==1)write(*,'(" Restart snapshot has levelmax=",I4)')nlevelmax_file

  ! Compute the proper level interval
  levelmin_max=MAX(r%levelmin,levelmin_file)
  nlevelmax_min=MIN(r%nlevelmax,nlevelmax_file)

  ! Allocate local variables
  allocate(noct_file(1:ncpu_file))
  allocate(noct_cum(1:ncpu_file))
  allocate(noct_skip(1:ncpu_file))
  allocate(ntarget_cum(1:g%ncpu))
  allocate(bound_key_target(1:nhilbert,0:g%ncpu))
  allocate(bound_key_target_tot(1:nhilbert,0:g%ncpu))

  ! Set some constants
  zero_key=0
  one_key=0
  one_key(1)=1

  !--------------------------
  ! Loop over relevant levels
  !--------------------------
  igrid=0
  do ilevel=levelmin_max,nlevelmax_min
     igrid_level=0

     !------------------------------
     ! Count number of octs in files
     !------------------------------
     if(g%myid==1)then
        noct_cum=0
        do icpu=1,ncpu_file
           call title(icpu,ncharcpu)
           file_amr='output_'//TRIM(nchar)//'/amr.out'//TRIM(ncharcpu)
           ilun=10
           noct_skip(icpu)=0
           open(unit=ilun,file=file_amr,access="stream"&
                & ,action="read",form='unformatted')
           do i=levelmin_file,ilevel-1
              ipos=13+4*(i-levelmin_file)
              read(ilun,POS=ipos)noct_tmp
              noct_skip(icpu)=noct_skip(icpu)+noct_tmp
           end do
           ipos=13+4*(ilevel-levelmin_file)
           read(ilun,POS=ipos)noct_file(icpu)
           if(icpu>1)then
              noct_cum(icpu)=noct_cum(icpu-1)+noct_file(icpu)
           else
              noct_cum(icpu)=noct_file(icpu)
           endif
           close(ilun)
        end do
     endif
#ifndef WITHOUTMPI
     call MPI_BCAST(noct_cum,ncpu_file,MPI_INTEGER,0,MPI_COMM_WORLD,info)
     call MPI_BCAST(noct_file,ncpu_file,MPI_INTEGER,0,MPI_COMM_WORLD,info)
     call MPI_BCAST(noct_skip,ncpu_file,MPI_INTEGER,0,MPI_COMM_WORLD,info)
#endif

     !----------------------------
     ! New numbers of octs per cpu
     !----------------------------
     do icpu=1,g%ncpu
        ntarget_cum(icpu)=int(dble(icpu)*dble(noct_cum(ncpu_file))/dble(g%ncpu))
     end do

     ! Compute interval of octs for current process
     if(g%myid>1)then
        nleft=ntarget_cum(g%myid-1)
     else
        nleft=0
     endif
     nright=ntarget_cum(g%myid)

     !-----------------------------------------------------
     ! Compute interval of file to open for current process
     !-----------------------------------------------------
     ileft=0
     iright=-1
     bound_key_target=0
     if(nright.GT.nleft)then
        do icpu=1,ncpu_file
           if(icpu>1)then
              if(noct_cum(icpu).GT.nleft.AND.noct_cum(icpu-1).LT.nright)then
                 if(ileft==0)ileft=icpu
                 iright=MAX(icpu,iright)
              endif
           else
              if(noct_cum(icpu).GT.nleft)then
                 if(ileft==0)ileft=icpu
                 iright=MAX(icpu,iright)
              endif
           endif
        end do
     endif

     !----------------------------
     ! Read octs data in files
     !----------------------------
     ! Loop over relevant files (if any)
     do icpu=ileft,iright
        if(icpu>1)then
           istart=MAX(nleft-noct_cum(icpu-1),0)+1
           iend=MIN(nright-noct_cum(icpu-1),noct_file(icpu))
        else
           istart=nleft+1
           iend=MIN(nright,noct_file(icpu))
        endif
        call title(icpu,ncharcpu)
        ! Prepare reading the AMR file
        file_amr='output_'//TRIM(nchar)//'/amr.out'//TRIM(ncharcpu)
        open(unit=10,file=file_amr,access="stream"&
             & ,action="read",form='unformatted')
        iskip_amr=13+4*(nlevelmax_file-levelmin_file+1)+&
             & (4*ndim+4*twotondim)*noct_skip(icpu)
        ! Prepare reading the HYDRO file
        if(r%hydro)then
           file_hydro='output_'//TRIM(nchar)//'/hydro.out'//TRIM(ncharcpu)
           open(unit=11,file=file_hydro,access="stream"&
                & ,action="read",form='unformatted')
           iskip_hydro=17+4*(nlevelmax_file-levelmin_file+1)+&
                & (8*twotondim*nvar)*noct_skip(icpu)
        endif
        ! Prepare reading the GRAV file
        if(r%poisson)then
           file_grav='output_'//TRIM(nchar)//'/grav.out'//TRIM(ncharcpu)
           open(unit=11,file=file_grav,access="stream"&
                & ,action="read",form='unformatted')
           iskip_grav=17+4*(nlevelmax_file-levelmin_file+1)+&
                & (8*twotondim*(ndim+1))*noct_skip(icpu)
        endif
        ! Loop over useful octs in file
        do i=istart,iend

           ! Read values from files
           ipos=iskip_amr+(4*ndim+4*twotondim)*(i-1)
           read(10,POS=ipos)ckey
           ipos=iskip_amr+(4*ndim+4*twotondim)*(i-1)+4*ndim
           read(10,POS=ipos)refined
           if(r%hydro)then
              ipos=iskip_hydro+(8*twotondim*nvar)*(i-1)
              read(11,POS=ipos)uold
           endif
           if(r%poisson)then
              ipos=iskip_grav+(8*twotondim*(ndim+1))*(i-1)
              read(11,POS=ipos)phi
              ipos=ipos+8*twotondim
              read(11,POS=ipos)f              
           endif
           
           ! Create new oct in memory
           igrid=igrid+1
           igrid_level=igrid_level+1
           if(igrid_level==1)m%head(ilevel)=igrid
           m%tail(ilevel)=igrid
           m%noct(ilevel)=m%noct(ilevel)+1

           ! Fill values from files
           m%grid(igrid)%lev=ilevel
           m%grid(igrid)%ckey=ckey
           m%grid(igrid)%refined=refined
#ifdef HYDRO
           if(r%hydro)then
              m%grid(igrid)%uold=uold
           endif
#endif
#ifdef GRAV
           if(r%poisson)then
              m%grid(igrid)%phi=phi
              m%grid(igrid)%f=f
           endif
#endif

           ! Set flag1 to preserve refinements
           do ind=1,twotondim
              if(m%grid(igrid)%refined(ind))then
                 m%grid(igrid)%flag1(ind)=1
              else
                 m%grid(igrid)%flag1(ind)=0
              endif
           end do

           ! Insert in hash table
           hash_key(0)=ilevel
           hash_key(1:ndim)=ckey
           call hash_set(m%grid_dict,hash_key,igrid)

           ! Compute Hilbert keys of new octs
           ix(1,1:ndim)=ckey(1:ndim)
           call hilbert_key(ix,hk,dummy_state,0,ilevel-1,1)
           m%grid(igrid)%hkey(1:nhilbert)=hk(1,1:nhilbert)
           bound_key_target(1:nhilbert,g%myid)=hk(1,1:nhilbert)+one_key
        end do
        close(10)
        if(r%hydro)then
           close(11)
        endif
     end do

     !--------------------------------
     ! Set up new domain decomposition
     !--------------------------------
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(bound_key_target,bound_key_target_tot,nhilbert*(g%ncpu+1),MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
     bound_key_target=bound_key_target_tot
#endif
     bound_key_target(1:nhilbert,0)=zero_key
     do icpu=1,g%ncpu
        if(gt_keys(bound_key_target(1:nhilbert,icpu-1),bound_key_target(1:nhilbert,icpu)))then
           bound_key_target(1:nhilbert,icpu)=bound_key_target(1:nhilbert,icpu-1)
        endif
     end do
     bound_key_target(1:nhilbert,g%ncpu)=m%hkey_max(1:nhilbert,ilevel)
     do icpu=0,g%ncpu
        m%domain(ilevel)%b(1:nhilbert,icpu)=bound_key_target(1:nhilbert,icpu)
     end do
     
  end do
  ! End loop over levels

  !-----------------------------
  ! Clean up final AMR structure
  !-----------------------------
  do ilevel=r%levelmin+1,r%nlevelmax
     m%head(ilevel)=m%head(ilevel-1)+m%noct(ilevel-1)
     m%tail(ilevel)=m%head(ilevel)+m%noct(ilevel)-1
  end do
  m%noct_used=m%tail(r%nlevelmax)

  ! Deallocate local arrays
  deallocate(noct_file)
  deallocate(noct_cum)
  deallocate(noct_skip)
  deallocate(ntarget_cum)
  deallocate(bound_key_target)
  deallocate(bound_key_target_tot)

  !-----------
  ! Super-octs
  !-----------
  do ilev=1,r%nlevelmax
     npatch(ilev)=twotondim**ilev
  end do
  do ilev=r%levelmin,r%nlevelmax
     n_same=0
     key_ref=0
     key_ref(1,1:r%nlevelmax)=-1
     do ioct=m%head(ilev),m%tail(ilev)
        m%grid(ioct)%superoct=1
        coarse_key(1:nhilbert)=m%grid(ioct)%hkey(1:nhilbert)
        do i=1,MIN(ilev-1,r%nsuperoct)
           coarse_key(1:nhilbert)=coarsen_key(coarse_key(1:nhilbert),ilev-1) ! ilev-1 used to speed up only
           if(eq_keys(coarse_key(1:nhilbert),key_ref(1:nhilbert,i)))then
              n_same(i)=n_same(i)+1
           else
              n_same(i)=1
              key_ref(1:nhilbert,i)=coarse_key(1:nhilbert)
           endif
           if(n_same(i).EQ.npatch(i))then
              m%grid(ioct-npatch(i)+1:ioct)%superoct=npatch(i)
           endif
        end do
     end do
  end do

  !---------------------
  ! Total number of octs
  !---------------------
  do ilev=r%levelmin,r%nlevelmax
     m%noct_tot(ilev)=m%noct(ilev)
     m%noct_min(ilev)=m%noct(ilev)
     m%noct_max(ilev)=m%noct(ilev)
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(m%noct(ilev),m%noct_tot(ilev),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(m%noct(ilev),m%noct_min(ilev),1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(m%noct(ilev),m%noct_max(ilev),1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif
  end do

  m%noct_used_max=m%noct_used
  m%noct_used_tot=m%noct_used
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(m%noct_used,m%noct_used_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(m%noct_used,m%noct_used_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif

end subroutine init_refine_restart_2
!################################################################
!################################################################
!################################################################
!################################################################
