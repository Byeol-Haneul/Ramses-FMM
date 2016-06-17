!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_basegrid
  use amr_commons
  use pm_commons
  use hilbert
  use hash, only:hash_set
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------
  ! This routine builds the coarse level grid
  !------------------------------------------
  integer::ilevel,i,j,k,igrid,ioct,ilev,info
  integer(kind=8)::ikey
  integer(kind=8),dimension(1:nvector,1:nhilbert)::hk=0
  integer(kind=8),dimension(1:nvector,1:ndim)::ix=0
  integer(kind=8),dimension(0:ndim)::hash_key,hash_test
  integer(kind=8),dimension(1:nhilbert,1:nlevelmax)::key_ref
  integer(kind=8),dimension(1:nhilbert)::coarse_key
  integer,dimension(1:nlevelmax)::n_same,npatch

  if(myid==1)write(*,*)'Building initial base grid'
  init=.true.

  ! Loop over the base level grid
  igrid=0
  do ikey=bound_key_level(1,myid-1,levelmin),bound_key_level(1,myid,levelmin)-1
     hk(1,1)=ikey
     ! Compute Cartesian index from Hilbert index
     call hilbert_reverse(ix,hk,levelmin-1,1)
     ! Insert new grid in main array
     igrid=igrid+1
     if(igrid==1)head(levelmin)=1
     tail(levelmin)=igrid
     noct(levelmin)=noct(levelmin)+1
     noct_used=noct_used+1
     grid(igrid)%lev=levelmin
     grid(igrid)%ckey(1:ndim)=ix(1,1:ndim)
     grid(igrid)%hkey(1:nhilbert)=hk(1,1:nhilbert)
     grid(igrid)%refined(1:twotondim)=.false.
     ! Insert new grid in hash table
     hash_key(0)=levelmin
     hash_key(1:ndim)=ix(1,1:ndim)
     call hash_set(grid_dict,hash_key,igrid)
  end do

  !-----------
  ! Super-octs
  !-----------
  do i=1,levelmin
     npatch(i)=twotondim**i
  end do
  ilev=levelmin
  n_same=0
  key_ref=0
  key_ref(1,1:nlevelmax)=-1
  do ioct=head(ilev),tail(ilev)
     grid(ioct)%superoct=1
     coarse_key(1:nhilbert)=grid(ioct)%hkey(1:nhilbert)
     do i=1,MIN(ilev-1,nsuperoct)
        coarse_key(1:nhilbert)=coarsen_key(coarse_key(1:nhilbert),ilev-1) ! ilev-1 used to speed up only
        if(eq_keys(coarse_key(1:nhilbert),key_ref(1:nhilbert,i)))then
           n_same(i)=n_same(i)+1
        else
           n_same(i)=1
           key_ref(1:nhilbert,i)=coarse_key(1:nhilbert)
        endif
        if(n_same(i).EQ.npatch(i))then
           grid(ioct-npatch(i)+1:ioct)%superoct=npatch(i)
        endif
     end do
  end do

  !---------------------
  ! Total number of octs
  !---------------------
  noct_tot(levelmin)=noct(levelmin)
  noct_min(levelmin)=noct(levelmin)
  noct_max(levelmin)=noct(levelmin)
  noct_used_max=noct_used
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(noct(levelmin),noct_tot(levelmin),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(noct(levelmin),noct_min(levelmin),1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(noct(levelmin),noct_max(levelmin),1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(noct_used,noct_used_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif

  if(hydro)call init_flow_fine(levelmin)
  if(pic)call rho_fine(levelmin)
  
end subroutine init_refine_basegrid
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_adaptive
  !--------------------------------------------------------------
  ! This routine builds additional refinements to the
  ! the initial AMR grid for filetype ne 'grafic'
  !--------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel,i,ivar,ilev

  if(myid==1)write(*,*)'Building initial adaptive grid'

  do i=levelmin,nlevelmax+1

     call refine_fine(levelmin)
     call load_balance(levelmin)

     do ilevel=nlevelmax,levelmin,-1
        if(hydro)then
           call init_flow_fine(ilevel)
           call upload_fine(ilevel)
        endif
     end do

     call rho_fine(levelmin)

     do ilevel=nlevelmax,levelmin,-1
        call flag_fine(ilevel,2)
     end do

  end do

  do ilevel=levelmin,nlevelmax
     call write_screen(ilevel)
  end do

  init=.false.
  
end subroutine init_refine_adaptive
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_restart
  !--------------------------------------------------------------
  ! This routine builds from a RAMSES restart file
  ! the initial AMR grid.
  !--------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use poisson_commons
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,ncpu_file,levelmin_file,nlevelmax_file
  integer::icpu,iskip_amr,iskip_hydro,iskip_grav,ilun,info
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
  integer(kind=8),dimension(1:nhilbert,1:nlevelmax)::key_ref
  integer(kind=8),dimension(1:nhilbert)::coarse_key,one_key,zero_key
  integer,dimension(1:nlevelmax)::n_same,npatch
  integer(kind=8)::ipos

  integer,dimension(1:ndim)::ckey
  logical,dimension(1:twotondim)::refined
  real(dp),dimension(1:twotondim,1:nvar)::uold
  real(dp),dimension(1:twotondim,1:ndim)::f
  real(dp),dimension(1:twotondim)::phi

  if(verbose)write(*,*)'Building adaptive grid from restart file',nrestart

  ! Read parameters from restart file
  call title(nrestart,nchar)
  file_params='output_'//TRIM(nchar)//'/params.out'
  call input_params(file_params,ncpu_file,levelmin_file,nlevelmax_file)
  if(myid==1)write(*,'(" Restart snapshot has levelmin=",I4)')levelmin_file
  if(myid==1)write(*,'(" Restart snapshot has levelmax=",I4)')nlevelmax_file

  ! Compute the proper level interval
  levelmin_max=MAX(levelmin,levelmin_file)
  nlevelmax_min=MIN(nlevelmax,nlevelmax_file)

  ! Allocate local variables
  allocate(noct_file(1:ncpu_file))
  allocate(noct_cum(1:ncpu_file))
  allocate(noct_skip(1:ncpu_file))
  allocate(ntarget_cum(1:ncpu))
  allocate(bound_key_target(1:nhilbert,0:ncpu))
  allocate(bound_key_target_tot(1:nhilbert,0:ncpu))

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
     if(myid==1)then
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
     do icpu=1,ncpu
        ntarget_cum(icpu)=int(dble(icpu)*dble(noct_cum(ncpu_file))/dble(ncpu))
     end do

     ! Compute interval of octs for current process
     if(myid>1)then
        nleft=ntarget_cum(myid-1)
     else
        nleft=0
     endif
     nright=ntarget_cum(myid)

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
        if(hydro)then
           file_hydro='output_'//TRIM(nchar)//'/hydro.out'//TRIM(ncharcpu)
           open(unit=11,file=file_hydro,access="stream"&
                & ,action="read",form='unformatted')
           iskip_hydro=17+4*(nlevelmax_file-levelmin_file+1)+&
                & (8*twotondim*nvar)*noct_skip(icpu)
        endif
        ! Prepare reading the GRAV file
        if(poisson)then
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
           if(hydro)then
              ipos=iskip_hydro+(8*twotondim*nvar)*(i-1)
              read(11,POS=ipos)uold
           endif
           if(poisson)then
              ipos=iskip_grav+(8*twotondim*(ndim+1))*(i-1)
              read(11,POS=ipos)phi
              ipos=ipos+8*twotondim
              read(11,POS=ipos)f              
           endif
           
           ! Create new oct in memory
           igrid=igrid+1
           igrid_level=igrid_level+1
           if(igrid_level==1)head(ilevel)=igrid
           tail(ilevel)=igrid
           noct(ilevel)=noct(ilevel)+1

           ! Fill values from files
           grid(igrid)%lev=ilevel
           grid(igrid)%ckey=ckey
           grid(igrid)%refined=refined
           if(hydro)then
              fluid(igrid)%uold=uold
           endif
           if(poisson)then
              grav(igrid)%phi=phi
              grav(igrid)%f=f
           endif

           ! Set flag1 to preserve refinements
           do ind=1,twotondim
              if(grid(igrid)%refined(ind))then
                 grid(igrid)%flag1(ind)=1
              else
                 grid(igrid)%flag1(ind)=0
              endif
           end do

           ! Insert in hash table
           hash_key(0)=ilevel
           hash_key(1:ndim)=ckey
           call hash_set(grid_dict,hash_key,igrid)

           ! Compute Hilbert keys of new octs
           ix(1,1:ndim)=ckey(1:ndim)
           call hilbert_key(ix,hk,dummy_state,0,ilevel-1,1)
           grid(igrid)%hkey(1:nhilbert)=hk(1,1:nhilbert)
           bound_key_target(1:nhilbert,myid)=hk(1,1:nhilbert)+one_key
        end do
        close(10)
        if(hydro)then
           close(11)
        endif
     end do

     !--------------------------------
     ! Set up new domain decomposition
     !--------------------------------
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(bound_key_target,bound_key_target_tot,nhilbert*(ncpu+1),MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
     bound_key_target=bound_key_target_tot
#endif
     bound_key_target(1:nhilbert,0)=zero_key
     do icpu=1,ncpu
        if(gt_keys(bound_key_target(1:nhilbert,icpu-1),bound_key_target(1:nhilbert,icpu)))then
           bound_key_target(1:nhilbert,icpu)=bound_key_target(1:nhilbert,icpu-1)
        endif
     end do
     bound_key_target(1:nhilbert,ncpu)=hkey_max(1:nhilbert,ilevel)
     do icpu=0,ncpu
        bound_key_level(1:nhilbert,icpu,ilevel)=bound_key_target(1:nhilbert,icpu)
     end do
     
  end do
  ! End loop over levels

  !-----------------------------
  ! Clean up final AMR structure
  !-----------------------------
  do ilevel=levelmin+1,nlevelmax
     head(ilevel)=head(ilevel-1)+noct(ilevel-1)
     tail(ilevel)=head(ilevel)+noct(ilevel)-1
  end do
  noct_used=tail(nlevelmax)

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
  do ilev=1,nlevelmax
     npatch(ilev)=twotondim**ilev
  end do
  do ilev=levelmin,nlevelmax
     n_same=0
     key_ref=0
     key_ref(1,1:nlevelmax)=-1
     do ioct=head(ilev),tail(ilev)
        grid(ioct)%superoct=1
        coarse_key(1:nhilbert)=grid(ioct)%hkey(1:nhilbert)
        do i=1,MIN(ilev-1,nsuperoct)
           coarse_key(1:nhilbert)=coarsen_key(coarse_key(1:nhilbert),ilev-1) ! ilev-1 used to speed up only
           if(eq_keys(coarse_key(1:nhilbert),key_ref(1:nhilbert,i)))then
              n_same(i)=n_same(i)+1
           else
              n_same(i)=1
              key_ref(1:nhilbert,i)=coarse_key(1:nhilbert)
           endif
           if(n_same(i).EQ.npatch(i))then
              grid(ioct-npatch(i)+1:ioct)%superoct=npatch(i)
           endif
        end do
     end do
  end do

  !---------------------
  ! Total number of octs
  !---------------------
  do ilev=levelmin,nlevelmax
     noct_tot(ilev)=noct(ilev)
     noct_min(ilev)=noct(ilev)
     noct_max(ilev)=noct(ilev)
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(noct(ilev),noct_tot(ilev),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(noct(ilev),noct_min(ilev),1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(noct(ilev),noct_max(ilev),1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif
  end do

  noct_used_max=noct_used
  noct_used_tot=noct_used
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(noct_used,noct_used_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(noct_used,noct_used_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif

end subroutine init_refine_restart
!################################################################
!################################################################
!################################################################
!################################################################
