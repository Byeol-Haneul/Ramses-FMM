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
  integer::ilevel,i,j,k,igrid,ipos,ioct,ilev,info
  integer(kind=8)::ikey
  integer(kind=8),dimension(1:nvector)::hk0,hk1,hk2
  integer(kind=8),dimension(1:nvector)::ix=0,iy=0,iz=0
  integer(kind=8),dimension(0:ndim)::hash_key,hash_test
  integer(kind=8),dimension(1:nlevelmax)::key_ref
  integer(kind=8)::coarse_key
  integer,dimension(1:nlevelmax)::n_same,npatch

  if(myid==1)write(*,*)'Building initial base grid'
  init=.true.

  ! Loop over the base level grid
  igrid=0
  do ikey=bound_key_level(myid-1,levelmin),bound_key_level(myid,levelmin)-1
     hk0(1)=ikey
     ! Compute Cartesian index from Hilbert index
     If(ndim==1)then
        call hilbert1d_reverse(ix,hk0,1)
     else if(ndim==2)then
        call hilbert2d_reverse(ix,iy,hk1,hk0,levelmin-1,1)
     else if (ndim==3)then
        call hilbert3d_reverse(ix,iy,iz,hk2,hk1,hk0,levelmin-1,1)
     end if
     ! Insert new grid in main array
     igrid=igrid+1
     if(igrid==1)head(levelmin)=1
!     write(*,'("PE ",5(I6,1X))')myid,igrid,ikey,ix(1),iy(1)
     tail(levelmin)=igrid
     noct(levelmin)=noct(levelmin)+1
     noct_used=noct_used+1
     grid(igrid)%lev=levelmin
     grid(igrid)%ckey(1)=ix(1)
#if NDIM>1
     grid(igrid)%ckey(2)=iy(1)
#endif
#if NDIM>2
     grid(igrid)%ckey(3)=iz(1)
#endif
     grid(igrid)%hkey=hk0(1)
     grid(igrid)%refined(1:twotondim)=.false.
     ! Insert new grid in hash table
     hash_key(0)=levelmin
     hash_key(1)=ix(1)
#if NDIM>1
     hash_key(2)=iy(1)
#endif
#if NDIM>2
     hash_key(3)=iz(1)
#endif
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
  key_ref=-1
  do ioct=head(ilev),tail(ilev)
     grid(ioct)%superoct=1
     coarse_key=grid(ioct)%hkey
     do i=1,MIN(ilev-1,nsuperoct)
        coarse_key=coarse_key/twotondim
        if(coarse_key.EQ.key_ref(i))then
           n_same(i)=n_same(i)+1
        else
           n_same(i)=1
           key_ref(i)=coarse_key
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
  integer::ilevel,i,ivar, ilev
  logical :: use_histograms 

  if(myid==1)write(*,*)'Building initial adaptive grid'

  do i=levelmin,nlevelmax+1

     call refine_fine(levelmin)

     do ilevel=nlevelmax,levelmin,-1
        if(hydro)then
           call init_flow_fine(ilevel)
           call upload_fine(ilevel)
        endif
     end do

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
