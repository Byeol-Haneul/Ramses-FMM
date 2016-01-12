!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_basegrid
  use amr_commons
  use pm_commons
  use hilbert
  use hash, only:hash_set,hash_get
  implicit none
  !-------------------------------------------
  ! This routine builds the initial AMR grid
  !-------------------------------------------
  integer::ilevel,i,j,k,igrid,ipos
  integer(kind=8)::ikey
  integer(kind=8),dimension(1:nvector)::hk0,hk1,hk2
  integer(kind=8),dimension(1:nvector)::ix=0,iy=0,iz=0
  integer(kind=8),dimension(0:ndim)::hash_key,hash_test

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
