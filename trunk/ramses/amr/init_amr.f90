subroutine init_amr
  use amr_commons
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif
  integer(kind=4)::ilevel,icpu
  integer(kind=8)::maxkey
  logical::ok

  if(verbose.and.myid==1)write(*,*)'Entering init_amr'

  ! Initial time step for each level
  dtold=0.0D0
  dtnew=0.0D0

  ! Allocate main oct array
  allocate(grid(1:ngridmax))

  ! Allocate oct hash table
  if(verbose.and.myid==1)write(*,*)'Initialize empty hash'
  call init_empty_hash(grid_dict,2*ngridmax)

  ! Set initial cpu boundaries
  ! Set maximum Cartesian key per level
  if(verbose.and.myid==1)write(*,*)'Initialize level cpu boundaries'
  allocate(bound_key_level(0:ncpu,levelmin:nlevelmax))
  allocate(ckey_max(levelmin:nlevelmax))
  do ilevel=levelmin,nlevelmax
     ckey_max(ilevel)=2**(ilevel-1)
     maxkey=2**((ilevel-1)*ndim)
     do icpu=1,ncpu-1
        bound_key_level(icpu,ilevel) = (icpu*maxkey)/ncpu
     end do
     bound_key_level(0,ilevel) = 0
     bound_key_level(ncpu,ilevel) = maxkey
  end do

  ! Allocate head, tail and numbers for each level
  if(verbose.and.myid==1)write(*,*)'Initialize oct decomposition'
  allocate(head(levelmin:nlevelmax))
  allocate(tail(levelmin:nlevelmax))
  allocate(noct(levelmin:nlevelmax))
  head=0       ! Head oct in the level
  tail=0       ! Tail oct in the level
  noct=0       ! Number of oct in the level
  noct_used=0  ! Total number of oct used

end subroutine init_amr



