#ifdef GRAV
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Mask restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_mask_2(r,g,m,ifinelevel,allmasked)
  use amr_parameters, only: dp,nvector,nhilbert,ndim,twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  use cache_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
  integer::info
  logical::allmasked_tot
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer, intent(in) :: ifinelevel
  logical, intent(inout) :: allmasked

  integer(kind=8),dimension(0:ndim) :: hash_key
  integer :: ichild,ind,igrid,icell
  integer :: parent_cell, get_parent_cell_2
  real(dp) :: ngpmask, mask_max
  real(dp) :: dtwotondim = (twotondim)
  
  ! Initialize volume fraction to zero at coarse level
  do igrid=m%head_mg(ifinelevel-1),m%tail_mg(ifinelevel-1)
     do ind=1,twotondim
        m%grid(igrid)%f(ind,3)=0d0
     end do
  end do
  
  hash_key(0)=ifinelevel
  
  call open_cache_2(r,g,m,operation_restrict_mask,domain_decompos_mg)
  
  ! Loop over grids
  do ichild=m%head_mg(ifinelevel),m%tail_mg(ifinelevel)

     ! Loop over cells
     do ind=1,twotondim        

        hash_key(1:ndim)=m%grid(ichild)%ckey(1:ndim)

        ! Get parent cell using write-only cache
        parent_cell=get_parent_cell_2(r,g,m,hash_key,m%mg_dict,.true.,.false.)
        igrid=(parent_cell-1)/twotondim+1
        icell=parent_cell-(igrid-1)*twotondim

        ! Convert mask value to volume fraction
        ngpmask=(1d0+m%grid(ichild)%f(ind,3))/2d0/dtwotondim
        m%grid(igrid)%f(icell,3)=m%grid(igrid)%f(icell,3)+ngpmask

     end do
  end do

  call close_cache_2(r,g,m,m%mg_dict)

  ! Convert volume fraction back to to mask value for coarse level
  do igrid=m%head_mg(ifinelevel-1),m%tail_mg(ifinelevel-1)
     do ind=1,twotondim
        m%grid(igrid)%f(ind,3)=2d0*m%grid(igrid)%f(ind,3)-1d0
     end do
  end do
  
  ! Check mask state at coarse level
  mask_max=-1.0
  do igrid=m%head_mg(ifinelevel-1),m%tail_mg(ifinelevel-1)
     do ind=1,twotondim
        mask_max=MAX(mask_max,m%grid(igrid)%f(ind,3))
     end do
  end do
  allmasked=(mask_max<=0d0)
  
  ! Allreduce on mask state
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(allmasked, allmasked_tot, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, info)
  allmasked=allmasked_tot
#endif

end subroutine restrict_mask_2

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Residual computation
! ------------------------------------------------------------------------

subroutine cmp_residual_mg_2(r,g,m,hash_dict, ilevel)
  use amr_parameters, only: dp,nvector,nhilbert,ndim,twondim,twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  use cache_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer, intent(in) :: ilevel
  type(hash_table) :: hash_dict

  ! Computes the residual for MG levels, and stores it into grid(igrid)%f(ind,1)
    
  integer :: get_grid_2
  integer, dimension(1:3,1:2,1:8) :: iii, jjj
  real(dp),dimension(1:twotondim,0:twondim),save::phi_nbor,dis_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer(kind=8),dimension(0:ndim) :: hash_nbor
  real(dp) :: dx, oneoverdx2, phi_c, dis_c, nb_sum
  integer  :: igrid, ind, inbor, idim, igridn, id, ig
  real(dp) :: dtwondim = (twondim)
  
  ! Set constants
  dx = r%boxlen/2**ilevel
  oneoverdx2 = 1.0d0/(dx*dx)
  
  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  
  call open_cache_2(r,g,m,operation_mg,domain_decompos_mg)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=m%head_mg(ilevel),m%tail_mg(ilevel)

     ! Get central oct potential
     do ind=1,twotondim
        phi_nbor(ind,0)=m%grid(igrid)%phi(ind)
        dis_nbor(ind,0)=m%grid(igrid)%f(ind,3)
     end do

     ! Get neighboring octs potential
     do inbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=m%grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=m%ckey_max(ilevel)-1
           if(hash_nbor(idim)==m%ckey_max(ilevel))hash_nbor(idim)=0
        enddo

        ! Get neighbouring grid using read-only cache
        igridn=get_grid_2(r,g,m,hash_nbor,hash_dict,.false.,.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=m%grid(igridn)%phi(ind)
              dis_nbor(ind,inbor)=m%grid(igridn)%f(ind,3)
           end do

        ! Otherwise set to zero and outside
        else
           do ind=1,twotondim
              phi_nbor(ind,inbor)=0.0
              dis_nbor(ind,inbor)=-1.0
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim

        ! Compute residual using 6 neighbors potential
        phi_c=m%grid(igrid)%phi(ind)
        dis_c=m%grid(igrid)%f(ind,3)

        nb_sum=0.0

        if(.not. btest(m%grid(igrid)%flag2(ind),0))then ! No scan needed

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 nb_sum=nb_sum+phi_nbor(id,ig)
              end do
           end do

        else ! Scan is required

           ! If cell is outside, set residual to zero
           if(dis_c<=0.0)then
              m%grid(igrid)%f(ind,1)=0.0
              cycle
           else

              ! Loop over neighbours
              do inbor=1,2
                 do idim=1,ndim
                    id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                    if(dis_nbor(id,ig)<=0.0)then
                       nb_sum=nb_sum+phi_c/dis_c*dis_nbor(id,ig)
                    else
                       nb_sum=nb_sum+phi_nbor(id,ig)
                    endif
                 end do
              end do

           endif

        endif

        ! Store residual in f(ind,1)
        m%grid(igrid)%f(ind,1)=-oneoverdx2*( nb_sum - dtwondim*phi_c )+m%grid(igrid)%f(ind,2)

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache_2(r,g,m,hash_dict)

end subroutine cmp_residual_mg_2

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Residual computation using fast communication method
! ------------------------------------------------------------------------
#ifdef FAST
subroutine cmp_residual_mg_fast(hash_dict, ilevel)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
  integer::info
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer, intent(in) :: ilevel
  type(hash_table) :: hash_dict

  ! Computes the residual for MG levels, and stores it into grid(igrid)%f(ind,1)
    
  integer :: get_grid
  integer, dimension(1:3,1:2,1:8) :: iii, jjj
  real(dp),dimension(1:twotondim,0:twondim),save::phi_nbor,dis_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer(kind=8),dimension(0:ndim) :: hash_nbor
  real(dp) :: dx, oneoverdx2, phi_c, dis_c, nb_sum
  integer  :: igrid, ind, inbor, idim, igridn, id, ig
  real(dp) :: dtwondim = (twondim)
  integer  :: icpu,i,istart,nbuffer,countrecv,countsend,tag=101
  integer,dimension(ncpu) :: reqsend,reqrecv  

  ! Set constants
  dx = boxlen/2**ilevel
  oneoverdx2 = 1.0d0/(dx*dx)
  
  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  
#ifndef WITHOUTMPI

  ! Update boundary conditions
  countrecv=0
  do icpu=1,ncpu
     nbuffer=buffer_mg(ilevel)%send_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=buffer_mg(ilevel)%send_oft(icpu)*twotondim+1
        call MPI_IRECV(buffer_mg(ilevel)%phi_send_buf(istart),nbuffer*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do

  do ind=1,twotondim
     do i=1,buffer_mg(ilevel)%recv_tot
        igrid=buffer_mg(ilevel)%grid_recv_buf(i)
        istart=(i-1)*twotondim+ind
        buffer_mg(ilevel)%phi_recv_buf(istart)=grid(igrid)%phi(ind)
     end do
  end do

  countsend=0
  do icpu=1,ncpu
     nbuffer=buffer_mg(ilevel)%recv_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=buffer_mg(ilevel)%recv_oft(icpu)*twotondim+1
        call MPI_ISEND(buffer_mg(ilevel)%phi_recv_buf(istart),nbuffer*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  do ind=1,twotondim
     do i=1,buffer_mg(ilevel)%send_tot
        istart=(i-1)*twotondim+ind
        buffer_mg(ilevel)%phi_remote(i,ind)=buffer_mg(ilevel)%phi_send_buf(istart)
     end do
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

  ! Loop over grids
  do igrid=head_mg(ilevel),tail_mg(ilevel)

     ! Get central oct potential
     do ind=1,twotondim
        phi_nbor(ind,0)=grid(igrid)%phi(ind)
        dis_nbor(ind,0)=grid(igrid)%f(ind,3)
     end do

     ! Get neighboring octs potential
     do inbor=1,twondim

        ! Get neighbouring grid using read-only cache
        igridn=buffer_mg(ilevel)%nbor_indx(igrid,inbor)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=grid(igridn)%phi(ind)
              dis_nbor(ind,inbor)=grid(igridn)%f(ind,3)
           end do
        else if (igridn==0)then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=0.0
              dis_nbor(ind,inbor)=-1.0
           end do
        else
           do ind=1,twotondim
              phi_nbor(ind,inbor)=buffer_mg(ilevel)%phi_remote(-igridn,ind)
              dis_nbor(ind,inbor)=buffer_mg(ilevel)%dis_remote(-igridn,ind)
           end do           
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim

        ! Compute residual using 6 neighbors potential
        phi_c=grid(igrid)%phi(ind)
        dis_c=grid(igrid)%f(ind,3)

        nb_sum=0.0

        if(.not. btest(grid(igrid)%flag2(ind),0))then ! No scan needed

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 nb_sum=nb_sum+phi_nbor(id,ig)
              end do
           end do

        else ! Scan is required

           ! If cell is outside, set residual to zero
           if(dis_c<=0.0)then
              grid(igrid)%f(ind,1)=0.0
              cycle
           else

              ! Loop over neighbours
              do inbor=1,2
                 do idim=1,ndim
                    id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                    if(dis_nbor(id,ig)<=0.0)then
                       nb_sum=nb_sum+phi_c/dis_c*dis_nbor(id,ig)
                    else
                       nb_sum=nb_sum+phi_nbor(id,ig)
                    endif
                 end do
              end do

           endif

        endif

        ! Store residual in f(ind,1)
        grid(igrid)%f(ind,1)=-oneoverdx2*( nb_sum - dtwondim*phi_c )+grid(igrid)%f(ind,2)

     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine cmp_residual_mg_fast
#endif
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Gauss-Seidel Red-Black sweeps
! ------------------------------------------------------------------------

subroutine gauss_seidel_mg_2(r,g,m,hash_dict,ilevel,safe,redstep)
  use amr_parameters, only: dp,nvector,nhilbert,ndim,twondim,twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  use cache_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer, intent(in) :: ilevel
  logical, intent(in) :: safe
  logical, intent(in) :: redstep
  type(hash_table) :: hash_dict

  ! Perform a Gauss-Seidel update of grid(igrid)%phi(ind).
  ! The domain mask is also needed.
  
  integer :: get_grid_2
  integer, dimension(1:3,1:2,1:8) :: iii, jjj
  real(dp),dimension(1:twotondim,0:twondim),save::phi_nbor,dis_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer(kind=8),dimension(0:ndim) :: hash_nbor
  real(dp) :: phi_c, dis_c, dx2, nb_sum, weight
  integer  :: igrid, ind, inbor, idim, igridn, id, ig, ind0
  real(dp) :: dtwondim = (twondim)

  integer, dimension(1:4) :: ired, iblack
  
  ! Set constants
  dx2  = (r%boxlen/2**ilevel)**2
  
  ired  (1:4)=(/1,4,6,7/)
  iblack(1:4)=(/2,3,5,8/)
  
  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  
  call open_cache_2(r,g,m,operation_mg,domain_decompos_mg)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=m%head_mg(ilevel),m%tail_mg(ilevel)
     
     ! Loop over cells
     do ind=1,twotondim

        ! Get central oct potential and distance
        phi_nbor(ind,0)=m%grid(igrid)%phi(ind)
        dis_nbor(ind,0)=m%grid(igrid)%f(ind,3)

     end do

     ! Get neighboring octs potential
     do inbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=m%grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=m%ckey_max(ilevel)-1
           if(hash_nbor(idim)==m%ckey_max(ilevel))hash_nbor(idim)=0
        enddo

        ! Get neighbouring grid using a read-only cache
        igridn=get_grid_2(r,g,m,hash_nbor,hash_dict,.false.,.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=m%grid(igridn)%phi(ind)
              dis_nbor(ind,inbor)=m%grid(igridn)%f(ind,3)
           end do

        ! Otherwise set to zero and outside
        else
           
           do ind=1,twotondim
              phi_nbor(ind,inbor)=0.0
              dis_nbor(ind,inbor)=-1.0
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells, with red/black ordering
     do ind0=1,twotondim/2      ! Only half of the cells for a red or black sweep

        if(redstep) then
           ind = ired  (ind0)
        else
           ind = iblack(ind0)
        end if

        ! Compute residual using 6 neighbors potential
        phi_c=m%grid(igrid)%phi(ind)
        dis_c=m%grid(igrid)%f(ind,3)

        nb_sum=0.0

        if(.not. btest(m%grid(igrid)%flag2(ind),0))then ! No scan needed

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 nb_sum=nb_sum+phi_nbor(id,ig)
              end do
           end do

           ! Update the potential, solving for potential on icell_amr
           m%grid(igrid)%phi(ind)=(nb_sum-dx2*m%grid(igrid)%f(ind,2))/dtwondim

        else ! Scan is required

           ! If cell is outside, don't update phi
           if(dis_c<=0.0)cycle
           if(safe .and. dis_c<1.0)cycle

           weight=0.0d0   ! Central weight for "Solve G-S"

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 if(dis_nbor(id,ig)<=0.0)then
                    weight=weight+dis_nbor(id,ig)/dis_c
                 else
                    nb_sum=nb_sum+phi_nbor(id,ig)
                 endif
              end do
           end do

           ! Update the potential
           m%grid(igrid)%phi(ind)=(nb_sum-dx2*m%grid(igrid)%f(ind,2))/(dtwondim - weight)

        endif

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache_2(r,g,m,hash_dict)

end subroutine gauss_seidel_mg_2

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Gauss-Seidel Red-Black sweeps
! ------------------------------------------------------------------------
#ifdef FAST
subroutine gauss_seidel_mg_fast(hash_dict,ilevel,safe,redstep)
  use amr_commons
  use poisson_commons
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
  integer::info
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer, intent(in) :: ilevel
  logical, intent(in) :: safe
  logical, intent(in) :: redstep
  type(hash_table) :: hash_dict

  ! Perform a Gauss-Seidel update of grid(igrid)%phi(ind).
  ! The domain mask is also needed.
  
  integer :: get_grid
  integer, dimension(1:3,1:2,1:8) :: iii, jjj
  real(dp),dimension(1:twotondim,0:twondim),save::phi_nbor,dis_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer(kind=8),dimension(0:ndim) :: hash_nbor
  real(dp) :: oneoverdx2, phi_c, dis_c, dx2, nb_sum, weight
  integer  :: igrid, ind, inbor, idim, igridn, id, ig, ind0, ipos
  real(dp) :: dtwondim = (twondim)
  integer  :: icpu,i,istart,nbuffer,countrecv,countsend,tag=101
  integer,dimension(ncpu) :: reqsend,reqrecv  

  integer, dimension(1:4) :: ired, iblack
  
  ! Set constants
  dx2  = (boxlen/2**ilevel)**2
  
  ired  (1:4)=(/1,4,6,7/)
  iblack(1:4)=(/2,3,5,8/)
  
  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  
#ifndef WITHOUTMPI

  ! Update boundary conditions
  countrecv=0
  do icpu=1,ncpu
     nbuffer=buffer_mg(ilevel)%send_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=buffer_mg(ilevel)%send_oft(icpu)*twotondim+1
        call MPI_IRECV(buffer_mg(ilevel)%phi_send_buf(istart),nbuffer*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do

  do ind=1,twotondim
     do i=1,buffer_mg(ilevel)%recv_tot
        igrid=buffer_mg(ilevel)%grid_recv_buf(i)
        istart=(i-1)*twotondim+ind
        buffer_mg(ilevel)%phi_recv_buf(istart)=grid(igrid)%phi(ind)
     end do
  end do

  countsend=0
  do icpu=1,ncpu
     nbuffer=buffer_mg(ilevel)%recv_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=buffer_mg(ilevel)%recv_oft(icpu)*twotondim+1
        call MPI_ISEND(buffer_mg(ilevel)%phi_recv_buf(istart),nbuffer*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  do ind=1,twotondim
     do i=1,buffer_mg(ilevel)%send_tot
        istart=(i-1)*twotondim+ind
        buffer_mg(ilevel)%phi_remote(i,ind)=buffer_mg(ilevel)%phi_send_buf(istart)
     end do
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

  ! Loop over grids
  do igrid=head_mg(ilevel),tail_mg(ilevel)
     
     ! Loop over cells
     do ind=1,twotondim

        ! Get central oct potential and distance
        phi_nbor(ind,0)=grid(igrid)%phi(ind)
        dis_nbor(ind,0)=grid(igrid)%f(ind,3)

     end do

     ! Get neighboring octs potential
     do inbor=1,twondim

        ! Get neighbouring grid using read-only cache
        igridn=buffer_mg(ilevel)%nbor_indx(igrid,inbor)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=grid(igridn)%phi(ind)
              dis_nbor(ind,inbor)=grid(igridn)%f(ind,3)
           end do
        else if (igridn==0)then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=0.0
              dis_nbor(ind,inbor)=-1.0
           end do
        else
           do ind=1,twotondim
              phi_nbor(ind,inbor)=buffer_mg(ilevel)%phi_remote(-igridn,ind)
              dis_nbor(ind,inbor)=buffer_mg(ilevel)%dis_remote(-igridn,ind)
           end do           
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells, with red/black ordering
     do ind0=1,twotondim/2      ! Only half of the cells for a red or black sweep

        if(redstep) then
           ind = ired  (ind0)
        else
           ind = iblack(ind0)
        end if

        ! Compute residual using 6 neighbors potential
        phi_c=grid(igrid)%phi(ind)
        dis_c=grid(igrid)%f(ind,3)

        nb_sum=0.0

        if(.not. btest(grid(igrid)%flag2(ind),0))then ! No scan needed

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 nb_sum=nb_sum+phi_nbor(id,ig)
              end do
           end do

           ! Update the potential, solving for potential on icell_amr
           grid(igrid)%phi(ind)=(nb_sum-dx2*grid(igrid)%f(ind,2))/dtwondim

        else ! Scan is required

           ! If cell is outside, don't update phi
           if(dis_c<=0.0)cycle
           if(safe .and. dis_c<1.0)cycle

           weight=0.0d0   ! Central weight for "Solve G-S"

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 if(dis_nbor(id,ig)<=0.0)then
                    weight=weight+dis_nbor(id,ig)/dis_c
                 else
                    nb_sum=nb_sum+phi_nbor(id,ig)
                 endif
              end do
           end do

           ! Update the potential
           grid(igrid)%phi(ind)=(nb_sum-dx2*grid(igrid)%f(ind,2))/(dtwondim - weight)

        endif

     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine gauss_seidel_mg_fast
#endif
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Residual restriction
! ------------------------------------------------------------------------

subroutine restrict_residual_2(r,g,m,ifinelevel)
  use amr_parameters, only: dp,nvector,nhilbert,ndim,twondim,twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  use cache_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer, intent(in) :: ifinelevel

  ! Restrict the residual of the fine level (stored in grid(ichild)%f(ind,1))
  ! into the rhs of the coarse level (stored in grid(igrid)%f(icell,2))
  ! For interior coarse cell only (we need the mask stored in grid(igrid)%f(icell,3))
  
  integer :: ichild, ind, get_parent_cell_2
  integer :: igrid, icell, parent_cell
  real(dp) :: dtwotondim = (twotondim)
  integer(kind=8),dimension(0:ndim) :: hash_key
  
  ! Set rhs to zero in coarse cells
  do igrid=m%head_mg(ifinelevel-1),m%tail_mg(ifinelevel-1)
     do ind=1,twotondim
        m%grid(igrid)%f(ind,2)=0.0
     end do
  end do

  hash_key(0)=ifinelevel

  call open_cache_2(r,g,m,operation_restrict_res,domain_decompos_mg)
  
  ! Loop over grids
  do ichild=m%head_mg(ifinelevel),m%tail_mg(ifinelevel)
     
     ! Loop over cells
     do ind=1,twotondim        
        
        ! Is fine cell masked?
        if(m%grid(ichild)%f(ind,3)<=0d0)cycle

        hash_key(1:ndim)=m%grid(ichild)%ckey(1:ndim)
        
        ! Get parent cell using read-write cache
        parent_cell=get_parent_cell_2(r,g,m,hash_key,m%mg_dict,.true.,.true.)
        igrid=(parent_cell-1)/twotondim+1
        icell=parent_cell-(igrid-1)*twotondim
        
        ! Is coarse cell masked?
        if(m%grid(igrid)%f(icell,3)<=0d0)cycle
        
        ! Stack fine cell residual in coarse cell rhs
        m%grid(igrid)%f(icell,2)=m%grid(igrid)%f(icell,2)+m%grid(ichild)%f(ind,1)/dtwotondim
        
     end do
  end do
  
  call close_cache_2(r,g,m,m%mg_dict)

end subroutine restrict_residual_2

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Interpolation and correction
! ------------------------------------------------------------------------

subroutine interpolate_and_correct_2(r,g,m,ifinelevel)
  use amr_parameters, only: dp,nvector,nhilbert,ndim,twondim,twotondim,threetondim
  use amr_commons, only: run_t,global_t,mesh_t
  use cache_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer, intent(in) :: ifinelevel
  
  ! Interpolate the solution of the coarse level (stored in grid(igrid)%phi(icell))
  ! and corrct the solutionn of the fine level (stored in grid(ichild)%phi(ind))
  
  integer(kind=8),dimension(0:ndim) :: hash_key
  integer,dimension(1:threetondim) :: igrid_nbor,ind_nbor
  integer  :: ichild, ind
  real(dp) :: aa, bb, cc, dd, coeff
  real(dp), dimension(1:8)     :: bbb
  integer,  dimension(1:8,1:8) :: ccc
  integer::ind_average,ind_father
  integer::igrid_nbr,ind_nbr
  real(dp),dimension(1:twotondim)::corr
  
  ! Local constants
  aa = 1.0D0/4.0D0**ndim
  bb = 3.0D0*aa
  cc = 9.0D0*aa
  dd = 27.D0*aa 
  bbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)
  
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)
  
  if(r%verbose)write(*,*)'entering interpolate  and correct ',ifinelevel

  call open_cache_2(r,g,m,operation_phi,domain_decompos_mg)

  hash_key(0)=ifinelevel

  ! Loop over grids
  do ichild=m%head_mg(ifinelevel),m%tail_mg(ifinelevel)

     ! For fine level, correction is interpolated from coarser level solution
     hash_key(1:ndim)=m%grid(ichild)%ckey(1:ndim)
     
     ! Get 3**ndim neighbouring parent cell using a read-only cache
     call get_threetondim_nbor_parent_cell_2(r,g,m,hash_key,m%mg_dict,igrid_nbor,ind_nbor,.false.,.true.)
     
     ! Loop over cells
     do ind=1,twotondim

        ! Set correction to zero
        corr(ind)=0d0

        ! Fine cell is masked as "outside": no correction
        if(m%grid(ichild)%f(ind,3)<=0.0)cycle

        ! Loop over relevant parent cells
        do ind_average=1,twotondim
           ind_father=ccc(ind_average,ind)
           coeff=bbb(ind_average)
           igrid_nbr=igrid_nbor(ind_father)
           ind_nbr=ind_nbor(ind_father)
           if (igrid_nbr>0) then
              corr(ind)=corr(ind)+coeff*m%grid(igrid_nbr)%phi(ind_nbr)
           endif
        end do

     end do
     ! End loop over cells
     
     do ind=1,threetondim
        call unlock_cache_2(r,g,m,igrid_nbor(ind))
     end do
     
     ! Add correction to fine level solution
     do ind=1,twotondim
        m%grid(ichild)%phi(ind)=m%grid(ichild)%phi(ind)+corr(ind)
     end do

  end do
  ! End loop over grids

  call close_cache_2(r,g,m,m%mg_dict)

end subroutine interpolate_and_correct_2

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Flag settings used to speed-up the sweeps
! ------------------------------------------------------------------------

subroutine set_scan_flag_2(r,g,m,hash_dict,ilevel)
  use amr_parameters, only: dp,nvector,nhilbert,ndim,twondim,twotondim,threetondim
  use amr_commons, only: run_t,global_t,mesh_t
  use cache_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer, intent(in) :: ilevel
  type(hash_table) :: hash_dict
  !
  integer :: get_grid_2
  integer :: ind, igrid, igridn, inbor, idim, id, ig
  integer, dimension(1:3,1:2,1:8)::iii, jjj
  real(dp),dimension(1:twotondim,0:twondim)::dis_nbor
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  real(dp)::dis_c
  
  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  
  call open_cache_2(r,g,m,operation_scan,domain_decompos_mg)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=m%head_mg(ilevel),m%tail_mg(ilevel)

     ! Get central oct potential
     do ind=1,twotondim
        dis_nbor(ind,0)=m%grid(igrid)%f(ind,3)
     end do

     ! Get neighboring octs potential
     do inbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=m%grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=m%ckey_max(ilevel)-1
           if(hash_nbor(idim)==m%ckey_max(ilevel))hash_nbor(idim)=0
        enddo

        ! Get neighbouring grid using read-only cache
        igridn=get_grid_2(r,g,m,hash_nbor,hash_dict,.false.,.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              dis_nbor(ind,inbor)=m%grid(igridn)%f(ind,3)
           end do

        ! Otherwise set to "outside"
        else
           do ind=1,twotondim
              dis_nbor(ind,inbor)=-1.0
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim

        ! Compute residual using 6 neighbors potential
        dis_c=m%grid(igrid)%f(ind,3)

        ! If cell is entirely inside, set flag tentatively to 0 (no scan needed)
        if(dis_c==1.0)then
           m%grid(igrid)%flag2(ind)=0

           ! Loop over neighbours
           do inbor=1,2
              do idim=1,ndim
                 id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
                 ! If one neighbour is outside, then scan needed
                 if(dis_nbor(id,ig)<=0.0)then
                    m%grid(igrid)%flag2(ind)=max(m%grid(igrid)%flag2(ind),1)
                 endif
              end do
           end do

         ! If cell is even partially outside, then scan needed
        else
           m%grid(igrid)%flag2(ind)=1           
        endif

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache_2(r,g,m,hash_dict)

end subroutine set_scan_flag_2

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Compute norm of residual 
! ------------------------------------------------------------------------

subroutine cmp_residual_norm2_2(r,m,ilevel, norm2)
  use amr_parameters, only: dp,ndim,twotondim
  use amr_commons, only: run_t,mesh_t
  implicit none
  type(run_t)::r
  type(mesh_t)::m
  integer,  intent(in)  :: ilevel
  real(kind=8), intent(out) :: norm2
  
  real(kind=8) :: dx2
  integer  :: ind, igrid
  
  ! Set constants
  dx2  = (r%boxlen/2**ilevel)**2
  norm2 = 0.0d0
  
  ! Loop over grids
  do igrid=m%head_mg(ilevel),m%tail_mg(ilevel)
     
     ! Loop over cells
     do ind=1,twotondim
        if(m%grid(igrid)%f(ind,3)<=0.0)cycle      ! Do not count masked cells
        norm2 = norm2 + m%grid(igrid)%f(ind,1)**2
     end do
     ! End loop over cells
     
  end do
  ! End loop over grids
  
  norm2 = dx2*norm2
  
end subroutine cmp_residual_norm2_2

#endif

