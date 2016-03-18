!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine force_fine(ilevel,icount)
  use amr_commons
  use poisson_commons
   implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ind,i,ngrid,info,idim,nstride
  real(dp)::dx,fact,fourpi
  real(kind=8)::rho_loc,rho_all,epot_loc,epot_all
  real(dp),dimension(1:nvector,1:ndim),save::xx,ff
 
  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel
#ifdef GRAV  
  !-------------------------------------
  ! Compute analytical gravity force
  !-------------------------------------
  if(gravity_type>0)then 
     ! Mesh size at level ilevel in code units
     dx=boxlen/2**ilevel
     ! Loop over grids by vector sweeps
     do igrid=head(ilevel),tail(ilevel),nvector
        ngrid=MIN(nvector,tail(ilevel)-igrid+1)
        ! Loop over cells
        do ind=1,twotondim
           ! Compute cell centre position in code units
           do idim=1,ndim
              nstride=2**(idim-1)
              do i=1,ngrid
                 xx(i,idim)=(2*grid(igrid+i-1)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx
              end do
           end do
           ! Call analytical gravity routine
           call gravana(xx,ff,dx,ngrid)
           ! Scatter variables to main memory
           do idim=1,ndim
              do i=1,ngrid
                 grid(igrid+i-1)%uold(ind,idim)=ff(i,idim)
              end do
           end do
        end do
        ! End loop over cells
     end do
     ! End loop over grids

  !------------------------------
  ! Compute gradient of potential
  !------------------------------
  else
     call gradient_phi(ilevel,icount)
  endif

  !----------------------------------------------
  ! Compute gravity potential and maximum density
  !----------------------------------------------
  ! Mesh size at level ilevel in code units
  dx=boxlen/2**ilevel
  ! Initialise global variables
  rho_loc=0.0; rho_all=0.0
  epot_loc=0.0; epot_all=0.0
  fourpi=4.0D0*ACOS(-1.0D0)
  if(cosmo)fourpi=1.5D0*omega_m*aexp
  fact=-dx**ndim/fourpi/2.0D0

  ! Loop over myid grids by vector sweeps
  do igrid=head(ilevel),tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
        ! Loop over dimensions
        do idim=1,ndim
           if(.not.grid(igrid)%refined(ind))then
              epot_loc=epot_loc+fact*grid(igrid)%f(ind,idim)**2
           endif
        end do
        rho_loc=MAX(rho_loc,dble(abs(grid(igrid)%rho(ind))))
     end do
     ! End loop over cells
  end do
  ! End loop over grids

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(epot_loc,epot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(rho_loc ,rho_all ,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
     epot_loc=epot_all
     rho_loc=rho_all
#endif
     epot_tot=epot_tot+epot_loc
     rho_max(ilevel)=rho_loc

#endif  
111 format('   Entering force_fine for level ',I2)

end subroutine force_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gradient_phi(ilevel,icount)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ilevel,icount
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in the current level grids, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::get_grid
  integer::i_nbor,igrid,idim,ind,igridn
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer,dimension(1:3,1:4,1:8)::ggg,hhh
  integer,dimension(1:8,1:8)::ccc
  integer,dimension(1:threetondim),save::igrid_nbor,ind_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer(kind=8),dimension(0:ndim)::hash_nbor
  real(dp),dimension(1:8)::bbb
  real(dp)::dx,a,b,aa,bb,cc,dd,tfrac
  real(dp)::phi1,phi2,phi3,phi4
  real(dp),dimension(1:twotondim,0:twondim),save::phi_nbor
#ifdef GRAV
  ! Mesh size at level ilevel in code units
  dx=boxlen/2**ilevel

  ! Rescaling factor
  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  !   |dim
  !   | |node
  !   | | |cell
  !   v v v
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! CIC method constants
  aa = 1.0D0/4.0D0**ndim
  bb = 3.0D0*aa
  cc = 9.0D0*aa
  dd = 27.D0*aa
  bbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  ! Sampling positions in the 3x3x3 father cell cube
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  if (icount .ne. 1 .and. icount .ne. 2)then
     write(*,*)'icount has bad value'
     call clean_stop
  endif

  ! Compute fraction of time steps for interpolation
  if (dtold(ilevel-1)>0.0)then
     tfrac=dtnew(ilevel)/dtold(ilevel-1)*(icount-1)
  else
     tfrac=0.0
  end if

  call open_cache(operation_phi,domain_decompos_amr)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=head(ilevel),tail(ilevel)
     
     ! Get central oct potential
     do ind=1,twotondim
        phi_nbor(ind,0)=grid(igrid)%phi(ind)
     end do

     ! Get neighboring octs potential
     do i_nbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=grid(igrid)%ckey(1:ndim)+shift(1:ndim,i_nbor)
        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
           if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
        enddo
        igridn=get_grid(hash_nbor,grid_dict,.false.,.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,i_nbor)=grid(igridn)%phi(ind)
           end do

        ! Otherwise interpolate from coarser level
        else
           ! Get 3**ndim parent cell using read-only cache
           call get_threetondim_nbor_parent_cell(hash_nbor,grid_dict,igrid_nbor,ind_nbor,.false.,.true.)
           call interpol_phi(igrid_nbor,ind_nbor,ccc,bbb,tfrac,phi_nbor(1,i_nbor))
           do ind=1,threetondim
              call unlock_cache(igrid_nbor(ind))
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim

        ! Loop over dimensions
        do idim=1,ndim

           ! Gather nodes indices
           id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind)
           id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind)
           id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind)
           id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind)

           ! Gather potential
           phi1=phi_nbor(id1,ig1)
           phi2=phi_nbor(id2,ig2)
           phi3=phi_nbor(id3,ig3)
           phi4=phi_nbor(id4,ig4)

           ! Compute acceleration
           grid(igrid)%f(ind,idim)=a*(phi1-phi2)-b*(phi3-phi4)

        end do
        ! End loop over dimensions

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache(grid_dict)
#endif
end subroutine gradient_phi





