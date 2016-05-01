!###########################################################
!###########################################################
!###########################################################
!###########################################################
#ifdef GRAV
subroutine phi_fine_cg(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !=========================================================
  ! Iterative Poisson solver with Conjugate Gradient method 
  ! to solve A x = b
  ! r  : stored in f(1)
  ! p  : stored in f(2)
  ! A p: stored in f(3)
  ! x  : stored in phi
  ! b  : stored in rho
  !=========================================================
  integer::i,igrid,idim,info,ind,iter,itermax
  real(dp)::error,error_ini
  real(dp)::dx2,fourpi,oneoversix,fact,fact2
  real(dp)::r2_old,alpha_cg,beta_cg
  real(kind=8)::r2,pAp,rhs_norm,r2_all,pAp_all,rhs_norm_all

  if(gravity_type>0)return
  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set constants
  dx2=(boxlen/2**ilevel)**2
  fourpi=4.D0*ACOS(-1.0D0)
  if(cosmo)fourpi=1.5D0*omega_m*aexp
  oneoversix=1.0D0/dble(twondim)
  fact=oneoversix*fourpi*dx2
  fact2=fact*fact

  !===============================
  ! Compute initial phi
  !===============================
  call make_initial_phi(ilevel,icount)

  !===============================
  ! Compute right-hand side norm
  !===============================
  rhs_norm=0.d0
  do igrid=head(ilevel),tail(ilevel)
     do ind=1,twotondim
        rhs_norm=rhs_norm+fact2*(grid(igrid)%rho(ind)-rho_tot)**2
     end do
  end do
  ! Compute global norms
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rhs_norm,rhs_norm_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  rhs_norm=rhs_norm_all
#endif
  rhs_norm=DSQRT(rhs_norm/dble(twotondim*noct_tot(ilevel)))

  !==============================================
  ! Compute r = b - Ax and store it into f(i,1)
  ! Also set p = r and store it into f(i,2)
  !==============================================
  call cmp_residual_cg(ilevel,icount)

  call build_cg(ilevel)

  !====================================
  ! Main iteration loop
  !====================================
  iter=0; itermax=10000
  error=1.0D0; error_ini=1.0D0
  do while(error>epsilon*error_ini.and.iter<itermax)

     iter=iter+1

     !====================================
     ! Compute residual norm
     !====================================
     r2=0.0d0
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           r2=r2+grid(igrid)%f(ind,1)**2
        end do
     end do
     ! Compute global norm
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(r2,r2_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     r2=r2_all
#endif

     !====================================
     ! Compute beta factor
     !====================================
     if(iter==1)then
        beta_cg=0.
     else
        beta_cg=r2/r2_old
     end if
     r2_old=r2

     !====================================
     ! Recurrence on p
     !====================================
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           grid(igrid)%f(ind,2)=grid(igrid)%f(ind,1)+beta_cg*grid(igrid)%f(ind,2)
        end do
     end do

     !==============================================
     ! Compute z = Ap and store it into f(i,3)
     !==============================================
!     call cmp_Ap_cg(ilevel)
     call cmp_Ap_cg_fast(ilevel)

     !====================================
     ! Compute p.Ap scalar product
     !====================================
     pAp=0.0d0
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           pAp=pAp+grid(igrid)%f(ind,2)*grid(igrid)%f(ind,3)
        end do
     end do
     ! Compute global sum
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(pAp,pAp_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     pAp=pAp_all
#endif

     !====================================
     ! Compute alpha factor
     !====================================
     alpha_cg = r2/pAp

     !====================================
     ! Recurrence on x
     !====================================
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           grid(igrid)%phi(ind)=grid(igrid)%phi(ind)+alpha_cg*grid(igrid)%f(ind,2)
        end do
     end do

     !====================================
     ! Recurrence on r
     !====================================
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           grid(igrid)%f(ind,1)=grid(igrid)%f(ind,1)-alpha_cg*grid(igrid)%f(ind,3)
        end do
     end do

     ! Compute error
     error=DSQRT(r2/dble(twotondim*noct_tot(ilevel)))
     if(iter==1)error_ini=error
     if(verbose)write(*,112)iter,error/rhs_norm,error/error_ini

  end do
  ! End main iteration loop

  if(myid==1)write(*,115)ilevel,iter,error/rhs_norm,error/error_ini
  if(iter>=itermax)then
     if(myid==1)write(*,*)'Poisson failed to converge...'
  end if

  call clean_cg

111 format('   Entering phi_fine_cg for level ',I2)
112 format('   ==> Step=',i5,' Error=',2(1pe10.3,1x))
115 format('   ==> Level=',i5,' Step=',i5,' Error=',2(1pe10.3,1x))

end subroutine phi_fine_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_residual_cg(ilevel,icount)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ilevel,icount
  !------------------------------------------------------------------
  ! This routine computes the residual for the Conjugate Gradient
  ! Poisson solver. The residual is stored in f(i,1) and f(i,2).
  !------------------------------------------------------------------
  integer::get_grid
  integer::i_nbor,igrid,idim,ind,igridn
  integer::id1,id2,ig1,ig2
  integer,dimension(1:8,1:8)::ccc
  real(dp)::dx2,fourpi,oneoversix,fact,residu
  real(dp)::aa,bb,cc,dd,tfrac
  real(dp),dimension(1:8)::bbb
  integer,dimension(1:3,1:2,1:8)::iii,jjj
  real(dp),dimension(1:twotondim,0:twondim),save::phi_nbor
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer,dimension(1:threetondim),save::igrid_nbor,ind_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))

  ! Set constants
  dx2=(boxlen/2**ilevel)**2
  fourpi=4.D0*ACOS(-1.0D0)
  if(cosmo)fourpi=1.5D0*omega_m*aexp
  oneoversix=1.0D0/dble(twondim)
  fact=oneoversix*fourpi*dx2

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

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

  call open_cache(operation_interpol,domain_decompos_amr)

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

        ! Get neighbouring grid using a read-only cache
        igridn=get_grid(hash_nbor,grid_dict,.false.,.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,i_nbor)=grid(igridn)%phi(ind)
           end do

        ! Otherwise interpolate from coarser level
        else
           ! Get 3**ndim neighbouring paremt cell using a read-only cache
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

        ! Compute residual using 6 neighbors potential
        residu=grid(igrid)%phi(ind)
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           residu=residu-oneoversix*(phi_nbor(id1,ig1)+phi_nbor(id2,ig2))
        end do
        residu=residu+fact*(grid(igrid)%rho(ind)-rho_tot)

        ! Store results in f(ind,1)
        grid(igrid)%f(ind,1)=residu

        ! Store results in f(ind,2)
        grid(igrid)%f(ind,2)=residu

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache(grid_dict)

end subroutine cmp_residual_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_Ap_cg(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes Ap for the Conjugate Gradient
  ! Poisson Solver and store the result into f(i,3).
  !------------------------------------------------------------------
  integer::get_grid
  integer::inbor,igrid,idim,ind,igridn
  integer::id1,id2,ig1,ig2
  real(dp)::oneoversix,residu
  integer,dimension(1:3,1:2,1:8)::iii,jjj
  real(dp),dimension(1:twotondim,0:twondim),save::phi_nbor
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))

  ! Set constants
  oneoversix=1.0D0/dble(twondim)

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  call open_cache(operation_cg,domain_decompos_amr)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=head(ilevel),tail(ilevel)

     ! Get central oct potential
     do ind=1,twotondim
        phi_nbor(ind,0)=grid(igrid)%f(ind,2)
     end do

     ! Get neighboring octs potential
     do inbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
           if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
        enddo

        ! Get neighbouring grid using read-only cache
        igridn=get_grid(hash_nbor,grid_dict,.false.,.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=grid(igridn)%f(ind,2)
           end do
        else
           do ind=1,twotondim
              phi_nbor(ind,inbor)=0.0
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim

        ! Compute Ap using neighbors potential
        residu=-grid(igrid)%f(ind,2)
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           residu=residu+oneoversix*(phi_nbor(id1,ig1)+phi_nbor(id2,ig2))
        end do

        ! Store results in f(ind,3)
        grid(igrid)%f(ind,3)=residu

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache(grid_dict)

end subroutine cmp_Ap_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_Ap_cg_fast(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes Ap for the Conjugate Gradient
  ! Poisson Solver and store the result into f(i,3).
  !------------------------------------------------------------------
  integer::get_grid
  integer::inbor,igrid,idim,ind,igridn
  integer::id1,id2,ig1,ig2,i,info
  integer::icpu,nbuffer,istart
  real(dp)::oneoversix,residu
  integer,dimension(1:3,1:2,1:8)::iii,jjj
  real(dp),dimension(1:twotondim,0:twondim),save::phi_nbor
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer::countrecv,countsend,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv

  ! Set constants
  oneoversix=1.0D0/dble(twondim)

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

#ifndef WITHOUTMPI
                               call timer('poisson - solver - mpi','start')
  ! Update boundary conditions
  countrecv=0
  do icpu=1,ncpu
     nbuffer=send_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=send_oft(icpu)*twotondim+1
        call MPI_IRECV(phi_send_buf(istart),nbuffer*twotondim, &
             & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do
  
  do ind=1,twotondim
     do i=1,recv_tot
        igrid=grid_recv_buf(i)
        istart=(i-1)*twotondim+ind
        phi_recv_buf(istart)=grid(igrid)%f(ind,2)
     end do
  end do

  countsend=0
  do icpu=1,ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=recv_oft(icpu)*twotondim+1
        call MPI_ISEND(phi_recv_buf(istart),nbuffer*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  do ind=1,twotondim
     do i=1,send_tot
        istart=(i-1)*twotondim+ind
        phi_remote(i,ind)=phi_send_buf(istart)
     end do
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif
                               call timer('poisson - solver','start')
  ! Loop over grids
  do igrid=head(ilevel),tail(ilevel)

     ! Get central oct potential
     do ind=1,twotondim
        phi_nbor(ind,0)=grid(igrid)%f(ind,2)
     end do

     ! Get neighboring octs potential
     do inbor=1,twondim

        ! Get neighbouring grid using read-only cache
        igridn=nbor_indx(igrid,inbor)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=grid(igridn)%f(ind,2)
           end do
        else if (igridn==0) then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=0.0
           end do
        else
           do ind=1,twotondim
              phi_nbor(ind,inbor)=phi_remote(-igridn,ind)
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim

        ! Compute Ap using neighbors potential
        residu=-grid(igrid)%f(ind,2)
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           residu=residu+oneoversix*(phi_nbor(id1,ig1)+phi_nbor(id2,ig2))
        end do

        ! Store results in f(ind,3)
        grid(igrid)%f(ind,3)=residu

     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine cmp_Ap_cg_fast
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine make_initial_phi(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel,icount
  !
  !
  !
  integer::igrid,idim,ind,igridn
  integer,dimension(1:8,1:8)::ccc
  real(dp)::aa,bb,cc,dd,tfrac
  real(dp),dimension(1:8)::bbb
  integer(kind=8),dimension(0:ndim)::hash_key
  integer,dimension(1:threetondim),save::igrid_nbor,ind_nbor
  real(dp),dimension(1:twotondim),save::phi_int

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

  call open_cache(operation_interpol,domain_decompos_amr)

  hash_key(0)=ilevel

  ! Loop over grids
  do igrid=head(ilevel),tail(ilevel)

     ! By default, initial phi is equal to zero

     ! Loop over cells
     do ind=1,twotondim
        grid(igrid)%phi(ind)=0.0d0
        do idim=1,ndim
           grid(igrid)%f(ind,idim)=0.0
        end do
     end do
     ! End loop over cells

     ! For fine levels, initial phi is interpolated from coarser level
     if(ilevel.GT.levelmin)then
        
        hash_key(1:ndim)=grid(igrid)%ckey(1:ndim)
        ! Get 3**ndim neghbouring parent cell using read-only cache
        call get_threetondim_nbor_parent_cell(hash_key,grid_dict,igrid_nbor,ind_nbor,.false.,.true.)
        call interpol_phi(igrid_nbor,ind_nbor,ccc,bbb,tfrac,phi_int)
        do ind=1,threetondim
           call unlock_cache(igrid_nbor(ind))
        end do

        ! Loop over cells
        do ind=1,twotondim
           grid(igrid)%phi(ind)=phi_int(ind)
        end do
        ! End loop over cells

     end if

  end do
  ! End loop over grids

  call close_cache(grid_dict)

end subroutine make_initial_phi
!###########################################################
!###########################################################
!###########################################################
!###########################################################

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
  
subroutine build_cg(ilevel)
  use amr_commons
  use poisson_commons
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  
  integer,intent(in)::ilevel
  !
  integer::get_grid
  integer::icoarselevel,igrid,inbor,idim,ipos,ichild,icpu,grid_cpu,ind,info
  integer::i,igridn,iremote
  integer(kind=8),dimension(0:ndim)::hash_key,hash_father,hash_nbor
  integer,dimension(1:ndim)::cart_key
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer(kind=4),dimension(1:nvector),save::dummy_state
  integer(kind=8),dimension(1:nvector,1:nhilbert),save::hk
  integer(kind=8),dimension(1:nvector,1:ndim),save::ix

#ifndef WITHOUTMPI

  allocate(nremote(1:ncpu))  
  allocate(nbor_indx(head(ilevel):tail(ilevel),1:twondim))

  hash_nbor(0)=ilevel
  
  nremote=0
  call open_cache(operation_cg,domain_decompos_amr)
  
  ! Loop over grids
  do igrid=head(ilevel),tail(ilevel)
     
     ! Gather twondim neighboring grids
     do inbor=1,twondim

        hash_nbor(1:ndim)=grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
           if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
        enddo
        
        ! Get neighbouring grid using read-only cache
        igridn=get_grid(hash_nbor,grid_dict,.false.,.true.)
        
        ! If grid exists, determine if remote
        if(igridn<=ngridmax)then
           nbor_indx(igrid,inbor)=igridn
        else
           ! Compute Cartesian keys of new oct
           cart_key(1:ndim)=hash_nbor(1:ndim)
           
           ! Compute Hilbert keys of new octs
           ix(1,1:ndim)=cart_key(1:ndim)
           call hilbert_key(ix,hk,dummy_state,0,ilevel-1,1)

           ! Determine parent processor and increment counter
           do icpu=1,ncpu
              if(    ge_keys(hk(1,1:nhilbert),bound_key_level(1:nhilbert,icpu-1,ilevel)).AND. &
                   & gt_keys(bound_key_level(1:nhilbert,icpu,ilevel),hk(1,1:nhilbert)))then
                 grid_cpu=icpu
              end if
           end do
           nremote(grid_cpu)=nremote(grid_cpu)+1
        end if

     end do
     ! End loop over neighbors
     
  end do
  ! End loop over grids
  
  call close_cache(grid_dict)

  ! Build communicator
  allocate(nalltoall(1:ncpu,1:ncpu))
  allocate(nalltoall_tot(1:ncpu,1:ncpu))
  nalltoall=0; nalltoall_tot=0
  do icpu=1,ncpu
     nalltoall(myid,icpu)=nremote(icpu)
  end do
  call MPI_ALLREDUCE(nalltoall,nalltoall_tot,ncpu*ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nalltoall=nalltoall_tot
  deallocate(nalltoall_tot)
  allocate(send_cnt(1:ncpu))
  allocate(send_oft(1:ncpu))
  allocate(recv_cnt(1:ncpu))
  allocate(recv_oft(1:ncpu))

  send_cnt=0; send_oft=0; send_tot=0
  recv_cnt=0; recv_oft=0; recv_tot=0
  do icpu=1,ncpu
     send_cnt(icpu)=nalltoall(myid,icpu)
     recv_cnt(icpu)=nalltoall(icpu,myid)
     send_tot=send_tot+send_cnt(icpu)
     recv_tot=recv_tot+recv_cnt(icpu)
     if(icpu<ncpu)then
        send_oft(icpu+1)=send_oft(icpu)+nalltoall(myid,icpu)
        recv_oft(icpu+1)=recv_oft(icpu)+nalltoall(icpu,myid)
     endif
  end do
  deallocate(nalltoall)

#if NDIM>0
  allocate(x_send_buf(1:send_tot))
#endif
#if NDIM>1
  allocate(y_send_buf(1:send_tot))
#endif
#if NDIM>2
  allocate(z_send_buf(1:send_tot))
#endif

  hash_nbor(0)=ilevel

  nremote=0
  call open_cache(operation_cg,domain_decompos_amr)

  ! Loop over grids
  do igrid=head(ilevel),tail(ilevel)

     ! Gather twondim neighboring grids
     do inbor=1,twondim

        hash_nbor(1:ndim)=grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
           if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
        enddo

        ! Get neighbouring grid using read-only cache
        igridn=get_grid(hash_nbor,grid_dict,.false.,.true.)

        ! If grid is remote
        if(igridn>ngridmax)then

           ! Compute Cartesian keys of new oct
           cart_key(1:ndim)=hash_nbor(1:ndim)

           ! Compute Hilbert keys of new octs
           ix(1,1:ndim)=cart_key(1:ndim)
           call hilbert_key(ix,hk,dummy_state,0,ilevel-1,1)

           ! Determine parent processor and increment counter
           do icpu=1,ncpu
              if(    ge_keys(hk(1,1:nhilbert),bound_key_level(1:nhilbert,icpu-1,ilevel)).AND. &
                   & gt_keys(bound_key_level(1:nhilbert,icpu,ilevel),hk(1,1:nhilbert)))then
                 grid_cpu=icpu
              end if
           end do
           nremote(grid_cpu)=nremote(grid_cpu)+1
           iremote=send_oft(grid_cpu)+nremote(grid_cpu)
#if NDIM>0
           x_send_buf(iremote)=cart_key(1)
#endif
#if NDIM>1
           y_send_buf(iremote)=cart_key(2)
#endif
#if NDIM>2
           z_send_buf(iremote)=cart_key(3)
#endif
           ! Store negative index
           nbor_indx(igrid,inbor)=-iremote
        end if

     end do
     ! End loop over neighbors

  end do
  ! End loop over grids

  call close_cache(grid_dict)

  deallocate(nremote)

#if NDIM>0
  allocate(x_recv_buf(1:recv_tot))
  call MPI_ALLTOALLV(x_send_buf,send_cnt,send_oft,MPI_INTEGER, &
       &             x_recv_buf,recv_cnt,recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  deallocate(x_send_buf)
#endif
#if NDIM>1
  allocate(y_recv_buf(1:recv_tot))
  call MPI_ALLTOALLV(y_send_buf,send_cnt,send_oft,MPI_INTEGER, &
       &             y_recv_buf,recv_cnt,recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  deallocate(y_send_buf)
#endif
#if NDIM>2
  allocate(z_recv_buf(1:recv_tot))
  call MPI_ALLTOALLV(z_send_buf,send_cnt,send_oft,MPI_INTEGER, &
       &             z_recv_buf,recv_cnt,recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  deallocate(z_send_buf)
#endif

  allocate(grid_recv_buf(1:recv_tot))

  hash_key(0)=ilevel
  do i=1,recv_tot
     hash_key(1)=x_recv_buf(i)
     hash_key(2)=y_recv_buf(i)
     hash_key(3)=z_recv_buf(i)
     ipos=hash_get(grid_dict,hash_key)
     grid_recv_buf(i)=ipos
  end do
  deallocate(x_recv_buf)
  deallocate(y_recv_buf)
  deallocate(z_recv_buf)

  allocate(phi_send_buf(1:send_tot*twotondim))
  allocate(phi_recv_buf(1:recv_tot*twotondim))
  allocate(phi_remote(1:send_tot,1:twotondim))
#endif

end subroutine build_cg

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine clean_cg
  use amr_commons
  use poisson_commons
  implicit none

#ifndef WITHOUTMPI
  
  deallocate(nbor_indx)
  deallocate(phi_remote)
  deallocate(phi_send_buf)
  deallocate(phi_recv_buf)
  deallocate(grid_recv_buf)
  deallocate(send_cnt)
  deallocate(send_oft)
  deallocate(recv_cnt)
  deallocate(recv_oft)

#endif

end subroutine clean_cg
#endif

