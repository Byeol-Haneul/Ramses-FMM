!###########################################################
!###########################################################
!###########################################################
!###########################################################
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
#ifdef GRAV
  do igrid=head(ilevel),tail(ilevel)
     do ind=1,twotondim
        rhs_norm=rhs_norm+fact2*(grid(igrid)%rho(ind)-rho_tot)**2
     end do
  end do
#endif
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
#ifdef GRAV
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           r2=r2+grid(igrid)%f(ind,1)**2
        end do
     end do
#endif
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
#ifdef GRAV
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           grid(igrid)%f(ind,2)=grid(igrid)%f(ind,1)+beta_cg*grid(igrid)%f(ind,2)
        end do
     end do
#endif

     !==============================================
     ! Compute z = Ap and store it into f(i,3)
     !==============================================
     call cmp_Ap_cg(ilevel)

     !====================================
     ! Compute p.Ap scalar product
     !====================================
     pAp=0.0d0
#ifdef GRAV
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           pAp=pAp+grid(igrid)%f(ind,2)*grid(igrid)%f(ind,3)
        end do
     end do
#endif
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
#ifdef GRAV
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           grid(igrid)%phi(ind)=grid(igrid)%phi(ind)+alpha_cg*grid(igrid)%f(ind,2)
        end do
     end do
#endif

     !====================================
     ! Recurrence on r
     !====================================
#ifdef GRAV
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           grid(igrid)%f(ind,1)=grid(igrid)%f(ind,1)-alpha_cg*grid(igrid)%f(ind,3)
        end do
     end do
#endif

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

  call open_cache(operation_phi,domain_decompos_amr)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=head(ilevel),tail(ilevel)

#ifdef GRAV
     ! Get central oct potential
     do ind=1,twotondim
        phi_nbor(ind,0)=grid(igrid)%phi(ind)
     end do
#endif

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
#ifdef GRAV
           do ind=1,twotondim
              phi_nbor(ind,i_nbor)=grid(igridn)%phi(ind)
           end do
#endif

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
#ifdef GRAV
        residu=grid(igrid)%phi(ind)
#endif
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           residu=residu-oneoversix*(phi_nbor(id1,ig1)+phi_nbor(id2,ig2))
        end do
#ifdef GRAV
        residu=residu+fact*(grid(igrid)%rho(ind)-rho_tot)
#endif

        ! Store results in f(ind,1)
#ifdef GRAV
        grid(igrid)%f(ind,1)=residu
#endif

        ! Store results in f(ind,2)
#ifdef GRAV
        grid(igrid)%f(ind,2)=residu
#endif

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
  integer::i_nbor,igrid,idim,ind,igridn
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
#ifdef GRAV
     do ind=1,twotondim
        phi_nbor(ind,0)=grid(igrid)%phi(ind)
     end do
#endif

     ! Get neighboring octs potential
     do i_nbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=grid(igrid)%ckey(1:ndim)+shift(1:ndim,i_nbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
           if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
        enddo

        ! Get neighbouring grid using read-only cache
        igridn=get_grid(hash_nbor,grid_dict,.false.,.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
#ifdef GRAV
           do ind=1,twotondim
              phi_nbor(ind,i_nbor)=grid(igridn)%f(ind,2)
           end do
#endif
        else
           do ind=1,twotondim
              phi_nbor(ind,i_nbor)=0.0
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim

        ! Compute Ap using neighbors potential
#ifdef GRAV
        residu=-grid(igrid)%f(ind,2)
#endif
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           residu=residu+oneoversix*(phi_nbor(id1,ig1)+phi_nbor(id2,ig2))
        end do

        ! Store results in f(ind,3)
#ifdef GRAV
        grid(igrid)%f(ind,3)=residu
#endif

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

  call open_cache(operation_phi,domain_decompos_amr)

  hash_key(0)=ilevel

  ! Loop over grids
  do igrid=head(ilevel),tail(ilevel)

     ! By default, initial phi is equal to zero

#ifdef GRAV
     ! Loop over cells
     do ind=1,twotondim
        grid(igrid)%phi(ind)=0.0d0
        do idim=1,ndim
           grid(igrid)%f(ind,idim)=0.0
        end do
     end do
     ! End loop over cells
#endif

     ! For fine levels, initial phi is interpolated from coarser level
     if(ilevel.GT.levelmin)then
        
        hash_key(1:ndim)=grid(igrid)%ckey(1:ndim)
        ! Get 3**ndim neghbouring parent cell using read-only cache
        call get_threetondim_nbor_parent_cell(hash_key,grid_dict,igrid_nbor,ind_nbor,.false.,.true.)
        call interpol_phi(igrid_nbor,ind_nbor,ccc,bbb,tfrac,phi_int)
        do ind=1,threetondim
           call unlock_cache(igrid_nbor(ind))
        end do

#ifdef GRAV
        ! Loop over cells
        do ind=1,twotondim
           grid(igrid)%phi(ind)=phi_int(ind)
        end do
        ! End loop over cells
#endif

     end if

  end do
  ! End loop over grids

  call close_cache(grid_dict)

end subroutine make_initial_phi
!###########################################################
!###########################################################
!###########################################################
!###########################################################



