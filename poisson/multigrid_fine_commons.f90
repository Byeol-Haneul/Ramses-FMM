! ------------------------------------------------------------------------
! Multigrid Poisson solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all generic fine multigrid routines, such as
!   * multigrid iterations @ MG levels
!   * MG workspace building
!
! Used variables:
!     -----------------------------------------------------------------
!     potential            phi     
!     physical RHS         rho     
!     residual             f(:,1)  
!     BC-modified RHS      f(:,2)  
!     mask                 f(:,3)  
!     scan flag            flag2   
!
! ------------------------------------------------------------------------

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################
#ifdef GRAV
! ------------------------------------------------------------------------
! Main multigrid routine, called by amr_step
! ------------------------------------------------------------------------
#ifdef TOTO
subroutine multigrid(ilevel,icount)
  use amr_commons
  use poisson_commons
  use poisson_parameters
  
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  
  integer, intent(in) :: ilevel,icount
  
  integer, parameter  :: MAXITER  = 10
  real(dp), parameter :: SAFE_FACTOR = 0.5
  
  integer  :: igrid, ifine, i, iter, info, icpu
  real(kind=8) :: res_norm2, i_res_norm2, i_res_norm2_tot, res_norm2_tot
  real(kind=8) :: err, last_err
  
  logical :: allmasked
  
  if(gravity_type>0)return
  if(noct_tot(ilevel)==0)return
  
  if(verbose) print '(A,I2)','Entering multigrid at level ',ilevel
  
  ! ---------------------------------------------------------------------
  ! Prepare first guess, mask and BCs at finest level
  ! ---------------------------------------------------------------------
  call make_initial_phi(ilevel,icount)  ! Initial guess
  call make_mask(ilevel)                ! Fill the fine level mask
  call make_bc_rhs(ilevel,icount)       ! Fill BC-modified RHS

  ! ---------------------------------------------------------------------
  ! Initialize Domain Decomposition and Hash Table for Multigrid
  ! ---------------------------------------------------------------------
  call init_mg(ilevel)
  
  ! ---------------------------------------------------------------------
  ! Build Multigrid hierarchy in memory
  ! ---------------------------------------------------------------------
  do ifine=ilevel,2,-1
     call build_mg(ifine)
  end do
  
  ! ---------------------------------------------------------------------
  ! Restrict mask up
  ! ---------------------------------------------------------------------
  levelmin_mg=1
  do ifine=ilevel,2,-1
     ! Restrict and communicate mask
     call restrict_mask(ifine,allmasked)
     if(allmasked) then ! Coarser level is fully masked: stop here
        levelmin_mg=ifine
        exit
     end if
  end do
  
  ! ---------------------------------------------------------------------
  ! Set scan flag (for optimisation)
  ! ---------------------------------------------------------------------
  call set_scan_flag(grid_dict,ilevel)
  do ifine=ilevel-1,levelmin_mg,-1
     call set_scan_flag(mg_dict,ifine)
  end do
  
  ! ---------------------------------------------------------------------
  ! Build communication buffer (for optimisation)
  ! ---------------------------------------------------------------------
  if(fast_solver)then
     allocate(buffer_mg(1:nlevelmax))
     call build_comm_mg(grid_dict,ilevel)
     do ifine=ilevel-1,levelmin_mg,-1
        call build_comm_mg(mg_dict,ifine)
     end do
  endif
  
  ! ---------------------------------------------------------------------
  ! Initiate solve at fine level
  ! ---------------------------------------------------------------------
  
  iter = 0
  err = 1.0d0
  main_iteration_loop: do

     iter=iter+1

     ! Pre-smoothing
     do i=1,ngs_fine
        if(fast_solver)then
           call gauss_seidel_mg_fast(grid_dict,ilevel,safe_mode(ilevel),.true. )  ! Red step
           call gauss_seidel_mg_fast(grid_dict,ilevel,safe_mode(ilevel),.false.)  ! Black step
        else
           call gauss_seidel_mg(grid_dict,ilevel,safe_mode(ilevel),.true. )  ! Red step
           call gauss_seidel_mg(grid_dict,ilevel,safe_mode(ilevel),.false.)  ! Black step
        endif
     end do
     
     ! Compute new residual
     if(fast_solver)then                               
        call cmp_residual_mg_fast(grid_dict,ilevel)
     else
        call cmp_residual_mg(grid_dict,ilevel)
     endif

     ! Compute initial residual norm
     if(iter==1) then
        call cmp_residual_norm2(ilevel,i_res_norm2)
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(i_res_norm2,i_res_norm2_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        i_res_norm2=i_res_norm2_tot
#endif
     end if
     
     if(ilevel>1) then

        ! Restrict residual to coarser level
        call restrict_residual(ilevel)

        ! Reset correction from upper level before solve
        do igrid=head_mg(ilevel-1),tail_mg(ilevel-1)
           grid(igrid)%phi(1:twotondim)=0.0d0
        end do
        
        ! Multigrid-solve the upper level
        call recursive_multigrid(ilevel-1, safe_mode(ilevel))
        
        ! Interpolate coarse solution and correct fine solution
        call interpolate_and_correct(ilevel)

     end if
     
     ! Post-smoothing
     do i=1,ngs_fine
        if(fast_solver)then
           call gauss_seidel_mg_fast(grid_dict,ilevel,safe_mode(ilevel),.true. )  ! Red step
           call gauss_seidel_mg_fast(grid_dict,ilevel,safe_mode(ilevel),.false.)  ! Black step
        else
           call gauss_seidel_mg(grid_dict,ilevel,safe_mode(ilevel),.true. )  ! Red step
           call gauss_seidel_mg(grid_dict,ilevel,safe_mode(ilevel),.false.)  ! Black step
        endif
     end do
     
     ! Update fine residual
     if(fast_solver)then
        call cmp_residual_mg_fast(grid_dict,ilevel)
     else
        call cmp_residual_mg(grid_dict,ilevel)
     endif

     ! Compute residual norm
     call cmp_residual_norm2(ilevel,res_norm2)
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(res_norm2,res_norm2_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     res_norm2=res_norm2_tot
#endif
     
     last_err = err
     err = sqrt(res_norm2/(i_res_norm2+1d-20*rho_tot**2))
     
     ! Verbosity
     if(verbose) print '(A,I5,A,1pE10.3)','   ==> Step=', &
          & iter,' Error=',err
     
     ! Converged?
     if(err<epsilon .or. iter>=MAXITER) exit
     
     ! Not converged, check error and possibly enable safe mode for the level
     if(err > last_err*SAFE_FACTOR .and. (.not. safe_mode(ilevel))) then
        if(verbose)print *,'CAUTION: Switching to safe MG mode for level ',ilevel
        safe_mode(ilevel) = .true.
     end if
     
  end do main_iteration_loop
  
  if(myid==1) print '(A,I5,A,I5,A,1pE10.3)','   ==> Level=',ilevel,' Step=',iter,' Error=',err
  if(myid==1 .and. iter==MAXITER) print *,'WARN: Fine multigrid Poisson failed to converge...'
  
  ! ---------------------------------------------------------------------
  ! Clean communication buffer
  ! ---------------------------------------------------------------------
  if(fast_solver)then
     call clean_comm_mg(ilevel)
     do ifine=ilevel-1,levelmin_mg,-1
        call clean_comm_mg(ifine)
     end do
     deallocate(buffer_mg)
  endif
  
  ! ---------------------------------------------------------------------
  ! Cleanup MG levels after solve complete
  ! ---------------------------------------------------------------------
  call cleanup_mg
  
end subroutine multigrid
#endif
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

subroutine multigrid_2(r,g,m,ilevel,icount)
  use amr_parameters, only: dp,twotondim
  use poisson_parameters, only: ngs_fine,ngs_coarse,ncycles_coarse_safe
  use amr_commons, only: run_t,global_t,mesh_t  
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer, intent(in) :: ilevel,icount
  
  integer, parameter  :: MAXITER  = 10
  real(dp), parameter :: SAFE_FACTOR = 0.5
  
  integer  :: igrid, ifine, i, iter, info, icpu
  real(kind=8) :: res_norm2, i_res_norm2, i_res_norm2_tot, res_norm2_tot
  real(kind=8) :: err, last_err
  
  logical :: allmasked
  
  if(r%gravity_type>0)return
  if(m%noct_tot(ilevel)==0)return
  
  if(r%verbose) print '(A,I2)','Entering multigrid at level ',ilevel
  
  ! ---------------------------------------------------------------------
  ! Prepare first guess, mask and BCs at finest level
  ! ---------------------------------------------------------------------
  call make_initial_phi_2(r,g,m,ilevel,icount)  ! Initial guess
  call make_mask_2(r,g,m,ilevel)                ! Fill the fine level mask
  call make_bc_rhs_2(r,g,m,ilevel,icount)       ! Fill BC-modified RHS

  ! ---------------------------------------------------------------------
  ! Initialize Domain Decomposition and Hash Table for Multigrid
  ! ---------------------------------------------------------------------
  call init_mg_2(r,g,m,ilevel)
  
  ! ---------------------------------------------------------------------
  ! Build Multigrid hierarchy in memory
  ! ---------------------------------------------------------------------
  do ifine=ilevel,2,-1
     call build_mg_2(r,g,m,ifine)
  end do
  
  ! ---------------------------------------------------------------------
  ! Restrict mask up
  ! ---------------------------------------------------------------------
  g%levelmin_mg=1
  do ifine=ilevel,2,-1
     ! Restrict and communicate mask
     call restrict_mask_2(r,g,m,ifine,allmasked)
     if(allmasked) then ! Coarser level is fully masked: stop here
        g%levelmin_mg=ifine
        exit
     end if
  end do
  
  ! ---------------------------------------------------------------------
  ! Set scan flag (for optimisation)
  ! ---------------------------------------------------------------------
  call set_scan_flag_2(r,g,m,m%grid_dict,ilevel)
  do ifine=ilevel-1,g%levelmin_mg,-1
     call set_scan_flag_2(r,g,m,m%mg_dict,ifine)
  end do
  
  ! ---------------------------------------------------------------------
  ! Initiate solve at fine level
  ! ---------------------------------------------------------------------
  
  iter = 0
  err = 1.0d0
  main_iteration_loop: do

     iter=iter+1

     ! Pre-smoothing
     do i=1,ngs_fine
        call gauss_seidel_mg_2(r,g,m,m%grid_dict,ilevel,g%safe_mode(ilevel),.true. )  ! Red step
        call gauss_seidel_mg_2(r,g,m,m%grid_dict,ilevel,g%safe_mode(ilevel),.false.)  ! Black step
     end do
     
     ! Compute new residual
     call cmp_residual_mg_2(r,g,m,m%grid_dict,ilevel)

     ! Compute initial residual norm
     if(iter==1) then
        call cmp_residual_norm2_2(r,g,m,ilevel,i_res_norm2)
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(i_res_norm2,i_res_norm2_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        i_res_norm2=i_res_norm2_tot
#endif
     end if
     
     if(ilevel>1) then

        ! Restrict residual to coarser level
        call restrict_residual_2(r,g,m,ilevel)

        ! Reset correction from upper level before solve
        do igrid=m%head_mg(ilevel-1),m%tail_mg(ilevel-1)
           m%grid(igrid)%phi(1:twotondim)=0.0d0
        end do
        
        ! Multigrid-solve the upper level
        call recursive_multigrid_2(r,g,m,ilevel-1, g%safe_mode(ilevel))
        
        ! Interpolate coarse solution and correct fine solution
        call interpolate_and_correct_2(r,g,m,ilevel)

     end if
     
     ! Post-smoothing
     do i=1,ngs_fine
        call gauss_seidel_mg_2(r,g,m,m%grid_dict,ilevel,g%safe_mode(ilevel),.true. )  ! Red step
        call gauss_seidel_mg_2(r,g,m,m%grid_dict,ilevel,g%safe_mode(ilevel),.false.)  ! Black step
     end do
     
     ! Update fine residual
     call cmp_residual_mg_2(r,g,m,m%grid_dict,ilevel)

     ! Compute residual norm
     call cmp_residual_norm2_2(r,g,m,ilevel,res_norm2)
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(res_norm2,res_norm2_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     res_norm2=res_norm2_tot
#endif
     
     last_err = err
     err = sqrt(res_norm2/(i_res_norm2+1d-20*g%rho_tot**2))
     
     ! Verbosity
     if(r%verbose) print '(A,I5,A,1pE10.3)','   ==> Step=', &
          & iter,' Error=',err
     
     ! Converged?
     if(err<r%epsilon .or. iter>=MAXITER) exit
     
     ! Not converged, check error and possibly enable safe mode for the level
     if(err > last_err*SAFE_FACTOR .and. (.not. g%safe_mode(ilevel))) then
        if(r%verbose)print *,'CAUTION: Switching to safe MG mode for level ',ilevel
        g%safe_mode(ilevel) = .true.
     end if
     
  end do main_iteration_loop
  
  if(g%myid==1) print '(A,I5,A,I5,A,1pE10.3)','   ==> Level=',ilevel,' Step=',iter,' Error=',err
  if(g%myid==1 .and. iter==MAXITER) print *,'WARN: Fine multigrid Poisson failed to converge...'
    
  ! ---------------------------------------------------------------------
  ! Cleanup MG levels after solve complete
  ! ---------------------------------------------------------------------
  call cleanup_mg_2(r,g,m)
  
end subroutine multigrid_2

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Recursive multigrid routine for coarse MG levels
! ------------------------------------------------------------------------
#ifdef TOTO
recursive subroutine recursive_multigrid(ifinelevel, safe)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  
  integer, intent(in) :: ifinelevel
  logical, intent(in) :: safe

  integer :: i, igrid, icpu, info, icycle, ncycle
  
  if(ifinelevel<=levelmin_mg) then
     ! Solve 'directly' :
     do i=1,2*ngs_coarse
        if(fast_solver)then
           call gauss_seidel_mg_fast(mg_dict,ifinelevel,safe,.true. )  ! Red step
           call gauss_seidel_mg_fast(mg_dict,ifinelevel,safe,.false.)  ! Black step
        else
           call gauss_seidel_mg(mg_dict,ifinelevel,safe,.true. )  ! Red step
           call gauss_seidel_mg(mg_dict,ifinelevel,safe,.false.)  ! Black step
        endif
     end do
     return
  end if
  
  if(safe) then
     ncycle=ncycles_coarse_safe
  else
     ncycle=1
  endif
  
  do icycle=1,ncycle
     
     ! Pre-smoothing
     do i=1,ngs_coarse
        if(fast_solver)then
           call gauss_seidel_mg_fast(mg_dict,ifinelevel,safe,.true. )  ! Red step
           call gauss_seidel_mg_fast(mg_dict,ifinelevel,safe,.false.)  ! Black step
        else
           call gauss_seidel_mg(mg_dict,ifinelevel,safe,.true. )  ! Red step
           call gauss_seidel_mg(mg_dict,ifinelevel,safe,.false.)  ! Black step
        endif
     end do     

     ! Compute residual and restrict into upper level RHS
     if(fast_solver)then
        call cmp_residual_mg_fast(mg_dict,ifinelevel)
     else
        call cmp_residual_mg(mg_dict,ifinelevel)
     endif

     ! Restrict residual to coarser level
     call restrict_residual(ifinelevel)
     
     ! Reset correction from upper level before solve
     do igrid=head_mg(ifinelevel-1),tail_mg(ifinelevel-1)
        grid(igrid)%phi(1:twotondim)=0.0d0
     end do
     
     ! Multigrid-solve the upper level
     call recursive_multigrid(ifinelevel-1, safe)
     
     ! Interpolate coarse solution and correct back into fine solution
     call interpolate_and_correct(ifinelevel)
     
     ! Post-smoothing
     do i=1,ngs_coarse
        if(fast_solver)then
           call gauss_seidel_mg_fast(mg_dict,ifinelevel,safe,.true. )  ! Red step
           call gauss_seidel_mg_fast(mg_dict,ifinelevel,safe,.false.)  ! Black step
        else
           call gauss_seidel_mg(mg_dict,ifinelevel,safe,.true. )  ! Red step
           call gauss_seidel_mg(mg_dict,ifinelevel,safe,.false.)  ! Black step
        endif
     end do
     
  end do
  
end subroutine recursive_multigrid
#endif
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

recursive subroutine recursive_multigrid_2(r,g,m,ifinelevel, safe)
  use amr_parameters, only: dp,twotondim
  use poisson_parameters, only: ngs_fine,ngs_coarse,ncycles_coarse_safe
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer, intent(in) :: ifinelevel
  logical, intent(in) :: safe

  integer :: i, igrid, icpu, info, icycle, ncycle
  
  if(ifinelevel<=g%levelmin_mg) then
     ! Solve 'directly' :
     do i=1,2*ngs_coarse
        call gauss_seidel_mg_2(r,g,m,m%mg_dict,ifinelevel,safe,.true. )  ! Red step
        call gauss_seidel_mg_2(r,g,m,m%mg_dict,ifinelevel,safe,.false.)  ! Black step
     end do
     return
  end if
  
  if(safe) then
     ncycle=ncycles_coarse_safe
  else
     ncycle=1
  endif
  
  do icycle=1,ncycle
     
     ! Pre-smoothing
     do i=1,ngs_coarse
        call gauss_seidel_mg_2(r,g,m,m%mg_dict,ifinelevel,safe,.true. )  ! Red step
        call gauss_seidel_mg_2(r,g,m,m%mg_dict,ifinelevel,safe,.false.)  ! Black step
     end do     

     ! Compute residual and restrict into upper level RHS
     call cmp_residual_mg_2(r,g,m,m%mg_dict,ifinelevel)

     ! Restrict residual to coarser level
     call restrict_residual_2(r,g,m,ifinelevel)
     
     ! Reset correction from upper level before solve
     do igrid=m%head_mg(ifinelevel-1),m%tail_mg(ifinelevel-1)
        m%grid(igrid)%phi(1:twotondim)=0.0d0
     end do
     
     ! Multigrid-solve the upper level
     call recursive_multigrid_2(r,g,m,ifinelevel-1, safe)
     
     ! Interpolate coarse solution and correct back into fine solution
     call interpolate_and_correct_2(r,g,m,ifinelevel)
     
     ! Post-smoothing
     do i=1,ngs_coarse
        call gauss_seidel_mg_2(r,g,m,m%mg_dict,ifinelevel,safe,.true. )  ! Red step
        call gauss_seidel_mg_2(r,g,m,m%mg_dict,ifinelevel,safe,.false.)  ! Black step
     end do
     
  end do
  
end subroutine recursive_multigrid_2

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Multigrid workspace initialisation
! ------------------------------------------------------------------------
#ifdef TOTO
subroutine init_mg(ilevel)
  use amr_commons
  use poisson_commons
  use hilbert
  implicit none  
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  integer::ilevel

  integer::ilev,idom
  
  allocate(head_mg(1:nlevelmax))
  allocate(tail_mg(1:nlevelmax))
  allocate(noct_mg(1:nlevelmax))
  allocate(noct_tot_mg(1:nlevelmax))

  ! Allocate and compute multigrid Hilbert key tick marks
  allocate(domain_mg(1:ilevel))
  call domain_mg(ilevel)%copy(domain(ilevel))
  do ilev=ilevel-1,1,-1
     call domain_mg(ilev)%copy(domain_mg(ilev+1))
     do idom=0,domain_mg(ilev)%n
        domain_mg(ilev)%b(1:nhilbert,idom)=coarsen_key(domain_mg(ilev+1)%b(1:nhilbert,idom),ilev)
     end do
  end do

  ! First level of multigrid hierarchy is current AMR level
  head_mg(ilevel)=head(ilevel)
  tail_mg(ilevel)=tail(ilevel)
  noct_mg(ilevel)=noct(ilevel)
  noct_tot_mg(ilevel)=noct_tot(ilevel)
  
  ! Save first free element in AMR grid array to restore state
  ifree_mg=noct_used+1
  
end subroutine init_mg
#endif
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

subroutine init_mg_2(r,g,m,ilevel)
  use amr_parameters, only: dp,nhilbert
  use amr_commons, only: run_t,global_t,mesh_t
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel

  integer::ilev,idom
  
  allocate(m%head_mg(1:r%nlevelmax))
  allocate(m%tail_mg(1:r%nlevelmax))
  allocate(m%noct_mg(1:r%nlevelmax))
  allocate(m%noct_tot_mg(1:r%nlevelmax))

  ! Allocate and compute multigrid Hilbert key tick marks
  allocate(m%domain_mg(1:ilevel))
  call m%domain_mg(ilevel)%copy(m%domain(ilevel))
  do ilev=ilevel-1,1,-1
     call m%domain_mg(ilev)%copy(m%domain_mg(ilev+1))
     do idom=0,m%domain_mg(ilev)%n
        m%domain_mg(ilev)%b(1:nhilbert,idom)=coarsen_key(m%domain_mg(ilev+1)%b(1:nhilbert,idom),ilev)
     end do
  end do

  ! First level of multigrid hierarchy is current AMR level
  m%head_mg(ilevel)=m%head(ilevel)
  m%tail_mg(ilevel)=m%tail(ilevel)
  m%noct_mg(ilevel)=m%noct(ilevel)
  m%noct_tot_mg(ilevel)=m%noct_tot(ilevel)
  
  ! Save first free element in AMR grid array to restore state
  m%ifree_mg=m%noct_used+1
  
end subroutine init_mg_2

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ---------------------------------------------------------------------
! Coarse grid MG activation for local grids
! ---------------------------------------------------------------------
#ifdef TOTO  
subroutine build_mg(ifinelevel)
  use amr_commons
  use poisson_commons
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  
  integer,intent(in)::ifinelevel
  
  integer::icoarselevel,igrid,inbor,idim,ipos,ichild,icpu,grid_cpu,ind,info
  integer(kind=8),dimension(0:ndim)::hash_key,hash_father,hash_nbor
  integer,dimension(1:ndim)::cart_key
  integer,dimension(1:3,1:8),save::shift_oct=reshape(&
       & (/-1,-1,-1,+1,-1,-1,-1,+1,-1,+1,+1,-1,&
       &   -1,-1,+1,+1,-1,+1,-1,+1,+1,+1,+1,+1/),(/3,8/))
  integer(kind=4),dimension(1:nvector),save::dummy_state
  integer(kind=8),dimension(1:nvector,1:nhilbert),save::hk
  integer(kind=8),dimension(1:nvector,1:ndim),save::ix
  integer(kind=8),dimension(1:nhilbert),save::hks

  icoarselevel=ifinelevel-1
  ifree=noct_used+1
  head_mg(icoarselevel)=ifree
  
  hash_father(0)=icoarselevel
  
  call open_cache(operation_build_mg,domain_decompos_mg)
  
  ! Loop over fine grids
  do igrid=head_mg(ifinelevel),tail_mg(ifinelevel)
     
     hash_key(1:ndim)=grid(igrid)%ckey(1:ndim)
     
     ! Gather twotondim neighboring father grids
     do inbor=1,twotondim

#ifndef WITHOUTMPI
        ! If counter is good, check on incoming messages and perform actions
        if(mail_counter==32)then
           call check_mail(MPI_REQUEST_NULL,mg_dict)
           mail_counter=0
        endif
        mail_counter=mail_counter+1
#endif

        hash_nbor(1:ndim)=hash_key(1:ndim)+shift_oct(1:ndim,inbor)
        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ifinelevel)-1
           if(hash_nbor(idim)==ckey_max(ifinelevel))hash_nbor(idim)=0
        enddo
        hash_father(1:ndim)=hash_nbor(1:ndim)/2
        
        ! Access hash table
        ipos=hash_get(mg_dict,hash_father)

        ! If grid does not exist, create it in memory
        if(ipos==0)then
           
           ! Compute Cartesian keys of new oct
           cart_key(1:ndim)=hash_father(1:ndim)
           
           ! Compute Hilbert keys of new octs
           ix(1,1:ndim)=cart_key(1:ndim)
           call hilbert_key(ix,hk,dummy_state,0,icoarselevel-1,1)
           
           ! Check if grid sits inside processor boundaries
           hks = hk(1,1:nhilbert)
           if( domain_mg(icoarselevel)%in_rank(hks)) then
              
              ! Set grid index to a virtual grid in local main memory
              ichild=ifree
              
              ! Go to next main memory free line
              ifree=ifree+1
              if(ifree.GT.ngridmax)then
                 write(*,*)'No more free memory'
                 write(*,*)'for multigrid...'
                 write(*,*)'Increase ngridmax'
                 call clean_abort
              end if
              ! Otherwise, determine parent processor and use the cache
           else
              hks = hk(1,1:nhilbert)
              grid_cpu = domain_mg(icoarselevel)%get_rank(hks)
              
              ! If next cache line is occupied, free it.
              if(occupied(free_cache))call destage(ngridmax+free_cache,mg_dict)
              ! Set grid index to a virtual grid in local cache memory
              ichild=ngridmax+free_cache
              occupied(free_cache)=.true.
              parent_cpu(free_cache)=grid_cpu
              dirty(free_cache)=.true.
              ! Go to next free cache line
              free_cache=free_cache+1
              ncache=ncache+1
              if(free_cache.GT.ncachemax)free_cache=1
              if(ncache.GT.ncachemax)ncache=ncachemax
           endif
           
           grid(ichild)%lev=icoarselevel
           grid(ichild)%ckey(1:ndim)=cart_key(1:ndim)
           grid(ichild)%hkey(1:nhilbert)=hk(1,1:nhilbert)
           grid(ichild)%refined(1:twotondim)=.true.
           grid(ichild)%flag1(1:twotondim)=0
           grid(ichild)%flag2(1:twotondim)=0
           grid(ichild)%superoct=1

           ! Intitialize gravity variables
           do ind=1,twotondim
              grid(ichild)%f(ind,1:ndim)=0
              grid(ichild)%phi(ind)=0
              grid(ichild)%phi_old(ind)=0
           enddo

           ! Insert new grid in hash table
           call hash_set(mg_dict,hash_father,ichild)
           
        end if
     end do
     ! End loop over coarse neighbors
     
  end do
  ! End loop over grids
  
  call close_cache(mg_dict)
  
  ! Multigrid oct statistics
  tail_mg(icoarselevel)=ifree-1
  noct_mg(icoarselevel)=tail_mg(icoarselevel)-head_mg(icoarselevel)+1
  noct_used=tail_mg(icoarselevel)
  noct_tot_mg(icoarselevel)=noct_mg(icoarselevel)
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(noct_mg(icoarselevel),noct_tot_mg(icoarselevel),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif

end subroutine build_mg
#endif
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

subroutine build_mg_2(r,g,m,ifinelevel)
  use amr_parameters, only: dp,nvector,nhilbert,ndim,twotondim
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
  integer,intent(in)::ifinelevel
  
  integer::icoarselevel,igrid,inbor,idim,ipos,ichild,icpu,grid_cpu,ind,info
  integer(kind=8),dimension(0:ndim)::hash_key,hash_father,hash_nbor
  integer,dimension(1:ndim)::cart_key
  integer,dimension(1:3,1:8),save::shift_oct=reshape(&
       & (/-1,-1,-1,+1,-1,-1,-1,+1,-1,+1,+1,-1,&
       &   -1,-1,+1,+1,-1,+1,-1,+1,+1,+1,+1,+1/),(/3,8/))
  integer(kind=4),dimension(1:nvector),save::dummy_state
  integer(kind=8),dimension(1:nvector,1:nhilbert),save::hk
  integer(kind=8),dimension(1:nvector,1:ndim),save::ix
  integer(kind=8),dimension(1:nhilbert),save::hks

  icoarselevel=ifinelevel-1
  m%ifree=m%noct_used+1
  m%head_mg(icoarselevel)=m%ifree
  
  hash_father(0)=icoarselevel
  
  call open_cache_2(r,g,m,operation_build_mg,domain_decompos_mg)
  
  ! Loop over fine grids
  do igrid=m%head_mg(ifinelevel),m%tail_mg(ifinelevel)
     
     hash_key(1:ndim)=m%grid(igrid)%ckey(1:ndim)
     
     ! Gather twotondim neighboring father grids
     do inbor=1,twotondim

#ifndef WITHOUTMPI
        ! If counter is good, check on incoming messages and perform actions
        if(mail_counter==32)then
           call check_mail_2(r,g,m,MPI_REQUEST_NULL,m%mg_dict)
           mail_counter=0
        endif
        mail_counter=mail_counter+1
#endif

        hash_nbor(1:ndim)=hash_key(1:ndim)+shift_oct(1:ndim,inbor)
        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=m%ckey_max(ifinelevel)-1
           if(hash_nbor(idim)==m%ckey_max(ifinelevel))hash_nbor(idim)=0
        enddo
        hash_father(1:ndim)=hash_nbor(1:ndim)/2
        
        ! Access hash table
        ipos=hash_get(m%mg_dict,hash_father)

        ! If grid does not exist, create it in memory
        if(ipos==0)then
           
           ! Compute Cartesian keys of new oct
           cart_key(1:ndim)=hash_father(1:ndim)
           
           ! Compute Hilbert keys of new octs
           ix(1,1:ndim)=cart_key(1:ndim)
           call hilbert_key(ix,hk,dummy_state,0,icoarselevel-1,1)
           
           ! Check if grid sits inside processor boundaries
           hks = hk(1,1:nhilbert)
           if( m%domain_mg(icoarselevel)%in_rank(hks)) then
              
              ! Set grid index to a virtual grid in local main memory
              ichild=m%ifree
              
              ! Go to next main memory free line
              m%ifree=m%ifree+1
              if(m%ifree.GT.r%ngridmax)then
                 write(*,*)'No more free memory'
                 write(*,*)'for multigrid...'
                 write(*,*)'Increase ngridmax'
                 call clean_abort
              end if
              ! Otherwise, determine parent processor and use the cache
           else
              hks = hk(1,1:nhilbert)
              grid_cpu = m%domain_mg(icoarselevel)%get_rank(hks)
              
              ! If next cache line is occupied, free it.
              if(m%occupied(m%free_cache))call destage_2(r,g,m,r%ngridmax+m%free_cache,m%mg_dict)
              ! Set grid index to a virtual grid in local cache memory
              ichild=r%ngridmax+m%free_cache
              m%occupied(m%free_cache)=.true.
              m%parent_cpu(m%free_cache)=grid_cpu
              m%dirty(m%free_cache)=.true.
              ! Go to next free cache line
              m%free_cache=m%free_cache+1
              m%ncache=m%ncache+1
              if(m%free_cache.GT.r%ncachemax)m%free_cache=1
              if(m%ncache.GT.r%ncachemax)m%ncache=r%ncachemax
           endif
           
           m%grid(ichild)%lev=icoarselevel
           m%grid(ichild)%ckey(1:ndim)=cart_key(1:ndim)
           m%grid(ichild)%hkey(1:nhilbert)=hk(1,1:nhilbert)
           m%grid(ichild)%refined(1:twotondim)=.true.
           m%grid(ichild)%flag1(1:twotondim)=0
           m%grid(ichild)%flag2(1:twotondim)=0
           m%grid(ichild)%superoct=1

           ! Intitialize gravity variables
           do ind=1,twotondim
              m%grid(ichild)%f(ind,1:ndim)=0
              m%grid(ichild)%phi(ind)=0
              m%grid(ichild)%phi_old(ind)=0
           enddo

           ! Insert new grid in hash table
           call hash_set(m%mg_dict,hash_father,ichild)
           
        end if
     end do
     ! End loop over coarse neighbors
     
  end do
  ! End loop over grids
  
  call close_cache_2(r,g,m,m%mg_dict)
  
  ! Multigrid oct statistics
  m%tail_mg(icoarselevel)=m%ifree-1
  m%noct_mg(icoarselevel)=m%tail_mg(icoarselevel)-m%head_mg(icoarselevel)+1
  m%noct_used=m%tail_mg(icoarselevel)
  m%noct_tot_mg(icoarselevel)=m%noct_mg(icoarselevel)
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(m%noct_mg(icoarselevel),m%noct_tot_mg(icoarselevel),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif

end subroutine build_mg_2

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Multigrid cleanup
! ------------------------------------------------------------------------
#ifdef TOTO
subroutine cleanup_mg
   use amr_commons
   use poisson_commons
   implicit none
   integer :: ilev

   ! Deallocate processor boundary array
   deallocate(head_mg,tail_mg,noct_mg,noct_tot_mg)
   do ilev=1,size(domain_mg)
      call domain_mg(ilev)%destroy
   end do
   deallocate(domain_mg)

  ! Reset the MG hash table
  call reset_entire_hash(mg_dict,.false.)
  
  ! Restore AMR grid array into its original state
  ifree=ifree_mg
  noct_used=ifree-1

end subroutine cleanup_mg
#endif
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

subroutine cleanup_mg_2(r,g,m)
  use amr_commons, only: run_t,global_t,mesh_t
  use hash
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m

  integer :: ilev

   ! Deallocate processor boundary array
   deallocate(m%head_mg,m%tail_mg,m%noct_mg,m%noct_tot_mg)
   do ilev=1,size(m%domain_mg)
      call m%domain_mg(ilev)%destroy
   end do
   deallocate(m%domain_mg)

  ! Reset the MG hash table
  call reset_entire_hash(m%mg_dict,.false.)
  
  ! Restore AMR grid array into its original state
  m%ifree=m%ifree_mg
  m%noct_used=m%ifree-1

end subroutine cleanup_mg_2

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Initialize mask at fine level into f(:,3)
! ------------------------------------------------------------------------
#ifdef TOTO
subroutine make_mask(ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer, intent(in) :: ilevel
  
  integer  :: igrid, ind
  
  ! Init mask to 1.0 on all fine level cells :
  do igrid=head(ilevel),tail(ilevel)
     do ind=1,twotondim
        grid(igrid)%f(ind,3)=1.0d0
     end do
  end do
  
end subroutine make_mask
#endif
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

subroutine make_mask_2(r,g,m,ilevel)
  use amr_parameters, only: twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer, intent(in) :: ilevel

  integer  :: igrid, ind
  
  ! Init mask to 1.0 on all fine level cells :
  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim
        m%grid(igrid)%f(ind,3)=1.0d0
     end do
  end do
  
end subroutine make_mask_2

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Preprocess the fine (AMR) level RHS to account for boundary conditions
!
!  _____#_____
! |     #     |      Cell I is INSIDE active domain (mask > 0)
! |  I  #  O  |      Cell O is OUTSIDE (mask <= 0 or nonexistent cell)
! |_____#_____|      # is the boundary
!       #
!
! phi(I) and phi(O) must BOTH be set at call time, if applicable
! phi(#) is computed from phi(I), phi(O) and the mask values
! If AMR cell O does not exist, phi(O) is computed by interpolation
!
! Sets BC-modified RHS    into f(:,2)
!
! ------------------------------------------------------------------------
#ifdef TOTO
subroutine make_bc_rhs(ilevel,icount)

  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer, intent(in) :: ilevel,icount
  
  integer, dimension(1:3,1:2,1:8) :: iii, jjj
  integer::igrid,idim,ind,igridn,inbor,ig,id
  integer::get_grid
  integer,dimension(1:8,1:8)::ccc
  real(dp)::aa,bb,cc,dd,tfrac
  real(dp),dimension(1:8)::bbb
  integer(kind=8),dimension(0:ndim)::hash_key,hash_nbor
  integer,dimension(1:threetondim),save::igrid_nbor,ind_nbor
  real(dp),dimension(1:twotondim),save::phi_int
  real(dp),dimension(1:twotondim,0:twondim),save::phi_nbor,dis_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))

  real(dp) :: dx, oneoverdx2, phi_b, nb_mask, nb_phi, w
  real(dp) :: scale, fourpi
  
  ! Set constants
  fourpi = 4.D0*ACOS(-1.0D0)
  if(cosmo) fourpi = 1.5D0*omega_m*aexp
  
  dx  = boxlen/2**ilevel
  oneoverdx2 = 1.0d0/(dx*dx)
  
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
        dis_nbor(ind,0)=grid(igrid)%f(ind,3)
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
              phi_nbor(ind,inbor)=grid(igridn)%phi(ind)
              dis_nbor(ind,inbor)=grid(igridn)%f(ind,3)
           end do

        ! Otherwise interpolate from coarser level
        else

           ! Get 3**ndim neighbouring parent cells using read-only cache
           call get_threetondim_nbor_parent_cell(hash_nbor,grid_dict,igrid_nbor,ind_nbor,.false.,.true.)
           call interpol_phi(igrid_nbor,ind_nbor,ccc,bbb,tfrac,phi_nbor(1,inbor))
           do ind=1,threetondim
              call unlock_cache(igrid_nbor(ind))
           end do
           do ind=1,twotondim
              dis_nbor(ind,inbor)=-1.0
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim
        
        ! Init BC-modified RHS to rho - rho_tot :
        grid(igrid)%f(ind,2) = fourpi*(grid(igrid)%rho(ind) - rho_tot)
        
        ! Do not process masked cells
        if(grid(igrid)%f(ind,3)<=0.0) cycle 
        
        ! Separate directions
        do idim=1,ndim

           ! Loop over the 2 neighbors
           do inbor=1,2
              id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
              
              nb_mask=dis_nbor(id,ig)
              if(nb_mask>0.0)cycle

              ! phi(#) interpolated with mask:
              nb_phi = phi_nbor(id,ig)
              w = nb_mask/(nb_mask-grid(igrid)%f(ind,3)) ! Linear parameter
              phi_b = ((1.0d0-w)*nb_phi + w*grid(igrid)%phi(ind))

              ! Increment correction for current cell
              grid(igrid)%f(ind,2) = grid(igrid)%f(ind,2) - 2.0d0*oneoverdx2*phi_b

           end do
        end do

     end do
  end do

  call close_cache(grid_dict)

end subroutine make_bc_rhs
#endif
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

subroutine make_bc_rhs_2(r,g,m,ilevel,icount)
  use amr_parameters, only: dp,ndim,twondim,twotondim,threetondim
  use amr_commons, only: run_t,global_t,mesh_t
  use cache_commons
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer, intent(in) :: ilevel,icount
  
  integer, dimension(1:3,1:2,1:8) :: iii, jjj
  integer::igrid,idim,ind,igridn,inbor,ig,id
  integer::get_grid_2
  integer,dimension(1:8,1:8)::ccc
  real(dp)::aa,bb,cc,dd,tfrac
  real(dp),dimension(1:8)::bbb
  integer(kind=8),dimension(0:ndim)::hash_key,hash_nbor
  integer,dimension(1:threetondim),save::igrid_nbor,ind_nbor
  real(dp),dimension(1:twotondim),save::phi_int
  real(dp),dimension(1:twotondim,0:twondim),save::phi_nbor,dis_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))

  real(dp) :: dx, oneoverdx2, phi_b, nb_mask, nb_phi, w
  real(dp) :: scale, fourpi
  
  ! Set constants
  fourpi = 4.D0*ACOS(-1.0D0)
  if(r%cosmo) fourpi = 1.5D0*g%omega_m*g%aexp
  
  dx  = r%boxlen/2**ilevel
  oneoverdx2 = 1.0d0/(dx*dx)
  
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
     call clean_stop(g)
  endif

  ! Compute fraction of time steps for interpolation
  if (g%dtold(ilevel-1)>0.0)then
     tfrac=g%dtnew(ilevel)/g%dtold(ilevel-1)*(icount-1)
  else
     tfrac=0.0
  end if

  call open_cache_2(r,g,m,operation_interpol,domain_decompos_amr)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)
     
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
        igridn=get_grid_2(r,g,m,hash_nbor,m%grid_dict,.false.,.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=m%grid(igridn)%phi(ind)
              dis_nbor(ind,inbor)=m%grid(igridn)%f(ind,3)
           end do

        ! Otherwise interpolate from coarser level
        else

           ! Get 3**ndim neighbouring parent cells using read-only cache
           call get_threetondim_nbor_parent_cell_2(r,g,m,hash_nbor,m%grid_dict,igrid_nbor,ind_nbor,.false.,.true.)
           call interpol_phi_2(r,g,m,igrid_nbor,ind_nbor,ccc,bbb,tfrac,phi_nbor(1,inbor))
           do ind=1,threetondim
              call unlock_cache_2(r,g,m,igrid_nbor(ind))
           end do
           do ind=1,twotondim
              dis_nbor(ind,inbor)=-1.0
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim
        
        ! Init BC-modified RHS to rho - rho_tot :
        m%grid(igrid)%f(ind,2) = fourpi*(m%grid(igrid)%rho(ind) - g%rho_tot)
        
        ! Do not process masked cells
        if(m%grid(igrid)%f(ind,3)<=0.0) cycle 
        
        ! Separate directions
        do idim=1,ndim

           ! Loop over the 2 neighbors
           do inbor=1,2
              id=jjj(idim,inbor,ind); ig=iii(idim,inbor,ind)
              
              nb_mask=dis_nbor(id,ig)
              if(nb_mask>0.0)cycle

              ! phi(#) interpolated with mask:
              nb_phi = phi_nbor(id,ig)
              w = nb_mask/(nb_mask-m%grid(igrid)%f(ind,3)) ! Linear parameter
              phi_b = ((1.0d0-w)*nb_phi + w*m%grid(igrid)%phi(ind))

              ! Increment correction for current cell
              m%grid(igrid)%f(ind,2) = m%grid(igrid)%f(ind,2) - 2.0d0*oneoverdx2*phi_b

           end do
        end do

     end do
  end do

  call close_cache_2(r,g,m,m%grid_dict)

end subroutine make_bc_rhs_2

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
#ifdef TOTO
subroutine build_comm_mg(hash_dict,ilevel)
  use amr_commons
  use poisson_commons
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  !
  integer,intent(in)::ilevel
  type(hash_table)::hash_dict
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
  integer(kind=8),dimension(1:nhilbert),save::hks
  integer(kind=8),dimension(1:nvector,1:ndim),save::ix
  integer::istart,nbuffer,countrecv,countsend,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv

#ifndef WITHOUTMPI

  allocate(nremote(1:ncpu))  
  allocate(buffer_mg(ilevel)%nbor_indx(head_mg(ilevel):tail_mg(ilevel),1:twondim))

  hash_nbor(0)=ilevel
  
  nremote=0
  call open_cache(operation_mg,domain_decompos_mg)
  
  ! Loop over grids
  do igrid=head_mg(ilevel),tail_mg(ilevel)
     
     ! Gather twondim neighboring grids
     do inbor=1,twondim

        hash_nbor(1:ndim)=grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
           if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
        enddo
        
        ! Get neighbouring grid using read-only cache
        igridn=get_grid(hash_nbor,hash_dict,.false.,.true.)
        
        ! If grid exists, determine if remote
        if(igridn<=ngridmax)then
           buffer_mg(ilevel)%nbor_indx(igrid,inbor)=igridn
        else
           ! Compute Cartesian keys of new oct
           cart_key(1:ndim)=hash_nbor(1:ndim)
           
           ! Compute Hilbert keys of new octs
           ix(1,1:ndim)=cart_key(1:ndim)
           call hilbert_key(ix,hk,dummy_state,0,ilevel-1,1)

           ! Determine parent processor and increment counter
           hks = hk(1,1:nhilbert)
           grid_cpu = domain_mg(ilevel)%get_rank(hks)
           nremote(grid_cpu)=nremote(grid_cpu)+1
        end if

     end do
     ! End loop over neighbors
     
  end do
  ! End loop over grids
  
  call close_cache(hash_dict)

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
  allocate(buffer_mg(ilevel)%send_cnt(1:ncpu))
  allocate(buffer_mg(ilevel)%send_oft(1:ncpu))
  allocate(buffer_mg(ilevel)%recv_cnt(1:ncpu))
  allocate(buffer_mg(ilevel)%recv_oft(1:ncpu))

  buffer_mg(ilevel)%send_cnt=0; buffer_mg(ilevel)%send_oft=0; buffer_mg(ilevel)%send_tot=0
  buffer_mg(ilevel)%recv_cnt=0; buffer_mg(ilevel)%recv_oft=0; buffer_mg(ilevel)%recv_tot=0
  do icpu=1,ncpu
     buffer_mg(ilevel)%send_cnt(icpu)=nalltoall(myid,icpu)
     buffer_mg(ilevel)%recv_cnt(icpu)=nalltoall(icpu,myid)
     buffer_mg(ilevel)%send_tot=buffer_mg(ilevel)%send_tot+buffer_mg(ilevel)%send_cnt(icpu)
     buffer_mg(ilevel)%recv_tot=buffer_mg(ilevel)%recv_tot+buffer_mg(ilevel)%recv_cnt(icpu)
     if(icpu<ncpu)then
        buffer_mg(ilevel)%send_oft(icpu+1)=buffer_mg(ilevel)%send_oft(icpu)+nalltoall(myid,icpu)
        buffer_mg(ilevel)%recv_oft(icpu+1)=buffer_mg(ilevel)%recv_oft(icpu)+nalltoall(icpu,myid)
     endif
  end do
  deallocate(nalltoall)

#if NDIM>0
  allocate(x_send_buf(1:buffer_mg(ilevel)%send_tot))
#endif
#if NDIM>1
  allocate(y_send_buf(1:buffer_mg(ilevel)%send_tot))
#endif
#if NDIM>2
  allocate(z_send_buf(1:buffer_mg(ilevel)%send_tot))
#endif

  hash_nbor(0)=ilevel

  nremote=0
  call open_cache(operation_mg,domain_decompos_mg)

  ! Loop over grids
  do igrid=head_mg(ilevel),tail_mg(ilevel)

     ! Gather twondim neighboring grids
     do inbor=1,twondim

        hash_nbor(1:ndim)=grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
           if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
        enddo

        ! Get neighbouring grid using read-only cache
        igridn=get_grid(hash_nbor,hash_dict,.false.,.true.)

        ! If grid is remote
        if(igridn>ngridmax)then

           ! Compute Cartesian keys of new oct
           cart_key(1:ndim)=hash_nbor(1:ndim)

           ! Compute Hilbert keys of new octs
           ix(1,1:ndim)=cart_key(1:ndim)
           call hilbert_key(ix,hk,dummy_state,0,ilevel-1,1)

           ! Determine parent processor and increment counter
           hks = hk(1,1:nhilbert)
           grid_cpu = domain_mg(ilevel)%get_rank(hks)
           nremote(grid_cpu)=nremote(grid_cpu)+1
           iremote=buffer_mg(ilevel)%send_oft(grid_cpu)+nremote(grid_cpu)
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
           buffer_mg(ilevel)%nbor_indx(igrid,inbor)=-iremote
        end if

     end do
     ! End loop over neighbors

  end do
  ! End loop over grids

  call close_cache(hash_dict)

  deallocate(nremote)

#if NDIM>0
  allocate(x_recv_buf(1:buffer_mg(ilevel)%recv_tot))
  call MPI_ALLTOALLV(x_send_buf,buffer_mg(ilevel)%send_cnt,buffer_mg(ilevel)%send_oft,MPI_INTEGER, &
       &             x_recv_buf,buffer_mg(ilevel)%recv_cnt,buffer_mg(ilevel)%recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  deallocate(x_send_buf)
#endif
#if NDIM>1
  allocate(y_recv_buf(1:buffer_mg(ilevel)%recv_tot))
  call MPI_ALLTOALLV(y_send_buf,buffer_mg(ilevel)%send_cnt,buffer_mg(ilevel)%send_oft,MPI_INTEGER, &
       &             y_recv_buf,buffer_mg(ilevel)%recv_cnt,buffer_mg(ilevel)%recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  deallocate(y_send_buf)
#endif
#if NDIM>2
  allocate(z_recv_buf(1:buffer_mg(ilevel)%recv_tot))
  call MPI_ALLTOALLV(z_send_buf,buffer_mg(ilevel)%send_cnt,buffer_mg(ilevel)%send_oft,MPI_INTEGER, &
       &             z_recv_buf,buffer_mg(ilevel)%recv_cnt,buffer_mg(ilevel)%recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  deallocate(z_send_buf)
#endif

  allocate(buffer_mg(ilevel)%grid_recv_buf(1:buffer_mg(ilevel)%recv_tot))

  hash_key(0)=ilevel
  do i=1,buffer_mg(ilevel)%recv_tot
     hash_key(1)=x_recv_buf(i)
     hash_key(2)=y_recv_buf(i)
     hash_key(3)=z_recv_buf(i)
     ipos=hash_get(hash_dict,hash_key)
     buffer_mg(ilevel)%grid_recv_buf(i)=ipos
  end do
  deallocate(x_recv_buf)
  deallocate(y_recv_buf)
  deallocate(z_recv_buf)

  allocate(buffer_mg(ilevel)%phi_send_buf(1:buffer_mg(ilevel)%send_tot*twotondim))
  allocate(buffer_mg(ilevel)%phi_recv_buf(1:buffer_mg(ilevel)%recv_tot*twotondim))
  allocate(buffer_mg(ilevel)%phi_remote(1:buffer_mg(ilevel)%send_tot,1:twotondim))
  allocate(buffer_mg(ilevel)%dis_remote(1:buffer_mg(ilevel)%send_tot,1:twotondim))

  ! Update boundary conditions for distance function
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
        buffer_mg(ilevel)%phi_recv_buf(istart)=grid(igrid)%f(ind,3)
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
        buffer_mg(ilevel)%dis_remote(i,ind)=buffer_mg(ilevel)%phi_send_buf(istart)
     end do
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

end subroutine build_comm_mg
#endif
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
#ifdef TOTO
subroutine clean_comm_mg(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ilevel

#ifndef WITHOUTMPI
  
  deallocate(buffer_mg(ilevel)%nbor_indx)
  deallocate(buffer_mg(ilevel)%phi_remote)
  deallocate(buffer_mg(ilevel)%dis_remote)
  deallocate(buffer_mg(ilevel)%phi_send_buf)
  deallocate(buffer_mg(ilevel)%phi_recv_buf)
  deallocate(buffer_mg(ilevel)%grid_recv_buf)
  deallocate(buffer_mg(ilevel)%send_cnt)
  deallocate(buffer_mg(ilevel)%send_oft)
  deallocate(buffer_mg(ilevel)%recv_cnt)
  deallocate(buffer_mg(ilevel)%recv_oft)

#endif

end subroutine clean_comm_mg
#endif
#endif

