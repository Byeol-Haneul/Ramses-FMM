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
#ifdef GRAV
! ------------------------------------------------------------------------
! Main multigrid routine, called by amr_step
! ------------------------------------------------------------------------

subroutine multigrid(r,g,m,p,mdl,ilevel,icount)
  use amr_parameters, only: dp,twotondim
  use poisson_parameters, only: ngs_fine,ngs_coarse,ncycles_coarse_safe
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer, intent(in) :: ilevel,icount
  
  integer, parameter  :: MAXITER  = 10
  real(dp), parameter :: SAFE_FACTOR = 0.5
  
  integer  :: igrid, ifine, i, iter, allmasked
  integer,dimension(1:4) :: input_array
  integer,dimension(1:4) :: output_array
  real(kind=8) :: res_norm2, i_res_norm2
  real(kind=8) :: err, last_err
  real(kind=8) :: i_res_norm2_tot, res_norm2_tot
  
  if(r%gravity_type>0)return
  if(m%noct_tot(ilevel)==0)return
  
  if(r%verbose) print '(A,I2)','Entering multigrid at level ',ilevel
  
  ! ---------------------------------------------------------------------
  ! Prepare first guess, mask and BCs at finest level
  ! ---------------------------------------------------------------------
  input_array(1)=ilevel
  input_array(2)=icount
  call r_make_initial_phi(r,g,m,p,mdl,mdl%ncpu,2,0,input_array)  ! Initial guess
  call r_make_mask(r,g,m,p,mdl,mdl%ncpu,1,0,ilevel)              ! Fill the fine level mask
  call r_make_bc_rhs(r,g,m,p,mdl,mdl%ncpu,2,0,input_array)       ! Fill BC-modified RHS

  ! ---------------------------------------------------------------------
  ! Initialize Domain Decomposition and Hash Table for Multigrid
  ! ---------------------------------------------------------------------
  call r_init_mg(r,g,m,p,mdl,mdl%ncpu,1,0,ilevel)
  
  ! ---------------------------------------------------------------------
  ! Build Multigrid hierarchy in memory
  ! ---------------------------------------------------------------------
  do ifine=ilevel,2,-1
     call r_build_mg(r,g,m,p,mdl,mdl%ncpu,1,0,ifine)
  end do
  
  ! ---------------------------------------------------------------------
  ! Restrict mask up
  ! ---------------------------------------------------------------------
  g%levelmin_mg=1
  do ifine=ilevel,2,-1
     ! Restrict and communicate mask
     call r_restrict_mask(r,g,m,p,mdl,mdl%ncpu,1,1,ifine,allmasked)
     if(allmasked==1) then ! Coarser level is fully masked: stop here
        g%levelmin_mg=ifine
        exit
     end if
  end do
  
  ! ---------------------------------------------------------------------
  ! Set scan flag (for optimisation)
  ! ---------------------------------------------------------------------
  input_array(1)=ilevel
  do ifine=ilevel,g%levelmin_mg,-1
     input_array(2)=ifine
     call r_set_scan_flag(r,g,m,p,mdl,mdl%ncpu,2,0,input_array)
  end do
  
  ! ---------------------------------------------------------------------
  ! Initiate solve at fine level
  ! ---------------------------------------------------------------------
  
  iter = 0
  err = 1.0d0
  main_iteration_loop: do

     iter=iter+1

     input_array(1)=ilevel
     input_array(2)=ilevel
     if(g%safe_mode(ilevel))then
        input_array(3)=1
     else
        input_array(3)=0
     endif
     
     ! Pre-smoothing
     do i=1,ngs_fine
        input_array(4)=1  ! Red step
        call r_gauss_seidel_mg(r,g,m,p,mdl,mdl%ncpu,4,0,input_array)
        input_array(4)=0  ! Black step
        call r_gauss_seidel_mg(r,g,m,p,mdl,mdl%ncpu,4,0,input_array)
     end do
     
     ! Compute new residual
     call r_cmp_residual_mg(r,g,m,p,mdl,mdl%ncpu,2,0,input_array)

     ! Compute initial residual norm
     if(iter==1) then
        call r_cmp_residual_norm2(r,g,m,p,mdl,mdl%ncpu,1,2,ilevel,output_array)
        i_res_norm2=transfer(output_array(1:2),i_res_norm2)
     end if

     if(ilevel>1) then

        ! Restrict residual to coarser level
        call r_restrict_residual(r,g,m,p,mdl,mdl%ncpu,1,0,ilevel)

        ! Reset correction from upper level before solve
        call r_reset_correction(r,g,m,p,mdl,mdl%ncpu,1,0,ilevel-1)
        
        ! Multigrid-solve the upper level
        call recursive_multigrid(r,g,m,p,mdl,ilevel-1, g%safe_mode(ilevel))
        
        ! Interpolate coarse solution and correct fine solution
        call r_interpolate_and_correct(r,g,m,p,mdl,mdl%ncpu,1,0,ilevel)

     end if

     ! Post-smoothing
     do i=1,ngs_fine
        input_array(4)=1  ! Red step
        call r_gauss_seidel_mg(r,g,m,p,mdl,mdl%ncpu,4,0,input_array)
        input_array(4)=0  ! Black step
        call r_gauss_seidel_mg(r,g,m,p,mdl,mdl%ncpu,4,0,input_array)
     end do
     
     ! Update fine residual
     call r_cmp_residual_mg(r,g,m,p,mdl,mdl%ncpu,2,0,input_array)

     ! Compute residual norm
     call r_cmp_residual_norm2(r,g,m,p,mdl,mdl%ncpu,1,2,ilevel,output_array)
     res_norm2=transfer(output_array(1:2),res_norm2)
     
     last_err = err
     err = sqrt(res_norm2/(i_res_norm2+1d-20*g%rho_tot**2))
     
     ! Verbosity
     if(r%verbose) print '(A,I5,A,1pE10.3)','   ==> Step=',iter,' Error=',err
     
     ! Converged?
     if(err<r%epsilon .or. iter>=MAXITER) exit
     
     ! Not converged, check error and possibly enable safe mode for the level
     if(err > last_err*SAFE_FACTOR .and. (.not. g%safe_mode(ilevel))) then
        if(r%verbose)print *,'CAUTION: Switching to safe MG mode for level ',ilevel
        g%safe_mode(ilevel) = .true.
     end if
     
  end do main_iteration_loop
  
  print '(A,I5,A,I5,A,1pE10.3)','   ==> Level=',ilevel,' Step=',iter,' Error=',err
  if(iter==MAXITER) print *,'WARN: Fine multigrid Poisson failed to converge...'
    
  ! ---------------------------------------------------------------------
  ! Cleanup MG levels after solve complete
  ! ---------------------------------------------------------------------
  call r_cleanup_mg(r,g,m,p,mdl,mdl%ncpu,0,0)
  
end subroutine multigrid

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Recursive multigrid routine for coarse MG levels
! ------------------------------------------------------------------------

recursive subroutine recursive_multigrid(r,g,m,p,mdl,ifinelevel, safe)
  use amr_parameters, only: dp,twotondim
  use poisson_parameters, only: ngs_fine,ngs_coarse,ncycles_coarse_safe
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer, intent(in) :: ifinelevel
  logical, intent(in) :: safe

  integer :: i, igrid, icycle, ncycle
  integer,dimension(1:4) :: input_array
  
  ! Set parameter array
  input_array(1)=ifinelevel+1
  input_array(2)=ifinelevel
  if(safe)then
     input_array(3)=1
  else
     input_array(3)=0
  endif
  
  if(ifinelevel<=g%levelmin_mg) then
     ! Solve 'directly' :
     do i=1,2*ngs_coarse
        input_array(4)=1  ! Red step
        call r_gauss_seidel_mg(r,g,m,p,mdl,mdl%ncpu,4,0,input_array)
        input_array(4)=0  ! Black step
        call r_gauss_seidel_mg(r,g,m,p,mdl,mdl%ncpu,4,0,input_array)
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
        input_array(4)=1  ! Red step
        call r_gauss_seidel_mg(r,g,m,p,mdl,mdl%ncpu,4,0,input_array)
        input_array(4)=0  ! Black step
        call r_gauss_seidel_mg(r,g,m,p,mdl,mdl%ncpu,4,0,input_array)
     end do     

     ! Compute residual and restrict into upper level RHS
     call r_cmp_residual_mg(r,g,m,p,mdl,mdl%ncpu,2,0,input_array)

     ! Restrict residual to coarser level
     call r_restrict_residual(r,g,m,p,mdl,mdl%ncpu,1,0,ifinelevel)
     
     ! Reset correction from upper level before solve
     call r_reset_correction(r,g,m,p,mdl,mdl%ncpu,1,0,ifinelevel-1)
     
     ! Multigrid-solve the upper level
     call recursive_multigrid(r,g,m,p,mdl,ifinelevel-1, safe)
     
     ! Interpolate coarse solution and correct back into fine solution
     call r_interpolate_and_correct(r,g,m,p,mdl,mdl%ncpu,1,0,ifinelevel)
     
     ! Post-smoothing
     do i=1,ngs_coarse
        input_array(4)=1  ! Red step
        call r_gauss_seidel_mg(r,g,m,p,mdl,mdl%ncpu,4,0,input_array)
        input_array(4)=0  ! Black step
        call r_gauss_seidel_mg(r,g,m,p,mdl,mdl%ncpu,4,0,input_array)
     end do
     
  end do
  
end subroutine recursive_multigrid

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Multigrid workspace initialisation
! ------------------------------------------------------------------------

recursive subroutine r_init_mg(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size
  integer::ilevel

  integer::next_range,next_cpu

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_INIT_MG,next_cpu,next_range,input_size,output_size,ilevel)
     call r_init_mg(r,g,m,p,mdl,next_range,input_size,output_size,ilevel)
  else
     call init_mg(r,m,ilevel)
  endif

end subroutine r_init_mg

subroutine init_mg(r,m,ilevel)
  use amr_parameters, only: dp,nhilbert
  use amr_commons, only: run_t,mesh_t
  use hilbert
  implicit none
  type(run_t)::r
  type(mesh_t)::m
  integer::ilevel

  integer::ilev,idom
  
  allocate(m%head_mg(1:r%nlevelmax))
  allocate(m%tail_mg(1:r%nlevelmax))
  allocate(m%noct_mg(1:r%nlevelmax))

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
  
  ! Save first free element in AMR grid array to restore state
  m%ifree_mg=m%noct_used+1
  
end subroutine init_mg

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ---------------------------------------------------------------------
! Coarse grid MG activation for local grids
! ---------------------------------------------------------------------

recursive subroutine r_build_mg(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size
  integer::ilevel

  integer::next_range,next_cpu

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_BUILD_MG,next_cpu,next_range,input_size,output_size,ilevel)
     call r_build_mg(r,g,m,p,mdl,next_range,input_size,output_size,ilevel)
  else
     call build_mg(r,g,m,ilevel)
  endif

end subroutine r_build_mg

subroutine build_mg(r,g,m,ifinelevel)
  use amr_parameters, only: dp,nhilbert,ndim,twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  use cache_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer,intent(in)::ifinelevel
  
  integer::icoarselevel,igrid,inbor,idim,ipos,ichild,grid_cpu,ind
  integer(kind=8),dimension(0:ndim)::hash_key,hash_father,hash_nbor
  integer(kind=4),dimension(1:ndim)::cart_key
  integer,dimension(1:3,1:8),save::shift_oct=reshape(&
       & (/-1,-1,-1,+1,-1,-1,-1,+1,-1,+1,+1,-1,&
       &   -1,-1,+1,+1,-1,+1,-1,+1,+1,+1,+1,+1/),(/3,8/))
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix

  icoarselevel=ifinelevel-1
  m%ifree=m%noct_used+1
  m%head_mg(icoarselevel)=m%ifree
  
  hash_father(0)=icoarselevel
  
  call open_cache(r,g,m,operation_build_mg,domain_decompos_mg)
  
  ! Loop over fine grids
  do igrid=m%head_mg(ifinelevel),m%tail_mg(ifinelevel)
     
     hash_key(1:ndim)=m%grid(igrid)%ckey(1:ndim)
     
     ! Gather twotondim neighboring father grids
     do inbor=1,twotondim

#ifndef WITHOUTMPI
        ! If counter is good, check on incoming messages and perform actions
        if(mail_counter==32)then
           call check_mail(r,g,m,MPI_REQUEST_NULL,m%mg_dict)
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
           cart_key(1:ndim)=int(hash_father(1:ndim),kind=4)
           
           ! Compute Hilbert keys of new octs
           ix(1:ndim)=cart_key(1:ndim)
           hk(1:nhilbert)=hilbert_key(ix,icoarselevel-1)
           
           ! Check if grid sits inside processor boundaries
           if( m%domain_mg(icoarselevel)%in_rank(hk)) then
              
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

           else
              ! Otherwise, determine parent processor and use the cache
              grid_cpu = m%domain_mg(icoarselevel)%get_rank(hk)
              
              ! If next cache line is occupied, free it.
              if(m%occupied(m%free_cache))call destage(r,g,m,r%ngridmax+m%free_cache,m%mg_dict)
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
           m%grid(ichild)%hkey(1:nhilbert)=hk(1:nhilbert)
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
  
  call close_cache(r,g,m,m%mg_dict)
  
  ! Multigrid oct statistics
  m%tail_mg(icoarselevel)=m%ifree-1
  m%noct_mg(icoarselevel)=m%tail_mg(icoarselevel)-m%head_mg(icoarselevel)+1
  m%noct_used=m%tail_mg(icoarselevel)

end subroutine build_mg

subroutine pack_flush_build_mg(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size)::msg_array

  integer::ind
  type(msg_small_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp(ind)=0.0d0
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_build_mg

subroutine unpack_flush_build_mg(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size)::msg_array

  integer::idim,ind
  type(msg_small_realdp)::msg

  msg=transfer(msg_array,msg)
  
  do ind=1,twotondim
     grid%refined(ind)=.false.
  end do
  
#ifdef GRAV
  do idim=1,ndim
     do ind=1,twotondim
        grid%f(ind,idim)=0.0d0
     end do
  end do
  do ind=1,twotondim
     grid%phi(ind)=0.0d0
     grid%phi_old(ind)=0.0d0
  end do
#endif

end subroutine unpack_flush_build_mg

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Multigrid cleanup
! ------------------------------------------------------------------------

recursive subroutine r_cleanup_mg(r,g,m,p,mdl,cpu_range,input_size,output_size)
  use amr_parameters, only: twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size

  integer::next_range,next_cpu
  integer::igrid
  
  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_CLEANUP_MG,next_cpu,next_range,input_size,output_size)
     call r_cleanup_mg(r,g,m,p,mdl,next_range,input_size,output_size)
  else
     call cleanup_mg(m)
  endif

end subroutine r_cleanup_mg

subroutine cleanup_mg(m)
  use amr_commons, only: mesh_t
  use hash
  implicit none
  type(mesh_t)::m

  integer :: ilev

   ! Deallocate processor boundary array
   deallocate(m%head_mg,m%tail_mg,m%noct_mg)
   do ilev=1,size(m%domain_mg)
      call m%domain_mg(ilev)%destroy
   end do
   deallocate(m%domain_mg)

  ! Reset the MG hash table
  call reset_entire_hash(m%mg_dict,.false.)
  
  ! Restore AMR grid array into its original state
  m%ifree=m%ifree_mg
  m%noct_used=m%ifree-1

end subroutine cleanup_mg

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Initialize mask at fine level into f(:,3)
! ------------------------------------------------------------------------

recursive subroutine r_make_mask(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size
  integer::ilevel

  integer::next_range,next_cpu

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_MAKE_MASK,next_cpu,next_range,input_size,output_size,ilevel)
     call r_make_mask(r,g,m,p,mdl,next_range,input_size,output_size,ilevel)
  else
     call make_mask(m,ilevel)
  endif

end subroutine r_make_mask

subroutine make_mask(m,ilevel)
  use amr_parameters, only: twotondim
  use amr_commons, only: mesh_t
  implicit none
  type(mesh_t)::m
  integer, intent(in) :: ilevel

  integer  :: igrid, ind
  
  ! Init mask to 1.0 on all fine level cells :
  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim
        m%grid(igrid)%f(ind,3)=1.0d0
     end do
  end do
  
end subroutine make_mask

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

recursive subroutine r_make_bc_rhs(r,g,m,p,mdl,cpu_range,input_size,output_size,input_array)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size
  integer,dimension(1:input_size)::input_array

  integer::next_range,next_cpu
  integer::ilevel,icount
  
  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_MAKE_BC_RHS,next_cpu,next_range,input_size,output_size,input_array)
     call r_make_bc_rhs(r,g,m,p,mdl,next_range,input_size,output_size,input_array)
  else
     ilevel=input_array(1)
     icount=input_array(2)
     call make_bc_rhs(r,g,m,ilevel,icount)
  endif

end subroutine r_make_bc_rhs

subroutine make_bc_rhs(r,g,m,ilevel,icount)
  use amr_parameters, only: dp,ndim,twondim,twotondim,threetondim
  use amr_commons, only: run_t,global_t,mesh_t
  use cache_commons
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m

  integer, intent(in) :: ilevel,icount
  
  integer, dimension(1:3,1:2,1:8) :: iii, jjj
  integer::igrid,idim,ind,igridn,inbor,ig,id
  integer::get_grid
  integer,dimension(1:8,1:8)::ccc
  real(dp)::aa,bb,cc,dd,tfrac
  real(dp),dimension(1:8)::bbb
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer,dimension(1:threetondim),save::igrid_nbor,ind_nbor
  real(dp),dimension(1:twotondim,0:twondim),save::phi_nbor,dis_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))

  real(dp) :: dx, oneoverdx2, phi_b, nb_mask, nb_phi, w
  real(dp) :: fourpi
  
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

  call open_cache(r,g,m,operation_interpol,domain_decompos_amr)

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
        igridn=get_grid(r,g,m,hash_nbor,m%grid_dict,.false.,.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=m%grid(igridn)%phi(ind)
              dis_nbor(ind,inbor)=m%grid(igridn)%f(ind,3)
           end do

        ! Otherwise interpolate from coarser level
        else

           ! Get 3**ndim neighbouring parent cells using read-only cache
           call get_threetondim_nbor_parent_cell(r,g,m,hash_nbor,m%grid_dict,igrid_nbor,ind_nbor,.false.,.true.)
           call interpol_phi(m,igrid_nbor,ind_nbor,ccc,bbb,tfrac,phi_nbor(1,inbor))
           do ind=1,threetondim
              call unlock_cache(r,g,m,igrid_nbor(ind))
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

  call close_cache(r,g,m,m%grid_dict)

end subroutine make_bc_rhs

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

#endif

