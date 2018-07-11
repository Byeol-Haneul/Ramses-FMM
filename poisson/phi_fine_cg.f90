#ifdef GRAV
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine m_phi_fine_cg(pst,ilevel,icount)
  use amr_parameters, only: ndim,twondim,twotondim,threetondim,nvector,dp
  use ramses_commons, only: pst_t
  use cache_commons
  implicit none
  type(pst_t)::pst
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
  integer::igrid,ind,iter,itermax
  integer,dimension(1)::dummy
  real(kind=8)::error,error_ini
  real(kind=8)::r2_old=0,alpha_cg,beta_cg
  real(kind=8)::r2,pAp,rhs_norm
  integer,dimension(1:4)::input_array
  integer,dimension(1:4)::output_array

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,mdl=>pst%s%mdl)
    
  if(r%gravity_type>0)return
  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,111)ilevel
111 format('   Entering phi_fine_cg for level ',I2)

  !===============================
  ! Compute initial phi
  !===============================
  input_array(1)=ilevel
  input_array(2)=icount
  call r_make_initial_phi(pst,input_array,2,dummy,0)

  !===============================
  ! Compute right-hand side norm
  !===============================
  call r_cmp_rhs_norm(pst,input_array,1,output_array,2)
  rhs_norm=transfer(output_array(1:2),rhs_norm)
  rhs_norm=DSQRT(rhs_norm/dble(twotondim*m%noct_tot(ilevel)))

  !==============================================
  ! Compute r = b - Ax and store it into f(i,1)
  ! Also set p = r and store it into f(i,2)
  !==============================================
  input_array(1)=ilevel
  input_array(2)=icount
  call r_cmp_residual_cg(pst,input_array,2,dummy,0)

  !====================================
  ! Main iteration loop
  !====================================
  iter=0; itermax=10000
  error=1.0D0; error_ini=1.0D0
  do while(error>r%epsilon*error_ini.and.iter<itermax)

     iter=iter+1

     !====================================
     ! Compute residual norm
     !====================================
     call r_cmp_r2_cg(pst,input_array,1,output_array,2)
     r2=transfer(output_array(1:2),r2)

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
     input_array(1)=ilevel
     input_array(2:3)=transfer(beta_cg,input_array(2:3))
     call r_recurrence_on_p(pst,input_array,3,dummy,0)

     !==============================================
     ! Compute z = Ap and store it into f(i,3)
     !==============================================
     call r_cmp_Ap_cg(pst,input_array,1,dummy,0)

     !====================================
     ! Compute p.Ap scalar product
     !====================================
     call r_cmp_pAp_cg(pst,input_array,1,output_array,2)
     pAp=transfer(output_array(1:2),pAp)

     !====================================
     ! Compute alpha factor
     !====================================
     alpha_cg = r2/pAp

     !====================================
     ! Recurrence on x and r
     !====================================
     input_array(1)=ilevel
     input_array(2:3)=transfer(alpha_cg,input_array(2:3))
     call r_recurrence_x_and_r(pst,input_array,3,dummy,0)

     !====================================
     ! Compute error
     !====================================
     error=DSQRT(r2/dble(twotondim*m%noct_tot(ilevel)))
     if(iter==1)error_ini=error
     if(r%verbose)write(*,112)iter,error/rhs_norm,error/error_ini
112  format('   ==> Step=',i5,' Error=',2(1pe10.3,1x))

  end do
  ! End main iteration loop

  write(*,115)ilevel,iter,error/rhs_norm,error/error_ini
115 format('   ==> Level=',i5,' Step=',i5,' Error=',2(1pe10.3,1x))
  if(iter>=itermax)then
     write(*,*)'Poisson failed to converge...'
  end if

  end associate
  
end subroutine m_phi_fine_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_recurrence_on_p(pst,input_array,input_size,output_array,output_size)
  use amr_parameters, only: twotondim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ind,igrid,ilevel
  real(kind=8)::beta_cg

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_RECURRENCE_ON_P,pst%iUpper+1,input_size,output_size,input_array)
     call r_recurrence_on_p(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     ilevel=input_array(1)
     beta_cg=transfer(input_array(2:3),beta_cg)
     do igrid=pst%s%m%head(ilevel),pst%s%m%tail(ilevel)
        do ind=1,twotondim
           pst%s%m%grid(igrid)%f(ind,2)=pst%s%m%grid(igrid)%f(ind,1)+beta_cg*pst%s%m%grid(igrid)%f(ind,2)
        end do
     end do
  endif

end subroutine r_recurrence_on_p
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_recurrence_x_and_r(pst,input_array,input_size,output_array,output_size)
  use amr_parameters, only: twotondim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ind,igrid,ilevel
  real(kind=8)::alpha_cg

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_RECURRENCE_X_AND_R,pst%iUpper+1,input_size,output_size,input_array)
     call r_recurrence_x_and_r(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     ilevel=input_array(1)
     alpha_cg=transfer(input_array(2:3),alpha_cg)
     ! Recurrence on x
     do igrid=pst%s%m%head(ilevel),pst%s%m%tail(ilevel)
        do ind=1,twotondim
           pst%s%m%grid(igrid)%phi(ind)=pst%s%m%grid(igrid)%phi(ind)+alpha_cg*pst%s%m%grid(igrid)%f(ind,2)
        end do
     end do
     ! Recurrence on r
     do igrid=pst%s%m%head(ilevel),pst%s%m%tail(ilevel)
        do ind=1,twotondim
           pst%s%m%grid(igrid)%f(ind,1)=pst%s%m%grid(igrid)%f(ind,1)-alpha_cg*pst%s%m%grid(igrid)%f(ind,3)
        end do
     end do
  endif

end subroutine r_recurrence_x_and_r
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_cmp_residual_cg(pst,input_array,input_size,output_array,output_size)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel,icount

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_CMP_RESIDUAL_CG,pst%iUpper+1,input_size,output_size,input_array)
     call r_cmp_residual_cg(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     ilevel=input_array(1)
     icount=input_array(2)
     call cmp_residual_cg(pst%s,ilevel,icount)
  endif

end subroutine r_cmp_residual_cg

subroutine cmp_residual_cg(s,ilevel,icount)
  use amr_parameters, only: ndim,twondim,twotondim,threetondim,nvector,dp
  use ramses_commons, only: ramses_t
  use cache_commons
  implicit none
  type(ramses_t)::s
  integer::ilevel,icount
  !------------------------------------------------------------------
  ! This routine computes the residual for the Conjugate Gradient
  ! Poisson solver. The residual is stored in f(i,1) and f(i,2).
  !------------------------------------------------------------------
  integer,external::get_grid
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

  associate(r=>s%r,g=>s%g,m=>s%m)
    
  ! Set constants
  dx2=(r%boxlen/2**ilevel)**2
  fourpi=4.D0*ACOS(-1.0D0)
  if(r%cosmo)fourpi=1.5D0*g%omega_m*g%aexp
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
     call mdl_abort
  endif

  ! Compute fraction of time steps for interpolation
  if (g%dtold(ilevel-1)>0.0)then
     tfrac=g%dtnew(ilevel)/g%dtold(ilevel-1)*(icount-1)
  else
     tfrac=0.0
  end if

  call open_cache(s,operation_interpol,domain_decompos_amr)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Get central oct potential
     do ind=1,twotondim
        phi_nbor(ind,0)=m%grid(igrid)%phi(ind)
     end do

     ! Get neighboring octs potential
     do i_nbor=1,twondim

        ! Get neighboring grid
        hash_nbor(1:ndim)=m%grid(igrid)%ckey(1:ndim)+shift(1:ndim,i_nbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=m%ckey_max(ilevel)-1
           if(hash_nbor(idim)==m%ckey_max(ilevel))hash_nbor(idim)=0
        enddo

        ! Get neighbouring grid using a read-only cache
        igridn=get_grid(s,hash_nbor,m%grid_dict,.false.,.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,i_nbor)=m%grid(igridn)%phi(ind)
           end do

        ! Otherwise interpolate from coarser level
        else
           ! Get 3**ndim neighbouring parent cell using a read-only cache
           call get_threetondim_nbor_parent_cell(s,hash_nbor,m%grid_dict,igrid_nbor,ind_nbor,.false.,.true.)
           call interpol_phi(m,igrid_nbor,ind_nbor,ccc,bbb,tfrac,phi_nbor(1,i_nbor))
           do ind=1,threetondim
              call unlock_cache(s,igrid_nbor(ind))
           end do
        endif

     end do
     ! End loop over neighboring octs

     ! Loop over cells
     do ind=1,twotondim

        ! Compute residual using 6 neighbors potential
        residu=m%grid(igrid)%phi(ind)
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           residu=residu-oneoversix*(phi_nbor(id1,ig1)+phi_nbor(id2,ig2))
        end do
        residu=residu+fact*(m%grid(igrid)%rho(ind)-g%rho_tot)

        ! Store results in f(ind,1)
        m%grid(igrid)%f(ind,1)=residu

        ! Store results in f(ind,2)
        m%grid(igrid)%f(ind,2)=residu

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache(s,m%grid_dict)

  end associate

end subroutine cmp_residual_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_cmp_Ap_cg(pst,input_array,input_size,output_array,output_size)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_CMP_AP_CG,pst%iUpper+1,input_size,output_size,input_array)
     call r_cmp_Ap_cg(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     ilevel=input_array(1)
     call cmp_Ap_cg(pst%s,ilevel)
  endif

end subroutine r_cmp_Ap_cg

subroutine cmp_Ap_cg(s,ilevel)
  use amr_parameters, only: ndim,twondim,twotondim,threetondim,nvector,dp
  use ramses_commons, only: ramses_t
  use cache_commons
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes Ap for the Conjugate Gradient
  ! Poisson Solver and store the result into f(i,3).
  !------------------------------------------------------------------
  integer,external::get_grid
  integer::inbor,igrid,idim,ind,igridn
  integer::id1,id2,ig1,ig2
  real(dp)::oneoversix,residu
  integer,dimension(1:3,1:2,1:8)::iii,jjj
  real(dp),dimension(1:twotondim,0:twondim),save::phi_nbor
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))

  associate(r=>s%r,g=>s%g,m=>s%m)

  ! Set constants
  oneoversix=1.0D0/dble(twondim)

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  call open_cache(s,operation_cg,domain_decompos_amr)

  hash_nbor(0)=ilevel

  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Get central oct potential
     do ind=1,twotondim
        phi_nbor(ind,0)=m%grid(igrid)%f(ind,2)
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
        igridn=get_grid(s,hash_nbor,m%grid_dict,.false.,.true.)

        ! If grid exists, then copy into array
        if(igridn>0)then
           do ind=1,twotondim
              phi_nbor(ind,inbor)=m%grid(igridn)%f(ind,2)
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
        residu=-m%grid(igrid)%f(ind,2)
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           residu=residu+oneoversix*(phi_nbor(id1,ig1)+phi_nbor(id2,ig2))
        end do

        ! Store results in f(ind,3)
        m%grid(igrid)%f(ind,3)=residu

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache(s,m%grid_dict)

  end associate

end subroutine cmp_Ap_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_make_initial_phi(pst,input_array,input_size,output_array,output_size)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::ilevel,icount
  
  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_MAKE_INITIAL_PHI,pst%iUpper+1,input_size,output_size,input_array)
     call r_make_initial_phi(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     ilevel=input_array(1)
     icount=input_array(2)
     call make_initial_phi(pst%s,ilevel,icount)
  endif

end subroutine r_make_initial_phi
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine make_initial_phi(s,ilevel,icount)
  use amr_parameters, only: ndim,twondim,twotondim,threetondim,nvector,dp
  use ramses_commons, only: ramses_t
  use cache_commons
  implicit none
  type(ramses_t)::s
  integer::ilevel,icount
  !
  !
  !
  integer::igrid,idim,ind
  integer,dimension(1:8,1:8)::ccc
  real(dp)::aa,bb,cc,dd,tfrac
  real(dp),dimension(1:8)::bbb
  integer(kind=8),dimension(0:ndim)::hash_key
  integer,dimension(1:threetondim),save::igrid_nbor,ind_nbor
  real(dp),dimension(1:twotondim),save::phi_int

  associate(r=>s%r,g=>s%g,m=>s%m)

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
     call mdl_abort
  endif

  ! Compute fraction of time steps for interpolation
  if (g%dtold(ilevel-1)>0.0)then
     tfrac=g%dtnew(ilevel)/g%dtold(ilevel-1)*(icount-1)
  else
     tfrac=0.0
  end if

  call open_cache(s,operation_interpol,domain_decompos_amr)

  hash_key(0)=ilevel

  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! By default, initial phi is equal to zero

     ! Loop over cells
     do ind=1,twotondim
        m%grid(igrid)%phi(ind)=0.0d0
        do idim=1,ndim
           m%grid(igrid)%f(ind,idim)=0.0
        end do
     end do
     ! End loop over cells

     ! For fine levels, initial phi is interpolated from coarser level
     if(ilevel.GT.r%levelmin)then
        
        hash_key(1:ndim)=m%grid(igrid)%ckey(1:ndim)
        ! Get 3**ndim neghbouring parent cell using read-only cache
        call get_threetondim_nbor_parent_cell(s,hash_key,m%grid_dict,igrid_nbor,ind_nbor,.false.,.true.)
        call interpol_phi(m,igrid_nbor,ind_nbor,ccc,bbb,tfrac,phi_int)
        do ind=1,threetondim
           call unlock_cache(s,igrid_nbor(ind))
        end do

        ! Loop over cells
        do ind=1,twotondim
           m%grid(igrid)%phi(ind)=phi_int(ind)
        end do
        ! End loop over cells

     end if

  end do
  ! End loop over grids

  call close_cache(s,m%grid_dict)

  end associate

end subroutine make_initial_phi

subroutine pack_fetch_interpol(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_three_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_three_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp_phi(ind)=grid%phi(ind)
     msg%realdp_phi_old(ind)=grid%phi_old(ind)
     msg%realdp_dis(ind)=grid%f(ind,3)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_interpol

subroutine unpack_fetch_interpol(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_three_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_three_realdp)::msg

  msg=transfer(msg_array,msg)
  
#ifdef GRAV
  do ind=1,twotondim
     grid%phi(ind)=msg%realdp_phi(ind)
     grid%phi_old(ind)=msg%realdp_phi_old(ind)
     grid%f(ind,3)=msg%realdp_dis(ind)
  end do
#endif

end subroutine unpack_fetch_interpol
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################
subroutine pack_fetch_cg(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp(ind)=grid%f(ind,2)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_cg
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################
subroutine unpack_fetch_cg(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

  msg=transfer(msg_array,msg)
  
#ifdef GRAV
  do ind=1,twotondim
     grid%f(ind,2)=msg%realdp(ind)
  end do
#endif

end subroutine unpack_fetch_cg
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################
recursive subroutine r_cmp_rhs_norm(pst,input_array,input_size,output_array,output_size)
  use amr_parameters, only: twotondim,twondim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array
  ! ------------------------------------------------------------------------
  ! Compute norm of residual 
  ! ------------------------------------------------------------------------
  integer::ind,igrid,ilevel
  integer,dimension(1:output_size)::next_output_array
  real(kind=8)::rhs_norm,next_rhs_norm
  real(kind=8)::dx2,fourpi,oneoversix,fact,fact2

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_CMP_RHS_NORM,pst%iUpper+1,input_size,output_size,input_array)
     call r_cmp_rhs_norm(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size,next_output_array)
     rhs_norm=transfer(output_array,rhs_norm)
     next_rhs_norm=transfer(next_output_array,next_rhs_norm)
     rhs_norm=rhs_norm+next_rhs_norm
     output_array=transfer(rhs_norm,output_array)
  else
     ilevel=input_array(1)
     ! Set constants
     dx2=(pst%s%r%boxlen/2.0d0**ilevel)**2
     fourpi=4.D0*ACOS(-1.0D0)
     if(pst%s%r%cosmo)fourpi=1.5D0*pst%s%g%omega_m*pst%s%g%aexp
     oneoversix=1.0D0/dble(twondim)
     fact=oneoversix*fourpi*dx2
     fact2=fact*fact
     rhs_norm=0.d0
     do igrid=pst%s%m%head(ilevel),pst%s%m%tail(ilevel)
        do ind=1,twotondim
           rhs_norm=rhs_norm+fact2*(pst%s%m%grid(igrid)%rho(ind)-pst%s%g%rho_tot)**2
        end do
     end do
     output_array=transfer(rhs_norm,output_array)
  endif

end subroutine r_cmp_rhs_norm
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################
recursive subroutine r_cmp_r2_cg(pst,input_array,input_size,output_array,output_size)
  use amr_parameters, only: twotondim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array
  ! ------------------------------------------------------------------------
  ! Compute norm of residual 
  ! ------------------------------------------------------------------------
  integer::ind,igrid,ilevel
  integer,dimension(1:output_size)::next_output_array
  real(kind=8)::r2,next_r2

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_CMP_R2_CG,pst%iUpper+1,input_size,output_size,input_array)
     call r_cmp_r2_cg(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size,next_output_array)
     r2=transfer(output_array,r2)
     next_r2=transfer(next_output_array,next_r2)
     r2=r2+next_r2
     output_array=transfer(r2,output_array)
  else
     ilevel=input_array(1)
     r2=0.0d0
     do igrid=pst%s%m%head(ilevel),pst%s%m%tail(ilevel)
        do ind=1,twotondim
           r2=r2+pst%s%m%grid(igrid)%f(ind,1)**2
        end do
     end do
     output_array=transfer(r2,output_array)
  endif

end subroutine r_cmp_r2_cg
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################
recursive subroutine r_cmp_pAp_cg(pst,input_array,input_size,output_array,output_size)
  use amr_parameters, only: twotondim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array
  ! ------------------------------------------------------------------------
  ! Compute norm of residual 
  ! ------------------------------------------------------------------------
  integer::igrid,ind,ilevel
  integer,dimension(1:output_size)::next_output_array
  real(kind=8)::pAp,next_pAp

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_CMP_PAP_CG,pst%iUpper+1,input_size,output_size,input_array)
     call r_cmp_pAp_cg(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size,next_output_array)
     pAp=transfer(output_array,pAp)
     next_pAp=transfer(next_output_array,next_pAp)
     pAp=pAp+next_pAp
     output_array=transfer(pAp,output_array)
  else
     ilevel=input_array(1)
     pAp=0.0d0
     do igrid=pst%s%m%head(ilevel),pst%s%m%tail(ilevel)
        do ind=1,twotondim
           pAp=pAp+pst%s%m%grid(igrid)%f(ind,2)*pst%s%m%grid(igrid)%f(ind,3)
        end do
     end do
     output_array=transfer(pAp,output_array)
  endif

end subroutine r_cmp_pAp_cg
! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

#endif

