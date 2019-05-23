module smooth_module
contains
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_smooth_fine(pst,ilevel,input_size,noct,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:1)::ilevel,noct

  integer,dimension(1:1)::next_noct
  integer::nflag
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SMOOTH_FINE,pst%iUpper+1,input_size,output_size,ilevel)
     call r_smooth_fine(pst%pLower,ilevel,input_size,noct,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_noct)
     noct=noct+next_noct
  else
     call smooth_fine(pst%s,ilevel(1),nflag)
     noct(1)=nflag
  endif

end subroutine r_smooth_fine
!############################################################
!############################################################
!############################################################
!############################################################
subroutine smooth_fine(s,ilevel,nflag)
  use amr_parameters, only: ndim,twotondim,twondim
  use ramses_commons, only: ramses_t
  use cache_commons
  use amr_commons, only: nbor
  use nbors_utils_p
  implicit none
  type(ramses_t)::s
  integer::ilevel,nflag
  ! -------------------------------------------------------------------
  ! Dilatation operator.
  ! This routine makes one cell width cubic buffer around flag1 cells 
  ! at level ilevel by following these 3 steps:
  ! step 1: flag1 cells with at least 1 flag1 neighbors (if ndim > 0) 
  ! step 2: flag1 cells with at least 2 flag1 neighbors (if ndim > 1) 
  ! step 3: flag1 cells with at least 2 flag1 neighbors (if ndim > 2) 
  ! Array flag2 is used as temporary workspace.
  ! -------------------------------------------------------------------
  integer::ismooth,count_nbor,ig,in
  integer::igrid,idim,ind,i_nbor,igrid_nbor,icell_nbor
  integer,dimension(1:3),save::n_nbor=(/1,2,2/)
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer,dimension(0:twondim)::igridn
  type(nbor),dimension(0:twondim)::gridn

  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,&
       &   0,-1,0,0,1,0,&
       &   0,0,-1,0,0,1/),(/3,6/))
  integer,dimension(1:8,1:6),save::ggg=reshape(&
       & (/1,0,1,0,1,0,1,0,&
       &   0,2,0,2,0,2,0,2,&
       &   3,3,0,0,3,3,0,0,&
       &   0,0,4,4,0,0,4,4,&
       &   5,5,5,5,0,0,0,0,&
       &   0,0,0,0,6,6,6,6/),(/8,6/))
  integer,dimension(1:8,1:6),save::hhh=reshape(&
       & (/2,1,4,3,6,5,8,7,&
       &   2,1,4,3,6,5,8,7,&
       &   3,4,1,2,7,8,5,6,&
       &   3,4,1,2,7,8,5,6,&
       &   5,6,7,8,1,2,3,4,&
       &   5,6,7,8,1,2,3,4/),(/8,6/))

  associate(r=>s%r,g=>s%g,m=>s%m)

  hash_nbor(0)=ilevel

  ! Loop over steps
  do ismooth=1,ndim

     ! Initialize flag2 to 0
     do igrid=m%head(ilevel),m%tail(ilevel)
        do ind=1,twotondim
           m%grid(igrid)%flag2(ind)=0
        end do
     end do

     call open_cache(s,operation_smooth,domain_decompos_amr)

     ! Count neighbors and set flag2 accordingly
     do igrid=m%head(ilevel),m%tail(ilevel)

        ! Get neighboring octs
        igridn(0)=igrid
        gridn(0)%p => m%grid(igrid)
        do i_nbor=1,twondim
           hash_nbor(1:ndim)=m%grid(igrid)%ckey(1:ndim)+shift(1:ndim,i_nbor)
           ! Periodic boundary conditons
           do idim=1,ndim
              if(hash_nbor(idim)<0)hash_nbor(idim)=m%ckey_max(ilevel)-1
              if(hash_nbor(idim)==m%ckey_max(ilevel))hash_nbor(idim)=0
           enddo
           call get_grid_p(s,hash_nbor,m%grid_dict,gridn(i_nbor)%p,.false.,.true.)
           call lock_cache_p(s,gridn(i_nbor)%p)
        end do

        ! Count neighbors and set flag2 accordingly        
        do ind=1,twotondim
           count_nbor=0
           do in=1,twondim
              ig=ggg(ind,in)
              icell_nbor=hhh(ind,in)
              if(associated(gridn(ig)%p))then
                 count_nbor=count_nbor+gridn(ig)%p%flag1(icell_nbor)
              endif
           end do
           ! flag2 cell if necessary
           if(count_nbor>=n_nbor(ismooth))then
              m%grid(igrid)%flag2(ind)=1
           endif
        end do

        do i_nbor=1,twondim
           call unlock_cache_p(s,gridn(i_nbor)%p)
        end do

     end do
     ! End loop over grids

    call close_cache(s,m%grid_dict)

     ! Set flag1=1 for cells with flag2=1
     do igrid=m%head(ilevel),m%tail(ilevel)
        do ind=1,twotondim
           if(m%grid(igrid)%flag1(ind)==1)m%grid(igrid)%flag2(ind)=0
        end do
        do ind=1,twotondim
           if(m%grid(igrid)%flag2(ind)==1)then
              m%grid(igrid)%flag1(ind)=1
              g%nflag=g%nflag+1
           endif
        end do
     end do     

  end do
  ! End loop over steps

  nflag=g%nflag

  end associate
  
end subroutine smooth_fine
!############################################################
!############################################################
!############################################################
!############################################################
end module smooth_module
