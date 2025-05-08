module rt_flag_module
contains
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine rt_flag(s,ilevel)
  use amr_parameters, only: ndim,twotondim,twondim,dp
  use amr_commons, only: oct,nbor
  use ramses_commons, only: ramses_t
  use rt_parameters, only: nrtvar, nrtgrp
  use cache_commons
  use cache
  use nbors_utils
  use boundaries, only: init_bound_refine
  implicit none
  type(ramses_t)::s
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine flag for refinement cells that satisfies
  ! some user-defined physical criteria at the level ilevel. 
  ! -------------------------------------------------------------------
  integer,dimension(1:3,1:8),save::iii=reshape(&
       & (/0,0,0,1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1/),(/3,8/))
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer::igrid,ind,idim,ivar,i_nbor
  integer::icelld,icellg,icellp,igroup
  integer,dimension(1:twondim)::icelln
  integer(kind=8),dimension(0:ndim)::hash_key,hash_nbor
  real(dp),dimension(1:nrtvar)::uug,uum,uud
  logical::ok, do_rt_refine
  type(nbor),dimension(1:twondim)::gridn
  real(dp),dimension(1:twondim)::c_factor
  type(oct),pointer::gridp
  type(msg_realdp)::dummy_realdp

  associate(r=>s%r,g=>s%g,m=>s%m)

  do_rt_refine=.false.
  do igroup=1, nrtgrp
    if( r%rt_err_grad_cn(igroup) .ne. -1.0 ) do_rt_refine=.true.
  end do
  if(.not. do_rt_refine) return ! No refinement done on radiation vars

  c_factor(:)=g%rt_c(ilevel)

  hash_key(0)=ilevel+1

  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_realdp)/32,&
                     pack=pack_fetch_rt,unpack=unpack_fetch_rt,&
                     bound=init_bound_refine)

  ! Loop over active grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Loop over cells
     do ind=1,twotondim

        ! Compute cell hash key
        hash_key(1:ndim)=2*m%grid(igrid)%ckey(1:ndim)+iii(1:ndim,ind)

        ! Initialize refinement to false
        ok=.false.

        ! If a neighbor cell does not exist,
        ! replace it by its father cell
        do i_nbor=1,twondim
           hash_nbor(0)=hash_key(0)
           ! Periodic boundary conditions
           do idim=1,ndim
              hash_nbor(idim)=hash_key(idim)+shift(idim,i_nbor)
              if(r%periodic(idim))then
                 if(hash_nbor(idim)<m%box_ckey_min(idim,ilevel+1))hash_nbor(idim)=m%box_ckey_max(idim,ilevel+1)-1
                 if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel+1))hash_nbor(idim)=m%box_ckey_min(idim,ilevel+1)
              endif
           enddo
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icellp,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
           if(associated(gridp))then
              gridn(i_nbor)%p=>gridp
              icelln(i_nbor)=icellp
              c_factor(i_nbor) = g%rt_c(ilevel)
           else
              hash_nbor(0)=hash_nbor(0)-1
              hash_nbor(1:ndim)=hash_nbor(1:ndim)/2
              call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icellp,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
              gridn(i_nbor)%p=>gridp
              icelln(i_nbor)=icellp
              c_factor(i_nbor) = g%rt_c(ilevel-1)
           endif
        end do

        ! Loop over dimensions
        do idim=1,ndim
           ! Gather rt variables
           do ivar=1,nrtvar
              icellg=icelln(2*idim-1)
              icelld=icelln(2*idim  )
#ifdef RT
              uug(ivar)=gridn(2*idim-1)%p%rtuold(icellg,ivar)*c_factor(2*idim-1)
              uum(ivar)=m%grid(igrid)%rtuold(ind,ivar)*g%rt_c(ilevel)
              uud(ivar)=gridn(2*idim)%p%rtuold(icelld,ivar)*c_factor(2*idim)
#endif 
           end do
           call rt_refine(r,uug,uum,uud,ok)
        end do

        do i_nbor=1,twondim
           call unlock_cache(s,gridn(i_nbor)%p)
        end do

        ! Count only newly flagged cells
        if(m%grid(igrid)%flag1(ind)==0.and.ok)g%nflag=g%nflag+1
        if(ok)m%grid(igrid)%flag1(ind)=1

     end do
     ! End loop over cells
  end do
  ! End loop over grids

  call close_cache(s,m%grid_dict)

  end associate

end subroutine rt_flag
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine pack_fetch_rt(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use rt_parameters, only: nrtvar
  use amr_commons, only: oct
  use cache_commons, only: msg_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_realdp)::msg

  do ind=1,twotondim
     if(grid%refined(ind))then
        msg%int4(ind)=1
     else
        msg%int4(ind)=0
     endif
  end do

  do ivar=1,nrtvar
     do ind=1,twotondim
#ifdef RT  
        msg%realdp_rt(ind,ivar)=grid%rtuold(ind,ivar)
#endif
     end do
  end do

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_rt
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine unpack_fetch_rt(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use rt_parameters, only: nrtvar
  use amr_commons, only: oct
  use cache_commons, only: msg_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  type(msg_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     if(msg%int4(ind)==1)then
        grid%refined(ind)=.true.
     else
        grid%refined(ind)=.false.
     endif
  end do

  do ivar=1,nrtvar
     do ind=1,twotondim
#ifdef RT
        grid%rtuold(ind,ivar)=msg%realdp_rt(ind,ivar)
#endif
     end do
  end do

end subroutine unpack_fetch_rt
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
end module rt_flag_module
