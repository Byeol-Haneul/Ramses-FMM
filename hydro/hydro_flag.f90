module hydro_flag_module
contains
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine hydro_flag(s,ilevel)
  use amr_parameters, only: ndim,twotondim,twondim,dp
  use amr_commons, only: oct,nbor
  use ramses_commons, only: ramses_t
  use hydro_parameters, only: nvar
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
  integer::igridd,igridg,icelld,icellg,igridp,icellp
  integer,dimension(1:twondim)::igridn,icelln
  integer(kind=8),dimension(0:ndim)::hash_key,hash_nbor
  real(dp),dimension(1:nvar)::uug,uum,uud
#ifdef MHD
  real(dp),dimension(1:6)::bbg,bbm,bbd
#endif
  logical::ok
  type(nbor),dimension(1:twondim)::gridn
  type(oct),pointer::gridp
  type(msg_realdp)::dummy_realdp

#ifdef HYDRO

  associate(r=>s%r,g=>s%g,m=>s%m)

  if(    r%err_grad_d==-1.0.and.&
#ifdef MHD
       & r%err_grad_A==-1.0.and.&
       & r%err_grad_B==-1.0.and.&
       & r%err_grad_C==-1.0.and.&
       & r%err_grad_B2==-1.0.and.&
#endif       
       & r%err_grad_p==-1.0.and.&
       & r%err_grad_u==-1.0.and.&
       & r%err_grad_xHI==-1.0.and.&
       & r%err_grad_xHII==-1.0)return

  hash_key(0)=ilevel+1

  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_realdp)/32,&
                     pack=pack_fetch_hydro,unpack=unpack_fetch_hydro,&
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
           else
              hash_nbor(0)=hash_nbor(0)-1
              hash_nbor(1:ndim)=hash_nbor(1:ndim)/2
              call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icellp,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
              gridn(i_nbor)%p=>gridp
              icelln(i_nbor)=icellp
           endif
        end do

        ! Loop over dimensions
        do idim=1,ndim
           ! Gather hydro variables
           do ivar=1,nvar
              icellg=icelln(2*idim-1)
              icelld=icelln(2*idim  )
              uug(ivar)=gridn(2*idim-1)%p%uold(icellg,ivar)
              uum(ivar)=m%grid(igrid)%uold(ind,ivar)
              uud(ivar)=gridn(2*idim)%p%uold(icelld,ivar)
           end do
#ifdef MHD
           ! Gather MHD variables
           do ivar=1,6
              icellg=icelln(2*idim-1)
              icelld=icelln(2*idim  )
              bbg(ivar)=gridn(2*idim-1)%p%bold(icellg,ivar)
              bbm(ivar)=m%grid(igrid)%bold(ind,ivar)
              bbd(ivar)=gridn(2*idim)%p%bold(icelld,ivar)
           end do
           call hydro_refine(r,uug,uum,uud,bbg,bbm,bbd,ok)
#else
           call hydro_refine(r,uug,uum,uud,ok)
#endif
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

#endif

end subroutine hydro_flag
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine pack_fetch_hydro(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
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
  
#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        msg%realdp(ind,ivar)=grid%uold(ind,ivar)
     end do
  end do
#endif

#ifdef MHD
  msg%realdp_mhd=grid%bold
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_hydro
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine unpack_fetch_hydro(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
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
  
#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        grid%uold(ind,ivar)=msg%realdp(ind,ivar)
     end do
  end do
#endif

#ifdef MHD
  grid%bold=msg%realdp_mhd
#endif

end subroutine unpack_fetch_hydro
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
end module hydro_flag_module
