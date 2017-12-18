!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine hydro_flag(r,g,m,ilevel)
  use amr_parameters, only: ndim,twotondim,twondim,dp
  use amr_commons, only: run_t,global_t,mesh_t
  use hydro_parameters, only: nvar
  use cache_commons
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
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
  integer::parent_cell,get_parent_cell
  integer::igridd,igridg,indd,indg,igridp
  integer,dimension(1:twondim),save::indn
  integer(kind=8),dimension(0:ndim)::hash_key,hash_nbor
  real(dp),dimension(1:nvar),save::uug,uum,uud
  logical::ok

#ifdef HYDRO

  hash_key(0)=ilevel+1

  if(    r%err_grad_d==-1.0.and.&
       & r%err_grad_p==-1.0.and.&
       & r%err_grad_u==-1.0)return

  call open_cache(r,g,m,operation_hydro,domain_decompos_amr)

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
           ! Periodic boundary conditons
           do idim=1,ndim
              hash_nbor(idim)=hash_key(idim)+shift(idim,i_nbor)
              if(hash_nbor(idim)<0)hash_nbor(idim)=m%ckey_max(ilevel+1)-1
              if(hash_nbor(idim)==m%ckey_max(ilevel+1))hash_nbor(idim)=0
           enddo
           parent_cell=get_parent_cell(r,g,m,hash_nbor,m%grid_dict,.false.,.true.)
           if(parent_cell>0)then
              indn(i_nbor)=parent_cell
           else
              hash_nbor(0)=hash_nbor(0)-1
              hash_nbor(1:ndim)=hash_nbor(1:ndim)/2
              parent_cell=get_parent_cell(r,g,m,hash_nbor,m%grid_dict,.false.,.true.)
              indn(i_nbor)=parent_cell
           endif
           igridp=(indn(i_nbor)-1)/twotondim+1
           call lock_cache(r,g,m,igridp)
        end do

        ! Loop over dimensions
        do idim=1,ndim
           ! Gather hydro variables
           do ivar=1,nvar
              igridg=(indn(2*idim-1)-1)/twotondim+1
              igridd=(indn(2*idim  )-1)/twotondim+1
              indg=indn(2*idim-1)-(igridg-1)*twotondim
              indd=indn(2*idim  )-(igridd-1)*twotondim
              uug(ivar)=m%grid(igridg)%uold(indg,ivar)
              uum(ivar)=m%grid(igrid )%uold(ind ,ivar)
              uud(ivar)=m%grid(igridd)%uold(indd,ivar)
           end do
           call hydro_refine(r,uug,uum,uud,ok)
        end do
        
        do i_nbor=1,twondim
           igridp=(indn(i_nbor)-1)/twotondim+1
           call unlock_cache(r,g,m,igridp)
        end do

        ! Count only newly flagged cells
        if(m%grid(igrid)%flag1(ind)==0.and.ok)g%nflag=g%nflag+1
        if(ok)m%grid(igrid)%flag1(ind)=1

     end do
     ! End loop over cells
  end do
  ! End loop over grids

  call close_cache(r,g,m,m%grid_dict)

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
  integer,dimension(1:msg_size)::msg_array

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

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_hydro
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine unpack_fetch_hydro(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size)::msg_array

  integer::ind,ivar
  type(msg_realdp)::msg

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
     do ind=1,twotondin
        grid%uold(ind,ivar)=msg%realdp(ind,ivar)
     end do
  end do
#endif

end subroutine unpack_fetch_hydro
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
