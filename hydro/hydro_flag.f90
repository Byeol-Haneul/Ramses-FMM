#ifdef TOTO
subroutine hydro_flag(ilevel)
  use amr_commons
  use hydro_commons
  use hash
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine flag for refinement cells that satisfies
  ! some user-defined physical criteria at the level ilevel. 
  ! -------------------------------------------------------------------
  integer,dimension(1:3,1:8),save::iii=reshape(&
       & (/0,0,0,1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1/),(/3,8/))
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer::igrid,ind,idim,ngrid,ivar,i_nbor
  integer::parent_cell,get_parent_cell
  integer::igridd,igridg,indd,indg,igridp
  integer,dimension(1:twondim),save::indn
  integer(kind=8),dimension(0:ndim)::hash_key,hash_nbor
  real(dp),dimension(1:nvar),save::uug,uum,uud
  logical::ok

#ifdef HYDRO

  if(ilevel==nlevelmax)return
  if(noct_tot(ilevel)==0)return

  hash_key(0)=ilevel+1

  if(    err_grad_d==-1.0.and.&
       & err_grad_p==-1.0.and.&
       & err_grad_u==-1.0)return

  call open_cache(operation_hydro,domain_decompos_amr)

  ! Loop over active grids
  do igrid=head(ilevel),tail(ilevel)

     ! Loop over cells
     do ind=1,twotondim

        ! Compute cell hash key
        hash_key(1:ndim)=2*grid(igrid)%ckey(1:ndim)+iii(1:ndim,ind)

        ! Initialize refinement to false
        ok=.false.

        ! If a neighbor cell does not exist,
        ! replace it by its father cell
        do i_nbor=1,twondim
           hash_nbor(0)=hash_key(0)
           ! Periodic boundary conditons
           do idim=1,ndim
              hash_nbor(idim)=hash_key(idim)+shift(idim,i_nbor)
              if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel+1)-1
              if(hash_nbor(idim)==ckey_max(ilevel+1))hash_nbor(idim)=0
           enddo
           parent_cell=get_parent_cell(hash_nbor,grid_dict,.false.,.true.)
           if(parent_cell>0)then
              indn(i_nbor)=parent_cell
           else
              hash_nbor(0)=hash_nbor(0)-1
              hash_nbor(1:ndim)=hash_nbor(1:ndim)/2
              parent_cell=get_parent_cell(hash_nbor,grid_dict,.false.,.true.)
              indn(i_nbor)=parent_cell
           endif
           igridp=(indn(i_nbor)-1)/twotondim+1
           call lock_cache(igridp)
        end do

        ! Loop over dimensions
        do idim=1,ndim
           ! Gather hydro variables
           do ivar=1,nvar
              igridg=(indn(2*idim-1)-1)/twotondim+1
              igridd=(indn(2*idim  )-1)/twotondim+1
              indg=indn(2*idim-1)-(igridg-1)*twotondim
              indd=indn(2*idim  )-(igridd-1)*twotondim
              uug(ivar)=grid(igridg)%uold(indg,ivar)
              uum(ivar)=grid(igrid )%uold(ind ,ivar)
              uud(ivar)=grid(igridd)%uold(indd,ivar)
           end do
           call hydro_refine(uug,uum,uud,ok)
        end do
        
        do i_nbor=1,twondim
           igridp=(indn(i_nbor)-1)/twotondim+1
           call unlock_cache(igridp)
        end do

        ! Count only newly flagged cells
        if(grid(igrid)%flag1(ind)==0.and.ok)nflag=nflag+1
        if(ok)grid(igrid)%flag1(ind)=1

     end do
     ! End loop over cells
  end do
  ! End loop over grids

  call close_cache(grid_dict)

#endif

end subroutine hydro_flag
#endif
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine hydro_flag_2(r,g,m,ilevel)
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
  integer::igrid,ind,idim,ngrid,ivar,i_nbor
  integer::parent_cell,get_parent_cell_2
  integer::igridd,igridg,indd,indg,igridp
  integer,dimension(1:twondim),save::indn
  integer(kind=8),dimension(0:ndim)::hash_key,hash_nbor
  real(dp),dimension(1:nvar),save::uug,uum,uud
  logical::ok

#ifdef HYDRO

  if(ilevel==r%nlevelmax)return
  if(m%noct_tot(ilevel)==0)return

  hash_key(0)=ilevel+1

  if(    r%err_grad_d==-1.0.and.&
       & r%err_grad_p==-1.0.and.&
       & r%err_grad_u==-1.0)return

  call open_cache_2(r,g,m,operation_hydro,domain_decompos_amr)

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
           parent_cell=get_parent_cell_2(r,g,m,hash_nbor,m%grid_dict,.false.,.true.)
           if(parent_cell>0)then
              indn(i_nbor)=parent_cell
           else
              hash_nbor(0)=hash_nbor(0)-1
              hash_nbor(1:ndim)=hash_nbor(1:ndim)/2
              parent_cell=get_parent_cell_2(r,g,m,hash_nbor,m%grid_dict,.false.,.true.)
              indn(i_nbor)=parent_cell
           endif
           igridp=(indn(i_nbor)-1)/twotondim+1
           call lock_cache_2(r,g,m,igridp)
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
           call hydro_refine_2(r,uug,uum,uud,ok)
        end do
        
        do i_nbor=1,twondim
           igridp=(indn(i_nbor)-1)/twotondim+1
           call unlock_cache_2(r,g,m,igridp)
        end do

        ! Count only newly flagged cells
        if(m%grid(igrid)%flag1(ind)==0.and.ok)g%nflag=g%nflag+1
        if(ok)m%grid(igrid)%flag1(ind)=1

     end do
     ! End loop over cells
  end do
  ! End loop over grids

  call close_cache_2(r,g,m,m%grid_dict)

#endif

end subroutine hydro_flag_2
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
