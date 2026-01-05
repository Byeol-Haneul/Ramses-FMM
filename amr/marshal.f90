module marshal
contains
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine pack_fetch_refine(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: ndim, twotondim
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtvar
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::idim,ind,ivar
  type(msg_large_realdp)::msg

  do ind=1,twotondim
     if(mesh%grid(igrid)%refined(ind))then
        msg%int4(ind)=1
     else
        msg%int4(ind)=0
     endif
  end do
  
#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        msg%realdp_hydro(ind,ivar)=mesh%uold(ind,ivar,igrid)
     end do
  end do
#endif
  
#ifdef MHD
  do ivar=1,6
     do ind=1,twotondim
        msg%realdp_mhd(ind,ivar)=mesh%bold(ind,ivar,igrid)
     end do
  end do
#endif

#ifdef GRAV
  do idim=1,ndim
     do ind=1,twotondim
        msg%realdp_poisson(ind,idim)=mesh%f(ind,idim,igrid)
     end do
  end do
  do ind=1,twotondim
     msg%realdp_poisson(ind,ndim+1)=mesh%phi(ind,igrid)
     msg%realdp_poisson(ind,ndim+2)=mesh%phi_old(ind,igrid)
  end do
#endif

#ifdef RT
  do ivar=1,nrtvar
     do ind=1,twotondim
        msg%realdp_rt(ind,ivar)=mesh%rtuold(ind,ivar,igrid)
     end do
  end do
#endif
  
  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_refine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine unpack_fetch_refine(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtvar
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::idim,ind,ivar
  type(msg_large_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     if(msg%int4(ind)==1)then
        mesh%grid(igrid)%refined(ind)=.true.
     else
        mesh%grid(igrid)%refined(ind)=.false.
     endif
  end do
  
#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        mesh%uold(ind,ivar,igrid)=msg%realdp_hydro(ind,ivar)
     end do
  end do
#endif

#ifdef MHD
  do ivar=1,6
     do ind=1,twotondim
        mesh%bold(ind,ivar,igrid)=msg%realdp_mhd(ind,ivar)
     end do
  end do
#endif

#ifdef GRAV
  do idim=1,ndim
     do ind=1,twotondim
        mesh%f(ind,idim,igrid)=msg%realdp_poisson(ind,idim)
     end do
  end do
  do ind=1,twotondim
     mesh%phi(ind,igrid)=msg%realdp_poisson(ind,ndim+1)
     mesh%phi_old(ind,igrid)=msg%realdp_poisson(ind,ndim+2)
  end do
#endif

#ifdef RT
  do ivar=1,nrtvar
     do ind=1,twotondim
        mesh%rtuold(ind,ivar,igrid)=msg%realdp_rt(ind,ivar)
     end do
  end do
#endif


end subroutine unpack_fetch_refine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_flag(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_int4
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_int4)::msg

  do ind=1,twotondim
     msg%int4(ind)=mesh%flag1(ind,igrid)
  end do

  msg_array=transfer(msg,msg_array)
  
end subroutine pack_fetch_flag
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine unpack_fetch_flag(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_int4
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_int4)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     mesh%flag1(ind,igrid)=msg%int4(ind)
  end do
  
end subroutine unpack_fetch_flag
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_flag2(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_int4
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_int4)::msg

  do ind=1,twotondim
     msg%int4(ind)=mesh%flag2(ind,igrid)
  end do

  msg_array=transfer(msg,msg_array)
  
end subroutine pack_fetch_flag2
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine unpack_fetch_flag2(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_int4
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_int4)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     mesh%flag2(ind,igrid)=msg%int4(ind)
  end do
  
end subroutine unpack_fetch_flag2
!###############################################################
!###############################################################
!###############################################################
!###############################################################
end module marshal
