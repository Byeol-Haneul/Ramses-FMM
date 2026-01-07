module marshal
contains
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine pack_fetch_refine(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtvar
  use amr_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::idim,ind,ivar
  type(msg_large_realdp)::msg

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
        msg%realdp_hydro(ind,ivar)=grid%uold(ind,ivar)
     end do
  end do
#endif
  
#ifdef HYDRO
  msg%realdp_mflux=grid%mflux
#endif

#ifdef MHD
  msg%realdp_mhd=grid%bold
#endif

#ifdef GRAV
  do idim=1,ndim
     do ind=1,twotondim
        msg%realdp_poisson(ind,idim)=grid%f(ind,idim)
     end do
  end do
  do ind=1,twotondim
     msg%realdp_poisson(ind,ndim+1)=grid%phi(ind)
     msg%realdp_poisson(ind,ndim+2)=grid%phi_old(ind)
  end do
#endif

#ifdef RT
  do ivar=1,nrtvar
     do ind=1,twotondim
        msg%realdp_rt(ind,ivar)=grid%rtuold(ind,ivar)
     end do
  end do
#endif
  
  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_refine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine unpack_fetch_refine(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtvar
  use amr_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::idim,ind,ivar
  type(msg_large_realdp)::msg

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
        grid%uold(ind,ivar)=msg%realdp_hydro(ind,ivar)
     end do
  end do
#endif

#ifdef HYDRO
  grid%mflux=msg%realdp_mflux
#endif

#ifdef MHD
  grid%bold=msg%realdp_mhd
#endif

#ifdef GRAV
  do idim=1,ndim
     do ind=1,twotondim
        grid%f(ind,idim)=msg%realdp_poisson(ind,idim)
     end do
  end do
  do ind=1,twotondim
     grid%phi(ind)=msg%realdp_poisson(ind,ndim+1)
     grid%phi_old(ind)=msg%realdp_poisson(ind,ndim+2)
  end do
#endif

#ifdef RT
  do ivar=1,nrtvar
     do ind=1,twotondim
        grid%rtuold(ind,ivar)=msg%realdp_rt(ind,ivar)
     end do
  end do
#endif


end subroutine unpack_fetch_refine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_flag(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_int4)::msg

  do ind=1,twotondim
     msg%int4(ind)=grid%flag1(ind)
  end do

  msg_array=transfer(msg,msg_array)
  
end subroutine pack_fetch_flag
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine unpack_fetch_flag(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_int4)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     grid%flag1(ind)=msg%int4(ind)
  end do
  
end subroutine unpack_fetch_flag
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_flag2(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_int4)::msg

  do ind=1,twotondim
     msg%int4(ind)=grid%flag2(ind)
  end do

  msg_array=transfer(msg,msg_array)
  
end subroutine pack_fetch_flag2
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine unpack_fetch_flag2(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_int4
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_int4)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     grid%flag2(ind)=msg%int4(ind)
  end do
  
end subroutine unpack_fetch_flag2
!###############################################################
!###############################################################
!###############################################################
!###############################################################
end module marshal
