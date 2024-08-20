module rt_flag_module
contains
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
        msg%realdp_rt(ind,ivar)=grid%rtuold(ind,ivar)
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
        grid%rtuold(ind,ivar)=msg%realdp_rt(ind,ivar)
     end do
  end do

end subroutine unpack_fetch_rt
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
end module rt_flag_module