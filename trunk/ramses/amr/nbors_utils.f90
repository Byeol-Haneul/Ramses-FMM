!###############################################################
!###############################################################
!###############################################################
!###############################################################
integer function get_parent_cell(hash_key) result(parent_cell)
  use amr_commons
  use hash
  implicit none
  integer(kind=8),dimension(0:ndim)::hash_key
  !
  integer(kind=8),dimension(0:ndim)::hash_father
  integer(kind=8),dimension(1:ndim)::ii
  integer::ind,ipos,idim,get_grid

  hash_father(0)=hash_key(0)-1
  hash_father(1:ndim)=hash_key(1:ndim)/2
  ii(1:ndim)=hash_key(1:ndim)-2*hash_father(1:ndim)
  ind=1
  do idim=1,ndim
     ind=ind+2**(idim-1)*ii(idim)
  end do
  ipos=get_grid(hash_father)
  parent_cell=0
  if(ipos>0)parent_cell=(ipos-1)*twotondim+ind
end function get_parent_cell
!##############################################################
!##############################################################
!##############################################################
!##############################################################
integer function get_grid(hash_key) result(child_grid)
  use amr_commons
  use hash
  implicit none
  integer(kind=8),dimension(0:ndim)::hash_key
  !
  child_grid=hash_get(grid_dict,hash_key)
end function get_grid
!##############################################################
!##############################################################
!##############################################################
!##############################################################
