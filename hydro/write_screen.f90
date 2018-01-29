subroutine write_screen(m,ilevel)
  use amr_parameters, only: dp,ndim
  use amr_commons, only: mesh_t
  implicit none
  type(mesh_t)::m
  integer::ilevel

  ! Local variables
  integer::igrid,ind

#ifdef HYDRO

  if(ndim>1)return

  if(m%noct_tot(ilevel)>0)then
     write(*,*)'================================'
     write(*,'(" Level",I4," has ",I4, " grids.")')ilevel,m%noct(ilevel)
     do igrid=m%head(ilevel),m%tail(ilevel)
        do ind=1,2
           write(*,'(I4,1X,I8,1X,I4,1X,I4,1X,I4,1X,L,1X,2(I4,1X),6(1PE10.3))')&
                & igrid,m%grid(igrid)%ckey(1),ind,m%grid(igrid)%flag1(ind),&
                & m%grid(igrid)%lev,m%grid(igrid)%refined(ind),&
                & m%grid(igrid)%hkey(1),m%grid(igrid)%superoct, &
                & (2*m%grid(igrid)%ckey(1)+ind-0.5)/(2.*m%ckey_max(ilevel)),&
                & m%grid(igrid)%uold(ind,1),m%grid(igrid)%uold(ind,2)/m%grid(igrid)%uold(ind,1),&
                & m%grid(igrid)%uold(ind,3)
        end do
     end do
     write(*,*)'================================'
  endif

#endif

end subroutine write_screen
