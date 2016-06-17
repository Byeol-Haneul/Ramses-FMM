subroutine write_screen(ilevel)
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !
  integer::igrid,ilevel,ind

  if(ndim>1)return

  if(noct_tot(ilevel)>0)then
     write(*,*)'================================'
     write(*,'(" Level",I4," has ",I4, " grids.")')ilevel,noct(ilevel)
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,2
           write(*,'(I4,1X,I8,1X,I4,1X,I4,1X,I4,1X,L,1X,2(I4,1X),6(1PE10.3))')&
                & igrid,grid(igrid)%ckey(1),ind,grid(igrid)%flag1(ind),grid(igrid)%lev,grid(igrid)%refined(ind),&
                & grid(igrid)%hkey(1),grid(igrid)%superoct, &
                & (2*grid(igrid)%ckey(1)+ind-0.5)/(2.*ckey_max(ilevel)),&
                & fluid(igrid)%uold(ind,1),fluid(igrid)%uold(ind,2)/fluid(igrid)%uold(ind,1),fluid(igrid)%uold(ind,3)
        end do
     end do
     write(*,*)'================================'
  endif

111 format(2(1pe12.5,1x))
112 format(i3,1x,1pe10.3,1x,8(1pe10.3,1x))
113 format(i3,1x,1pe12.5,1x,9(1pe10.3,1x))
114 format(' Output ',i5,' cells')
115 format(' Output ',i5,' parts')

end subroutine write_screen
