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

#ifdef HYDRO

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
                & grid(igrid)%uold(ind,1),grid(igrid)%uold(ind,2)/grid(igrid)%uold(ind,1),grid(igrid)%uold(ind,3)
        end do
     end do
     write(*,*)'================================'
  endif

#endif

111 format(2(1pe12.5,1x))
112 format(i3,1x,1pe10.3,1x,8(1pe10.3,1x))
113 format(i3,1x,1pe12.5,1x,9(1pe10.3,1x))
114 format(' Output ',i5,' cells')
115 format(' Output ',i5,' parts')

end subroutine write_screen

subroutine write_screen_2(r,g,m,ilevel)
  use amr_parameters, only: dp,ndim
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
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
                & igrid,m%grid(igrid)%ckey(1),ind,m%grid(igrid)%flag1(ind),m%grid(igrid)%lev,m%grid(igrid)%refined(ind),&
                & m%grid(igrid)%hkey(1),m%grid(igrid)%superoct, &
                & (2*m%grid(igrid)%ckey(1)+ind-0.5)/(2.*m%ckey_max(ilevel)),&
                & m%grid(igrid)%uold(ind,1),m%grid(igrid)%uold(ind,2)/m%grid(igrid)%uold(ind,1),m%grid(igrid)%uold(ind,3)
        end do
     end do
     write(*,*)'================================'
  endif

#endif

111 format(2(1pe12.5,1x))
112 format(i3,1x,1pe10.3,1x,8(1pe10.3,1x))
113 format(i3,1x,1pe12.5,1x,9(1pe10.3,1x))
114 format(' Output ',i5,' cells')
115 format(' Output ',i5,' parts')

end subroutine write_screen_2
