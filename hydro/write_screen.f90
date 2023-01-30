subroutine write_screen(r,m)
  use amr_parameters, only: dp,ndim
  use amr_commons, only: mesh_t, run_t
  implicit none
  type(run_t)::r
  type(mesh_t)::m

  ! Local variables
  integer::ilevel,nleaf,i
  integer,dimension(:),allocatable::ll,ii
  real(kind=8),dimension(:),allocatable::xx,dd,uu,pp

  integer::igrid,ind
  logical::leaf
  
#ifdef HYDRO

  if(ndim>1)return

  ! Compute number of cells
  nleaf = 0
  do ilevel=r%levelmin,r%nlevelmax
     if(m%noct_tot(ilevel)>0)then
        do igrid=m%head(ilevel),m%tail(ilevel)
           do ind=1,2
              leaf = .not. m%grid(igrid)%refined(ind)
              if(leaf)nleaf=nleaf+1
           end do
        end do
     endif

  end do

  if(nleaf>0)then
     write(*,*)'================================================'
     write(*,114)nleaf

     allocate(ll(1:nleaf))
     allocate(xx(1:nleaf))
     allocate(dd(1:nleaf))
     allocate(uu(1:nleaf))
     allocate(pp(1:nleaf))
     allocate(ii(1:nleaf))

     nleaf=0
     do ilevel=r%levelmin,r%nlevelmax

        if(m%noct_tot(ilevel)>0)then
           do igrid=m%head(ilevel),m%tail(ilevel)
              do ind=1,2
                 leaf = .not. m%grid(igrid)%refined(ind)
                 if(leaf)then
                    nleaf=nleaf+1
                    ll(nleaf)=m%grid(igrid)%lev
                    xx(nleaf)=(2*(m%grid(igrid)%ckey(1)-m%box_ckey_min(1,ilevel))+ind-0.5)/(2.*m%ckey_max(ilevel))*r%boxlen
                    dd(nleaf)=m%grid(igrid)%uold(ind,1)
                    uu(nleaf)=m%grid(igrid)%uold(ind,2)/m%grid(igrid)%uold(ind,1)
                    pp(nleaf)=(r%gamma-1)*(m%grid(igrid)%uold(ind,3)-0.5*m%grid(igrid)%uold(ind,2)**2/m%grid(igrid)%uold(ind,1))
                 endif
              end do
           end do
        endif
     end do
     write(*,*)'lev      x           d          u          P'
     call quick_sort_dp(xx,ii,nleaf)
     do i=1,nleaf
        write(*,113) &
             & ll(ii(i)), &
             & xx(i), &
             & dd(ii(i)), &
             & uu(ii(i)), &
             & pp(ii(i))
     end do
     write(*,*)'================================================'

     deallocate(ll)
     deallocate(xx)
     deallocate(dd)
     deallocate(uu)
     deallocate(pp)
     deallocate(ii)

  endif

113 format(i3,1x,1pe12.5,1x,3(1pe10.3,1x))
114 format(' Output ',i5,' cells')

#endif

end subroutine write_screen
