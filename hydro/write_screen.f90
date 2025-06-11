subroutine write_screen(r,m)
  use amr_parameters, only: ndim
  use amr_commons, only: mesh_t, run_t
  use hydro_parameters, only: nener
  implicit none
  type(run_t)::r
  type(mesh_t)::m

  ! Local variables
  integer::ilevel,nleaf,i
  integer,dimension(:),allocatable::ll,ii
  real(kind=8),dimension(:),allocatable::xx,dd,uu,vv,ww,pp
#ifdef MHD
  real(kind=8),dimension(:),allocatable::AA,BB,CC
#endif
#if NENER>0
  integer::irad
  real(kind=8),dimension(:,:),allocatable::prad
#endif
  integer::igrid,ind
  logical::leaf
  real(kind=8)::ekin,emag,erad
  
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
     allocate(ll(1:nleaf))
     allocate(xx(1:nleaf))
     allocate(dd(1:nleaf))
     allocate(uu(1:nleaf))
     allocate(vv(1:nleaf))
     allocate(ww(1:nleaf))
#if NENER>0
     allocate(prad(1:nener,1:nleaf))
#endif
#ifdef MHD
     allocate(AA(1:nleaf))
     allocate(BB(1:nleaf))
     allocate(CC(1:nleaf))
#endif
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
                    vv(nleaf)=m%grid(igrid)%uold(ind,3)/m%grid(igrid)%uold(ind,1)
                    ww(nleaf)=m%grid(igrid)%uold(ind,4)/m%grid(igrid)%uold(ind,1)
                    ekin=0.5*dd(nleaf)*(uu(nleaf)**2+vv(nleaf)**2+ww(nleaf)**2)
                    emag=0.0
#ifdef MHD
                    AA(nleaf)=0.5*(m%grid(igrid)%bold(ind,1)+m%grid(igrid)%bold(ind,4))
                    BB(nleaf)=0.5*(m%grid(igrid)%bold(ind,2)+m%grid(igrid)%bold(ind,5))
                    CC(nleaf)=0.5*(m%grid(igrid)%bold(ind,3)+m%grid(igrid)%bold(ind,6))
                    emag=0.5*(AA(nleaf)**2+BB(nleaf)**2+CC(nleaf)**2)
#endif
                    erad=0.0d0
#if NENER>0
                    do irad=1,nener
                       erad=erad+m%grid(igrid)%uold(ind,5+irad)
                       prad(irad,nleaf)=(r%gamma_rad(irad)-1)*m%grid(igrid)%uold(ind,5+irad)
                    end do
#endif
                    pp(nleaf)=m%grid(igrid)%uold(ind,5)-ekin-emag-erad
                    pp(nleaf)=(r%gamma-1)*pp(nleaf)

                    if(ABS(uu(nleaf))<1d-99)uu(nleaf)=0.0D0
                    if(ABS(vv(nleaf))<1d-99)vv(nleaf)=0.0D0
                    if(ABS(ww(nleaf))<1d-99)ww(nleaf)=0.0D0
#ifdef MHD
                    if(ABS(AA(nleaf))<1d-99)uu(nleaf)=0.0D0
                    if(ABS(BB(nleaf))<1d-99)vv(nleaf)=0.0D0
                    if(ABS(CC(nleaf))<1d-99)ww(nleaf)=0.0D0
#endif
                 endif
              end do
           end do
        endif
     end do
     call quick_sort_dp(xx,ii,nleaf)

#ifdef MHD
     write(*,*)'================================================'
     write(*,114)nleaf
     write(*,*)'lev      x           d          u          v          w          P          A          B          C'
     do i=1,nleaf
        write(*,115)      &
             & ll(ii(i)), &
             & xx(i),     &
             & dd(ii(i)), &
             & uu(ii(i)), &
             & vv(ii(i)), &
             & ww(ii(i)), &
             & pp(ii(i)), &
             & AA(ii(i)), &
             & BB(ii(i)), &
             & CC(ii(i))
     end do
     write(*,*)'================================================'
#elif NENER>0
     write(*,*)'================================================'
     write(*,114)nleaf
     write(*,*)'lev      x           d          u          P        P_rad'
     do i=1,nleaf
        write(*,116)      &
             & ll(ii(i)), &
             & xx(i),     &
             & dd(ii(i)), &
             & uu(ii(i)), &
             & pp(ii(i)), &
             & prad(1,ii(i))
     end do
     write(*,*)'================================================'
#else
     write(*,*)'================================================'
     write(*,114)nleaf
     write(*,*)'lev      x           d          u          P'
     do i=1,nleaf
        write(*,113)      &
             & ll(ii(i)), &
             & xx(i),     &
             & dd(ii(i)), &
             & uu(ii(i)), &
             & pp(ii(i))
     end do
     write(*,*)'================================================'
#endif
     deallocate(ll)
     deallocate(xx)
     deallocate(dd)
     deallocate(uu)
     deallocate(vv)
     deallocate(ww)
     deallocate(pp)
     deallocate(ii)
#if NENER>0
     deallocate(prad)
#endif
#ifdef MHD
     deallocate(AA)
     deallocate(BB)
     deallocate(CC)
#endif

  endif

113 format(i3,1x,1pe12.5,1x,3(1pe10.3,1x))
114 format(' Output ',i5,' cells')
115 format(i3,1x,1pe12.5,1x,8(1pe10.3,1x))
116 format(i3,1x,1pe12.5,1x,4(1pe10.3,1x))

#endif

end subroutine write_screen
