subroutine write_screen(r,m)
  use amr_parameters, only: dp,ndim
  use amr_commons, only: mesh_t, run_t
  implicit none
  type(run_t)::r
  type(mesh_t)::m

  ! Local variables
  integer::ilevel
  integer::igrid,ind
  logical::leaf
  
#ifdef HYDRO

  if(ndim>1)return
  
  write(*,*)'==============================================='
  do ilevel=r%levelmin,r%nlevelmax
     
     if(m%noct_tot(ilevel)>0)then
!        write(*,*)'================================'
!        write(*,'(" Level",I4," has ",I4, " grids.")')ilevel,m%noct(ilevel)
        do igrid=m%head(ilevel),m%tail(ilevel)
           do ind=1,2
!!$              write(*,'(I4,1X,I8,1X,I4,1X,I4,1X,I4,1X,L,1X,2(I4,1X),6(1PE10.3))')&
!!$                   & igrid,m%grid(igrid)%ckey(1),ind,m%grid(igrid)%flag1(ind),&
!!$                   & m%grid(igrid)%lev,m%grid(igrid)%refined(ind),&
!!$                   & m%grid(igrid)%hkey(1),m%grid(igrid)%superoct, &
!!$                   & (2*m%grid(igrid)%ckey(1)+ind-0.5)/(2.*m%ckey_max(ilevel)),&
!!$                   & m%grid(igrid)%uold(ind,1),m%grid(igrid)%uold(ind,2)/m%grid(igrid)%uold(ind,1),&
!!$                   & (r%gamma-1)*(m%grid(igrid)%uold(ind,3)-0.5*m%grid(igrid)%uold(ind,2)**2/m%grid(igrid)%uold(ind,1))

              leaf = .not. m%grid(igrid)%refined(ind)
              if(leaf)then
                 write(*,'(1(I4,1X),4(1PE10.3,1X))')&
                      & m%grid(igrid)%lev,&
                      & (2*m%grid(igrid)%ckey(1)+ind-0.5)/(2.*m%ckey_max(ilevel))*r%boxlen,&
                      & m%grid(igrid)%uold(ind,1),m%grid(igrid)%uold(ind,2)/m%grid(igrid)%uold(ind,1),&
                      & (r%gamma-1)*(m%grid(igrid)%uold(ind,3)-0.5*m%grid(igrid)%uold(ind,2)**2/m%grid(igrid)%uold(ind,1))
              endif
              
           end do
        end do
!        write(*,*)'================================'
     endif

  end do
  write(*,*)'==============================================='

#endif

end subroutine write_screen
