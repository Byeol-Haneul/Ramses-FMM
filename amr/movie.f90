!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine output_frame_2(r,g,m,p)
  use amr_parameters, only: dp,ndim,nvector,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  ! Local variables
  integer::info
  integer,parameter::tag=100

  character(len=5) :: istep_str
  character(len=100) :: moviedir, moviecmd, infofile
  character(len=100),dimension(0:nvar+2) :: moviefiles
  integer::nlevelmax_frame,nstride
  integer::ilun,ind
  integer::imin,imax,jmin,jmax,ii,jj,kk
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::xcen,ycen,zcen,delx,dely,delz
  real(dp)::xleft_frame,xright_frame,yleft_frame,yright_frame,zleft_frame,zright_frame
  real(dp)::xleft,xright,yleft,yright,zleft,zright
  real(dp)::xxleft,xxright,yyleft,yyright
  real(dp)::dx_frame,dy_frame,dx,dx_loc,dx_min
  real(dp)::dx_cell,dy_cell,dz_cell,dvol
  logical::ok
  real(dp),dimension(1:ndim)::xx
  real(kind=8),dimension(:,:,:),allocatable::data_frame,data_frame_all
  real(kind=8),dimension(:,:),allocatable::dens,dens_all,vol,vol_all
  real(kind=4),dimension(:,:),allocatable::data_single
  real(kind=8)::temp,ekk
  integer::igrid,idim,ilevel
  integer::proj_ind,nh_temp,nw_temp
  real(kind=4)::ratio

  logical::opened

  character(len=1)::temp_string

  nh_temp = r%nh_frame
  nw_temp = r%nw_frame
  
 do proj_ind=1,LEN(trim(r%proj_axis)) 
  opened=.false.

#if NDIM > 1
  if(r%imov<1)r%imov=1
  if(r%imov>r%imovout)return

  ! Determine the filename, dir, etc
  if(g%myid==1)write(*,*)'Computing and dumping movie frame'

  call title(r%imov, istep_str)
  write(temp_string,'(I1)') proj_ind
  moviedir = 'movie_threadsafe'//trim(temp_string)//'/'
  moviecmd = 'mkdir -p '//trim(moviedir)
  if(g%myid==1) write(*,*) "Writing frame ", istep_str
  if(.not.g%withoutmkdir) then 
#ifdef NOSYSTEM
     if(g%myid==1)call PXFMKDIR(TRIM(moviedir),LEN(TRIM(moviedir)),O'755',info)  
#else
     if(g%myid==1)call system(moviecmd)
#endif
  endif
  
  infofile = trim(moviedir)//'info_'//trim(istep_str)//'.txt'
  if(g%myid==1)call output_info_2(r,g,infofile)
  
  moviefiles(0) = trim(moviedir)//'temp_'//trim(istep_str)//'.map'
  moviefiles(1) = trim(moviedir)//'dens_'//trim(istep_str)//'.map'
  moviefiles(2) = trim(moviedir)//'vx_'//trim(istep_str)//'.map'
  moviefiles(3) = trim(moviedir)//'vy_'//trim(istep_str)//'.map'
#if NDIM>2
  moviefiles(4) = trim(moviedir)//'vz_'//trim(istep_str)//'.map'
#endif
#if NDIM==2
  moviefiles(4) = trim(moviedir)//'pres_'//trim(istep_str)//'.map'
#endif
#if NDIM>2
  moviefiles(5) = trim(moviedir)//'pres_'//trim(istep_str)//'.map'
#endif
#if NVAR>5
  do ll=6,NVAR
    write(dummy,'(I3.1)') ll
    moviefiles(ll) = trim(moviedir)//'var'//trim(adjustl(dummy))//'_'//trim(istep_str)//'.map'
 end do
#endif
  moviefiles(NVAR+1) = trim(moviedir)//'dm_'//trim(istep_str)//'.map'
  moviefiles(NVAR+2) = trim(moviedir)//'stars_'//trim(istep_str)//'.map'

  if(r%levelmax_frame==0)then
     nlevelmax_frame=r%nlevelmax
  else if (r%levelmax_frame.gt.r%nlevelmax)then
     nlevelmax_frame=r%nlevelmax
  else
     nlevelmax_frame=r%levelmax_frame
  endif

  ! Conversion factor from user units to cgs units
  call units_2(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Compute frame boundaries
  if(r%proj_axis(proj_ind:proj_ind).eq.'x')then
     xcen=r%ycentre_frame(proj_ind*4-3)+r%ycentre_frame(proj_ind*4-2)*g%aexp+r%ycentre_frame(proj_ind*4-1)*g%aexp**2+r%ycentre_frame(proj_ind*4)*g%aexp**3
     ycen=r%zcentre_frame(proj_ind*4-3)+r%zcentre_frame(proj_ind*4-2)*g%aexp+r%zcentre_frame(proj_ind*4-1)*g%aexp**2+r%zcentre_frame(proj_ind*4)*g%aexp**3
     zcen=r%xcentre_frame(proj_ind*4-3)+r%xcentre_frame(proj_ind*4-2)*g%aexp+r%xcentre_frame(proj_ind*4-1)*g%aexp**2+r%xcentre_frame(proj_ind*4)*g%aexp**3
     delx=min(2*min(xcen,ycen,zcen,r%boxlen-xcen,r%boxlen-ycen,r%boxlen-zcen),r%deltay_frame(proj_ind*2-1)+r%deltay_frame(proj_ind*2)/g%aexp) !+r%deltax_frame(3)*g%aexp**2+r%deltax_frame(4)*g%aexp**3  !Essentially comoving or physical
     dely=min(2*min(xcen,ycen,zcen,r%boxlen-xcen,r%boxlen-ycen,r%boxlen-zcen),r%deltaz_frame(proj_ind*2-1)+r%deltaz_frame(proj_ind*2)/g%aexp) !+r%deltay_frame(3)*g%aexp**2+r%deltay_frame(4)*g%aexp**3
     delz=min(2*min(xcen,ycen,zcen,r%boxlen-xcen,r%boxlen-ycen,r%boxlen-zcen),r%deltax_frame(proj_ind*2-1)+r%deltax_frame(proj_ind*2)/g%aexp) !+r%deltaz_frame(3)*g%aexp**2+r%deltaz_frame(4)*g%aexp**3
  elseif(r%proj_axis(proj_ind:proj_ind).eq.'y')then
     xcen=r%xcentre_frame(proj_ind*4-3)+r%xcentre_frame(proj_ind*4-2)*g%aexp+r%xcentre_frame(proj_ind*4-1)*g%aexp**2+r%xcentre_frame(proj_ind*4)*g%aexp**3
     ycen=r%zcentre_frame(proj_ind*4-3)+r%zcentre_frame(proj_ind*4-2)*g%aexp+r%zcentre_frame(proj_ind*4-1)*g%aexp**2+r%zcentre_frame(proj_ind*4)*g%aexp**3
     zcen=r%ycentre_frame(proj_ind*4-3)+r%ycentre_frame(proj_ind*4-2)*g%aexp+r%ycentre_frame(proj_ind*4-1)*g%aexp**2+r%ycentre_frame(proj_ind*4)*g%aexp**3
     delx=min(2*min(xcen,ycen,zcen,r%boxlen-xcen,r%boxlen-ycen,r%boxlen-zcen),r%deltax_frame(proj_ind*2-1)+r%deltax_frame(proj_ind*2)/g%aexp) !+r%deltax_frame(3)*g%aexp**2+r%deltax_frame(4)*g%aexp**3  !Essentially comoving or physical
     dely=min(2*min(xcen,ycen,zcen,r%boxlen-xcen,r%boxlen-ycen,r%boxlen-zcen),r%deltaz_frame(proj_ind*2-1)+r%deltaz_frame(proj_ind*2)/g%aexp) !+r%deltay_frame(3)*g%aexp**2+r%deltay_frame(4)*g%aexp**3
     delz=min(2*min(xcen,ycen,zcen,r%boxlen-xcen,r%boxlen-ycen,r%boxlen-zcen),r%deltay_frame(proj_ind*2-1)+r%deltay_frame(proj_ind*2)/g%aexp) !+r%deltaz_frame(3)*g%aexp**2+r%deltaz_frame(4)*g%aexp**3
  else
     xcen=r%xcentre_frame(proj_ind*4-3)+r%xcentre_frame(proj_ind*4-2)*g%aexp+r%xcentre_frame(proj_ind*4-1)*g%aexp**2+r%xcentre_frame(proj_ind*4)*g%aexp**3
     ycen=r%ycentre_frame(proj_ind*4-3)+r%ycentre_frame(proj_ind*4-2)*g%aexp+r%ycentre_frame(proj_ind*4-1)*g%aexp**2+r%ycentre_frame(proj_ind*4)*g%aexp**3
     zcen=r%zcentre_frame(proj_ind*4-3)+r%zcentre_frame(proj_ind*4-2)*g%aexp+r%zcentre_frame(proj_ind*4-1)*g%aexp**2+r%zcentre_frame(proj_ind*4)*g%aexp**3
     delx=min(2*min(xcen,ycen,zcen,r%boxlen-xcen,r%boxlen-ycen,r%boxlen-zcen),r%deltax_frame(proj_ind*2-1)+r%deltax_frame(proj_ind*2)/g%aexp) !+r%deltax_frame(3)*g%aexp**2+r%deltax_frame(4)*g%aexp**3  !Essentially comoving or physical
     dely=min(2*min(xcen,ycen,zcen,r%boxlen-xcen,r%boxlen-ycen,r%boxlen-zcen),r%deltay_frame(proj_ind*2-1)+r%deltay_frame(proj_ind*2)/g%aexp) !+r%deltay_frame(3)*g%aexp**2+r%deltay_frame(4)*g%aexp**3
     delz=min(2*min(xcen,ycen,zcen,r%boxlen-xcen,r%boxlen-ycen,r%boxlen-zcen),r%deltaz_frame(proj_ind*2-1)+r%deltaz_frame(proj_ind*2)/g%aexp) !+r%deltaz_frame(3)*g%aexp**2+r%deltaz_frame(4)*g%aexp**3
  endif
  
  ratio = delx/dely
  if(ratio.gt.1)then
     r%nw_frame=nh_temp*ratio
  else
     r%nh_frame=nw_temp/ratio
  endif
  
  xleft_frame=xcen-delx/2.
  xright_frame=xcen+delx/2.
  yleft_frame=ycen-dely/2.
  yright_frame=ycen+dely/2.
  zleft_frame=zcen-delz/2.
  zright_frame=zcen+delz/2.
  
  ! Allocate image
  allocate(data_frame(1:r%nw_frame,1:r%nh_frame,0:NVAR+2))
  allocate(dens(1:r%nw_frame,1:r%nh_frame))
  allocate(vol(1:r%nw_frame,1:r%nh_frame))
  data_frame=0d0
  dens=0d0
  vol=0d0
  dx_frame=delx/dble(r%nw_frame)
  dy_frame=dely/dble(r%nh_frame)

  if(r%hydro) then
     ! Loop over levels
     do ilevel=r%levelmin,nlevelmax_frame
        
        ! Mesh size at level ilevel in coarse cell units
        dx=0.5D0**ilevel
        
        dx_loc=dx*r%boxlen
        dx_min=0.5D0**r%nlevelmax*r%boxlen
        
        ! Loop over grids by vector sweeps
        do igrid=m%head(ilevel),m%tail(ilevel)
           
           ! Loop over cells
           do ind=1,twotondim
              
              ! Compute cell centre position in code units
              do idim=1,ndim
                 nstride=2**(idim-1)
                 xx(idim)=(2*m%grid(igrid)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx_loc
              end do
              
              ! Check if cell is to be considered
              ok=(.NOT.m%grid(igrid)%refined(ind)).or.(ilevel==nlevelmax_frame)
   
              if(ok)then
                 ! Check if the cell intersect the domain
#if NDIM>2                 
                 if(r%proj_axis(proj_ind:proj_ind).eq.'x')then
                    xleft=xx(2)-dx_loc/2.
                    xright=xx(2)+dx_loc/2.
                    yleft=xx(3)-dx_loc/2.
                    yright=xx(3)+dx_loc/2.
                 elseif(r%proj_axis(proj_ind:proj_ind).eq.'y')then
                    xleft=xx(1)-dx_loc/2.
                    xright=xx(1)+dx_loc/2.
                    yleft=xx(3)-dx_loc/2.
                    yright=xx(3)+dx_loc/2.
                 else
                    xleft=xx(1)-dx_loc/2.
                    xright=xx(1)+dx_loc/2.
                    yleft=xx(2)-dx_loc/2.
                    yright=xx(2)+dx_loc/2.
                 endif
                 
                 if(r%proj_axis(proj_ind:proj_ind).eq.'x')then
                    zleft=xx(1)-dx_loc/2.
                    zright=xx(1)+dx_loc/2.
                 elseif(r%proj_axis(proj_ind:proj_ind).eq.'y')then
                    zleft=xx(2)-dx_loc/2.
                    zright=xx(2)+dx_loc/2.
                 else
                    zleft=xx(3)-dx_loc/2.
                    zright=xx(3)+dx_loc/2.
                 endif
                 if(    xright.lt.xleft_frame.or.xleft.ge.xright_frame.or.&
                      & yright.lt.yleft_frame.or.yleft.ge.yright_frame.or.&
                      & zright.lt.zleft_frame.or.zleft.ge.zright_frame)cycle
#else
                 xleft=xx(1)-dx_loc/2.
                 xright=xx(1)+dx_loc/2.
                 yleft=xx(2)-dx_loc/2.
                 yright=xx(2)+dx_loc/2.
                 
                 if(    xright.lt.xleft_frame.or.xleft.ge.xright_frame.or.&
                      & yright.lt.yleft_frame.or.yleft.ge.yright_frame)cycle
#endif
                 ! Compute map indices for the cell
                 if(xleft>xleft_frame)then
                    imin=min(int((xleft-xleft_frame)/dx_frame)+1,r%nw_frame)
                 else
                    imin=1
                 endif
                 imax=min(int((xright-xleft_frame)/dx_frame)+1,r%nw_frame)
                 if(yleft>yleft_frame)then
                    jmin=min(int((yleft-yleft_frame)/dy_frame)+1,r%nh_frame) ! change
                 else
                    jmin=1
                 endif
                 jmax=min(int((yright-yleft_frame)/dy_frame)+1,r%nh_frame) ! change
                 
                 ! Fill up map with projected mass
#if NDIM>2                 
                 dz_cell=min(zright_frame,zright)-max(zleft_frame,zleft) ! change
#endif
                 do ii=imin,imax
                    xxleft=xleft_frame+dble(ii-1)*dx_frame
                    xxright=xxleft+dx_frame
                    dx_cell=min(xxright,xright)-max(xxleft,xleft)
                    do jj=jmin,jmax
                       yyleft=yleft_frame+dble(jj-1)*dy_frame
                       yyright=yyleft+dy_frame
                       dy_cell=min(yyright,yright)-max(yyleft,yleft)
                       ! Intersection volume
                       dvol=dx_cell*dy_cell
#if NDIM>2                 
                       dvol=dvol*dz_cell
#endif
#ifdef HYDRO
                       dens(ii,jj)=dens(ii,jj)+dvol*max(m%grid(igrid)%uold(ind,1),r%smallr)
                       vol(ii,jj)=vol(ii,jj)+dvol                       
                       data_frame(ii,jj,1)=data_frame(ii,jj,1)+dvol*max(m%grid(igrid)%uold(ind,1),r%smallr)**2
                       do kk=2,NVAR
                          if(r%movie_vars(kk).eq.1) data_frame(ii,jj,kk)=data_frame(ii,jj,kk)+dvol*m%grid(igrid)%uold(ind,kk)
                       end do
#endif
                       if (r%movie_vars(0).eq.1)then
#ifdef HYDRO
                          ! Get temperature
                          ekk=0.0d0
                          do idim=1,3
                             ekk=ekk+0.5*m%grid(igrid)%uold(ind,idim+1)**2/max(m%grid(igrid)%uold(ind,1),r%smallr)
                          enddo
                          temp=(r%gamma-1.0)*(m%grid(igrid)%uold(ind,ndim+2)-ekk) !pressure
                          temp=max(temp/max(m%grid(igrid)%uold(ind,1),r%smallr),r%smallc**2)*scale_T2 !temperature in K                          
                          data_frame(ii,jj,0)=data_frame(ii,jj,0)+dvol*max(m%grid(igrid)%uold(ind,1),r%smallr)*temp !mass weighted temperature
#endif
                       end if
                    end do
                 end do
              end if
   
           end do
           ! End loop over cells
   
        end do
        ! End loop over grids
     end do
  ! End loop over levels
  end if

#ifndef WITHOUTMPI
  allocate(data_frame_all(1:r%nw_frame,1:r%nh_frame,0:NVAR+2))
  call MPI_ALLREDUCE(data_frame,data_frame_all,r%nw_frame*r%nh_frame*(NVAR+2+1),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  allocate(dens_all(1:r%nw_frame,1:r%nh_frame))
  call MPI_ALLREDUCE(dens,dens_all,r%nw_frame*r%nh_frame,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  allocate(vol_all(1:r%nw_frame,1:r%nh_frame))
  call MPI_ALLREDUCE(vol,vol_all,r%nw_frame*r%nh_frame,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  data_frame=data_frame_all
  dens=dens_all
  vol=vol_all
  deallocate(data_frame_all)
  deallocate(dens_all)
  deallocate(vol_all)
#endif
  ! Convert into mass weighted                                                                                                         
  do ii=1,r%nw_frame
    do jj=1,r%nh_frame
      do kk=0,NVAR
        if(r%movie_vars(kk).eq.1) data_frame(ii,jj,kk)=data_frame(ii,jj,kk)/dens(ii,jj)
      end do
    end do
  end do
  deallocate(dens)
  deallocate(vol)

  if(g%myid==1)then
     ilun=10
     allocate(data_single(1:r%nw_frame,1:r%nh_frame))
     ! Output mass weighted density
     do kk=0, NVAR+2
       if (r%movie_vars(kk).eq.1)then
         open(ilun,file=TRIM(moviefiles(kk)),form='unformatted')
         data_single=data_frame(:,:,kk)
         rewind(ilun)  
         if(r%tendmov>0)then
            write(ilun)g%t,delx,dely,delz
         else
            write(ilun)g%aexp,delx,dely,delz
         endif
         write(ilun)r%nw_frame,r%nh_frame
         write(ilun)data_single
         close(ilun)
       end if
     end do
     
     deallocate(data_single)
  endif
  
  deallocate(data_frame)
#endif

  ! Update counter
  if(proj_ind.eq.len(trim(r%proj_axis))) then 
     ! Increase counter and skip frames if timestep is large
     r%imov=r%imov+1
     do while((r%amovout(r%imov)<g%aexp.or.r%tmovout(r%imov)<g%t).and.(r%imov.lt.r%imovout))
        r%imov=r%imov+1
     end do
  endif

  r%nw_frame = nw_temp
  r%nh_frame = nh_temp

 enddo

end subroutine output_frame_2
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine set_movie_vars_2(r)
  use amr_commons, only: run_t
  ! This routine sets the movie vars from textual form
  type(run_t)::r
  
#if NVAR>5
  integer::ll
  character(LEN=5)::dummy
#endif

  if(ANY(r%movie_vars_txt=='temp '))r%movie_vars(0)=1
  if(ANY(r%movie_vars_txt=='dens '))r%movie_vars(1)=1
  if(ANY(r%movie_vars_txt=='vx   '))r%movie_vars(2)=1
  if(ANY(r%movie_vars_txt=='vy   '))r%movie_vars(3)=1
#if NDIM>2
  if(ANY(r%movie_vars_txt=='vz   '))r%movie_vars(4)=1
#endif
#if NDIM==2
  if(ANY(r%movie_vars_txt=='pres '))r%movie_vars(4)=1
#endif
#if NDIM>2
  if(ANY(r%movie_vars_txt=='pres '))r%movie_vars(5)=1
#endif
#if NVAR>5
  do ll=6,NVAR
    write(dummy,'(I3.1)') ll
    if(ANY(r%movie_vars_txt=='var'//trim(adjustl(dummy))//' '))r%movie_vars(ll)=1
 end do
#endif
#ifdef SOLVERmhd
  if(ANY(r%movie_vars_txt=='bxl  '))r%movie_vars(6)=1
  if(ANY(r%movie_vars_txt=='byl  '))r%movie_vars(7)=1
  if(ANY(r%movie_vars_txt=='bzl  '))r%movie_vars(8)=1
  if(ANY(r%movie_vars_txt=='bxr  '))r%movie_vars(NVAR+1)=1
  if(ANY(r%movie_vars_txt=='byr  '))r%movie_vars(NVAR+2)=1
  if(ANY(r%movie_vars_txt=='bzr  '))r%movie_vars(NVAR+3)=1
  if(ANY(r%movie_vars_txt=='pmag '))r%movie_vars(NVAR+4)=1
  if(ANY(r%movie_vars_txt=='dm   '))r%movie_vars(NVAR+5)=1
  if(ANY(r%movie_vars_txt=='stars'))r%movie_vars(NVAR+6)=1
#else
  if(ANY(r%movie_vars_txt=='dm   '))r%movie_vars(NVAR+1)=1
  if(ANY(r%movie_vars_txt=='stars'))r%movie_vars(NVAR+2)=1
#endif
end subroutine set_movie_vars_2
