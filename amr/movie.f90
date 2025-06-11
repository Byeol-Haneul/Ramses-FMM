module movie_module
contains
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine m_output_frame(pst)
  use amr_parameters, only: ndim, nvector, twotondim, flen
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtgrp
  use ramses_commons, only: pst_t
  use output_amr_module, only: output_info
  use mdl_module, only: mdl_mkdir
  implicit none
  type(pst_t)::pst

  ! Local variables
  character(len=1)::temp_string
  character(len=5)::istep_str
  character(len=flen)::moviedir,moviecmd,infofile
  character(len=flen),dimension(0:nvar+2+nrtgrp)::moviefiles
  integer::ilun,info,kk,ind_proj,input_size,output_size
  integer,dimension(:),allocatable::input_array,output_array
  real(kind=8)::delx,dely,delz
  real(kind=8),dimension(:),allocatable::data_frame
  real(kind=8),dimension(:),allocatable::dens
  real(kind=4),dimension(:),allocatable::data_single
  integer::ll
  character(LEN=5)::dummy

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)

  do ind_proj=1,LEN(trim(r%proj_axis)) 
     
#if NDIM > 1
     if(r%imov<1)r%imov=1
     if(r%imov>r%imovout)return
     
     ! Determine the filename, dir, etc
     write(*,*)'Computing and dumping movie frame'
     
     call title(r%imov, istep_str)
     write(temp_string,'(I1)') ind_proj
     moviedir = 'movie'//trim(temp_string)//'/'
!     moviecmd = 'mkdir -p '//trim(moviedir)
     write(*,*) "Writing frame ", istep_str
     if(.not.g%withoutmkdir) then 
       call mdl_mkdir(mdl,moviedir)
! #ifdef NOSYSTEM
!         call PXFMKDIR(TRIM(moviedir),LEN(TRIM(moviedir)),O'755',info)
! #else
!         call system(moviecmd)
! #endif
     endif
     
     infofile = trim(moviedir)//'info_'//trim(istep_str)//'.txt'
     call output_info(r,g,infofile)
     
     moviefiles(1) = trim(moviedir)//'dens_'//trim(istep_str)//'.map'
     moviefiles(2) = trim(moviedir)//'vx_'//trim(istep_str)//'.map'
     moviefiles(3) = trim(moviedir)//'vy_'//trim(istep_str)//'.map'
     moviefiles(4) = trim(moviedir)//'vz_'//trim(istep_str)//'.map'
     moviefiles(5) = trim(moviedir)//'temp_'//trim(istep_str)//'.map'
#if NVAR>5
     do ll=6,NVAR
        write(dummy,'(I3.1)') ll
        moviefiles(ll) = trim(moviedir)//'var'//trim(adjustl(dummy))//'_'//trim(istep_str)//'.map'
     end do
#endif
     do ll=nvar+1,nvar+nrtgrp
        write(dummy,'(I3.1)') ll-nvar
        moviefiles(ll) = trim(moviedir)//'Fp'//trim(adjustl(dummy))//'_'//trim(istep_str)//'.map'
     end do
     moviefiles(NVAR+nrtgrp+1) = trim(moviedir)//'dm_'//trim(istep_str)//'.map'
     moviefiles(NVAR+nrtgrp+2) = trim(moviedir)//'stars_'//trim(istep_str)//'.map'
          
     ! Allocate image
     allocate(data_single(1:r%nw_frame*r%nh_frame))
     allocate(data_frame (1:r%nw_frame*r%nh_frame))
     allocate(dens       (1:r%nw_frame*r%nh_frame))
     input_size=2
     output_size=2*r%nw_frame*r%nh_frame
     allocate(input_array(1:input_size))
     allocate(output_array(1:output_size))     
     input_array(1)=ind_proj

     ! Compute image size
     if(r%proj_axis(ind_proj:ind_proj).eq.'x')then
        delx=r%deltay_frame(ind_proj*2-1)+r%deltay_frame(ind_proj*2)/g%aexp
        dely=r%deltaz_frame(ind_proj*2-1)+r%deltaz_frame(ind_proj*2)/g%aexp
        delz=r%deltax_frame(ind_proj*2-1)+r%deltax_frame(ind_proj*2)/g%aexp
     elseif(r%proj_axis(ind_proj:ind_proj).eq.'y')then
        delx=r%deltax_frame(ind_proj*2-1)+r%deltax_frame(ind_proj*2)/g%aexp
        dely=r%deltaz_frame(ind_proj*2-1)+r%deltaz_frame(ind_proj*2)/g%aexp
        delz=r%deltay_frame(ind_proj*2-1)+r%deltay_frame(ind_proj*2)/g%aexp
     else
        delx=r%deltax_frame(ind_proj*2-1)+r%deltax_frame(ind_proj*2)/g%aexp
        dely=r%deltay_frame(ind_proj*2-1)+r%deltay_frame(ind_proj*2)/g%aexp
        delz=r%deltaz_frame(ind_proj*2-1)+r%deltaz_frame(ind_proj*2)/g%aexp
     endif

     if(r%hydro) then
        
        ! Compute column density
        input_array(2)=0
        call r_output_frame(pst,input_array,input_size,output_array,output_size)
        dens=transfer(output_array,dens)
        
        do kk=1,NVAR+nrtgrp
           if(r%movie_vars(kk).eq.1)then
              ! Compute mass-weighted projected quantities
              input_array(2)=kk
              call r_output_frame(pst,input_array,input_size,output_array,output_size)
              data_frame=transfer(output_array,data_frame)
              ! Divide by column density
              data_frame=data_frame/dens
              ! Store image to file
              ilun=10
              open(ilun,file=TRIM(moviefiles(kk)),form='unformatted')
              data_single=real(data_frame,kind=4)
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
          
     endif

     deallocate(data_single)
     deallocate(data_frame)
     deallocate(dens)
     deallocate(input_array)
     deallocate(output_array)
     
     ! Update counter
     if(ind_proj.eq.len(trim(r%proj_axis))) then 
        ! Increase counter and skip frames if timestep is large
        r%imov=r%imov+1
        if(r%aendmov>0)then
           do while((r%aendmov*dble(r%imov)/dble(r%imovout)<g%aexp).and.(r%imov.lt.r%imovout))
              r%imov=r%imov+1
           end do
        end if
        if(r%tendmov>0)then
           do while((r%tendmov*dble(r%imov)/dble(r%imovout)<g%t).and.(r%imov.lt.r%imovout))
              r%imov=r%imov+1
           end do
        end if
     endif
     
#endif

  enddo

  end associate

end subroutine m_output_frame
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
recursive subroutine r_output_frame(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  use hilbert
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer,dimension(:),allocatable::next_output_array
  
  integer::ind_proj,ind_var
  real(kind=8),dimension(:),allocatable::map,next_map
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_OUTPUT_FRAME,pst%iUpper+1,input_size,output_size,input_array)
     call r_output_frame(pst%pLower,input_array,input_size,output_array,output_size)
     allocate(next_output_array(1:output_size))
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output_array)
     allocate(map(1:pst%s%r%nw_frame*pst%s%r%nh_frame))
     allocate(next_map(1:pst%s%r%nw_frame*pst%s%r%nh_frame))
     map=transfer(output_array,map)
     next_map=transfer(next_output_array,next_map)
     map=map+next_map
     output_array=transfer(map,output_array)
     deallocate(next_output_array)
     deallocate(map)
     deallocate(next_map)
  else
     ind_proj=input_array(1)
     ind_var=input_array(2)
     allocate(map(1:pst%s%r%nw_frame*pst%s%r%nh_frame))
     map=0d0
     call output_frame(pst%s%r,pst%s%g,pst%s%m,ind_proj,ind_var,pst%s%r%nw_frame*pst%s%r%nh_frame,map)
     output_array=transfer(map,output_array)
     deallocate(map)
  endif

end subroutine r_output_frame
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine output_frame(r,g,m,ind_proj,ind_var,map_size,map)
  use amr_parameters, only:ndim, nvector, twotondim
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtgrp
  use amr_commons, only: run_t, global_t, mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ind_proj,map_size,ind_var
  real(kind=8),dimension(1:map_size)::map
  
  ! Local variables
  integer::nlevelmax_frame,nstride
  integer::ilun,ind,ind_map
  integer::imin,imax,jmin,jmax,ii,jj
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::xcen,ycen,zcen
  real(kind=8)::xleft_frame,xright_frame,yleft_frame,yright_frame,zleft_frame,zright_frame
  real(kind=8)::xleft,xright,yleft,yright,zleft,zright
  real(kind=8)::xxleft,xxright,yyleft,yyright
  real(kind=8)::delx,dely,delz
  real(kind=8)::dx_frame,dy_frame,dx
  real(kind=8)::dx_cell,dy_cell,dz_cell,dvol
  logical::ok
  real(kind=8),dimension(1:ndim)::xx
  real(kind=8)::temp,ekk
  integer::igrid,idim,ilevel

  if(r%levelmax_frame==0)then
     nlevelmax_frame=r%nlevelmax
  else if (r%levelmax_frame.gt.r%nlevelmax)then
     nlevelmax_frame=r%nlevelmax
  else
     nlevelmax_frame=r%levelmax_frame
  endif
  
  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  ! Compute frame centre
  if(r%proj_axis(ind_proj:ind_proj).eq.'x')then
     xcen=r%ycentre_frame(ind_proj*4-3)+r%ycentre_frame(ind_proj*4-2)*g%aexp+r%ycentre_frame(ind_proj*4-1)*g%aexp**2+r%ycentre_frame(ind_proj*4)*g%aexp**3
     ycen=r%zcentre_frame(ind_proj*4-3)+r%zcentre_frame(ind_proj*4-2)*g%aexp+r%zcentre_frame(ind_proj*4-1)*g%aexp**2+r%zcentre_frame(ind_proj*4)*g%aexp**3
     zcen=r%xcentre_frame(ind_proj*4-3)+r%xcentre_frame(ind_proj*4-2)*g%aexp+r%xcentre_frame(ind_proj*4-1)*g%aexp**2+r%xcentre_frame(ind_proj*4)*g%aexp**3
  elseif(r%proj_axis(ind_proj:ind_proj).eq.'y')then
     xcen=r%xcentre_frame(ind_proj*4-3)+r%xcentre_frame(ind_proj*4-2)*g%aexp+r%xcentre_frame(ind_proj*4-1)*g%aexp**2+r%xcentre_frame(ind_proj*4)*g%aexp**3
     ycen=r%zcentre_frame(ind_proj*4-3)+r%zcentre_frame(ind_proj*4-2)*g%aexp+r%zcentre_frame(ind_proj*4-1)*g%aexp**2+r%zcentre_frame(ind_proj*4)*g%aexp**3
     zcen=r%ycentre_frame(ind_proj*4-3)+r%ycentre_frame(ind_proj*4-2)*g%aexp+r%ycentre_frame(ind_proj*4-1)*g%aexp**2+r%ycentre_frame(ind_proj*4)*g%aexp**3
  else
     xcen=r%xcentre_frame(ind_proj*4-3)+r%xcentre_frame(ind_proj*4-2)*g%aexp+r%xcentre_frame(ind_proj*4-1)*g%aexp**2+r%xcentre_frame(ind_proj*4)*g%aexp**3
     ycen=r%ycentre_frame(ind_proj*4-3)+r%ycentre_frame(ind_proj*4-2)*g%aexp+r%ycentre_frame(ind_proj*4-1)*g%aexp**2+r%ycentre_frame(ind_proj*4)*g%aexp**3
     zcen=r%zcentre_frame(ind_proj*4-3)+r%zcentre_frame(ind_proj*4-2)*g%aexp+r%zcentre_frame(ind_proj*4-1)*g%aexp**2+r%zcentre_frame(ind_proj*4)*g%aexp**3
  endif

  ! Compute frame size
  if(r%proj_axis(ind_proj:ind_proj).eq.'x')then
     delx=r%deltay_frame(ind_proj*2-1)+r%deltay_frame(ind_proj*2)/g%aexp !+r%deltax_frame(3)*g%aexp**2+r%deltax_frame(4)*g%aexp**3
     dely=r%deltaz_frame(ind_proj*2-1)+r%deltaz_frame(ind_proj*2)/g%aexp !+r%deltay_frame(3)*g%aexp**2+r%deltay_frame(4)*g%aexp**3
     delz=r%deltax_frame(ind_proj*2-1)+r%deltax_frame(ind_proj*2)/g%aexp !+r%deltaz_frame(3)*g%aexp**2+r%deltaz_frame(4)*g%aexp**3
  elseif(r%proj_axis(ind_proj:ind_proj).eq.'y')then
     delx=r%deltax_frame(ind_proj*2-1)+r%deltax_frame(ind_proj*2)/g%aexp !+r%deltax_frame(3)*g%aexp**2+r%deltax_frame(4)*g%aexp**3
     dely=r%deltaz_frame(ind_proj*2-1)+r%deltaz_frame(ind_proj*2)/g%aexp !+r%deltay_frame(3)*g%aexp**2+r%deltay_frame(4)*g%aexp**3
     delz=r%deltay_frame(ind_proj*2-1)+r%deltay_frame(ind_proj*2)/g%aexp !+r%deltaz_frame(3)*g%aexp**2+r%deltaz_frame(4)*g%aexp**3
  else
     delx=r%deltax_frame(ind_proj*2-1)+r%deltax_frame(ind_proj*2)/g%aexp !+r%deltax_frame(3)*g%aexp**2+r%deltax_frame(4)*g%aexp**3
     dely=r%deltay_frame(ind_proj*2-1)+r%deltay_frame(ind_proj*2)/g%aexp !+r%deltay_frame(3)*g%aexp**2+r%deltay_frame(4)*g%aexp**3
     delz=r%deltaz_frame(ind_proj*2-1)+r%deltaz_frame(ind_proj*2)/g%aexp !+r%deltaz_frame(3)*g%aexp**2+r%deltaz_frame(4)*g%aexp**3
  endif

  ! Size of pixel
  dx_frame=delx/dble(r%nw_frame)
  dy_frame=dely/dble(r%nh_frame)

  ! Careful with box boundaries
  xcen=min(max(xcen,delx/2.0),r%boxlen-delx/2.0)
  ycen=min(max(ycen,dely/2.0),r%boxlen-dely/2.0)
  zcen=min(max(zcen,delz/2.0),r%boxlen-delz/2.0)

  xleft_frame =xcen-delx/2.0
  xright_frame=xcen+delx/2.0
  yleft_frame =ycen-dely/2.0
  yright_frame=ycen+dely/2.0
  zleft_frame =zcen-delz/2.0
  zright_frame=zcen+delz/2.0
  
  ! Loop over levels
  do ilevel=r%levelmin,nlevelmax_frame
     
     ! Mesh size at level ilevel in coarse cell units
     dx=r%boxlen/2**ilevel
     
     ! Loop over grids by vector sweeps
     do igrid=m%head(ilevel),m%tail(ilevel)
        
        ! Loop over cells
        do ind=1,twotondim
           
           ! Compute cell centre position in code units
           do idim=1,ndim
              nstride=2**(idim-1)
              xx(idim)=(2*m%grid(igrid)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5d0)*dx
           end do
           
           ! Check if cell is to be considered
           ok=(.NOT.m%grid(igrid)%refined(ind)).or.(ilevel==nlevelmax_frame)
           
           if(ok)then
              ! Check if the cell intersect the domain
#if NDIM>2                 
              if(r%proj_axis(ind_proj:ind_proj).eq.'x')then
                 xleft =xx(2)-dx/2.0d0
                 xright=xx(2)+dx/2.0d0
                 yleft =xx(3)-dx/2.0d0
                 yright=xx(3)+dx/2.0d0
              elseif(r%proj_axis(ind_proj:ind_proj).eq.'y')then
                 xleft =xx(1)-dx/2.0d0
                 xright=xx(1)+dx/2.0d0
                 yleft =xx(3)-dx/2.0d0
                 yright=xx(3)+dx/2.0d0
              else
                 xleft =xx(1)-dx/2.0d0
                 xright=xx(1)+dx/2.0d0
                 yleft =xx(2)-dx/2.0d0
                 yright=xx(2)+dx/2.0d0
              endif
              
              if(r%proj_axis(ind_proj:ind_proj).eq.'x')then
                 zleft =xx(1)-dx/2.
                 zright=xx(1)+dx/2.
              elseif(r%proj_axis(ind_proj:ind_proj).eq.'y')then
                 zleft =xx(2)-dx/2.
                 zright=xx(2)+dx/2.
              else
                 zleft =xx(3)-dx/2.
                 zright=xx(3)+dx/2.
              endif
              if(    xright.lt.xleft_frame.or.xleft.ge.xright_frame.or.&
                   & yright.lt.yleft_frame.or.yleft.ge.yright_frame.or.&
                   & zright.lt.zleft_frame.or.zleft.ge.zright_frame)cycle
#else
              xleft =xx(1)-dx/2.
              xright=xx(1)+dx/2.
#if NDIM>1
              yleft =xx(2)-dx/2.
              yright=xx(2)+dx/2.
#endif              
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
                    ind_map=ii+(jj-1)*r%nw_frame
#ifdef HYDRO
                    if(ind_var==0)then
                       ! Compute column density map
                       map(ind_map)=map(ind_map)+dvol*max(dble(m%grid(igrid)%uold(ind,1)),r%smallr)
                    else if(ind_var==1)then
                       ! Compute mass-weighted mean density
                       map(ind_map)=map(ind_map)+dvol*max(dble(m%grid(igrid)%uold(ind,1)),r%smallr)**2
                    else if(ind_var==5)then
                       ! Compute mass-weighted mean temperature
                       ! Kinetic energy
                       ekk=0.0d0
                       do idim=1,3
                          ekk=ekk+0.5d0*m%grid(igrid)%uold(ind,idim+1)**2/max(dble(m%grid(igrid)%uold(ind,1)),r%smallr)
                       enddo
                       ! Pressure
                       temp=(r%gamma-1.0d0)*(m%grid(igrid)%uold(ind,5)-ekk)
                       ! Temperature in K
                       temp=max(temp/max(dble(m%grid(igrid)%uold(ind,1)),r%smallr),r%smallc**2)*scale_T2
                       map(ind_map)=map(ind_map)+dvol*max(dble(m%grid(igrid)%uold(ind,1)),r%smallr)*temp
#ifdef RT
                    else if(ind_var.ge.nvar+1 .and. ind_var.lt.nvar+1+nrtgrp)then
                       ! Photon flux
                       map(ind_map) = map(ind_map) &
                                    + dvol * max(dble(m%grid(igrid)%uold(ind,1)),r%smallr) &
                                    * m%grid(igrid)%rtuold(ind,1+(ind_var-nvar-1)*(ndim+1)) * g%rt_c(ilevel)
#endif
                    else
                       ! Other variables
                       map(ind_map)=map(ind_map)+dvol*m%grid(igrid)%uold(ind,ind_var)
                    end if
#endif
                 end do
              end do
           end if
           
        end do
        ! End loop over cells
        
     end do
     ! End loop over grids
  end do
  ! End loop over levels
  
end subroutine output_frame
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
subroutine set_movie_vars(r)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtgrp
  ! This routine sets the movie vars from textual form
  type(run_t)::r
  integer::tmp
  character(LEN=5)::dummy
  integer::ll

  if(ANY(r%movie_vars_txt=='dens '))r%movie_vars(1)=1
  if(ANY(r%movie_vars_txt=='vx   '))r%movie_vars(2)=1
  if(ANY(r%movie_vars_txt=='vy   '))r%movie_vars(3)=1
  if(ANY(r%movie_vars_txt=='vz   '))r%movie_vars(4)=1
  if(ANY(r%movie_vars_txt=='temp '))r%movie_vars(5)=1
#if NVAR>5
  do ll=6,nvar
     write(dummy,'(I3.1)') ll
     if(ANY(r%movie_vars_txt=='var'//trim(adjustl(dummy))//' '))r%movie_vars(ll)=1
  end do
#endif
  do ll=nvar+1,nvar+nrtgrp
     tmp=ll-nvar
     write(dummy,'(I3.1)') tmp
     if(ANY(r%movie_vars_txt=='Fp'//trim(adjustl(dummy))//' '))r%movie_vars(ll)=1
  end do
  if(ANY(r%movie_vars_txt=='dm   '))r%movie_vars(NVAR+nrtgrp+1)=1
  if(ANY(r%movie_vars_txt=='stars'))r%movie_vars(NVAR+nrtgrp+2)=1
  
end subroutine set_movie_vars
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
end module movie_module
