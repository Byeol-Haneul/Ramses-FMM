!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flow  
  use amr_commons
  implicit none

  integer::ilevel,ivar
  
  if(verbose)write(*,*)'Entering init_flow'
  do ilevel=nlevelmax,levelmin,-1
     call init_flow_fine(ilevel)
     call upload_fine(ilevel)
  end do
  if(verbose)write(*,*)'Complete init_flow'

end subroutine init_flow
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flow_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  
  integer::igrid,ngrid,ind,idim,nstride,i,ivar
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp)::dx
  logical::ok_file1,ok_file2,ok_file
  character(LEN=80)::filename

  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  !--------------------------------------
  ! Compute initial conditions from files
  !--------------------------------------
  filename=TRIM(initfile(ilevel))//'/ic_d'
  INQUIRE(file=filename,exist=ok_file1)
  if(multiple)then
     filename=TRIM(initfile(ilevel))//'/dir_deltab/ic_deltab.00001'
     INQUIRE(file=filename,exist=ok_file2)
  else
     filename=TRIM(initfile(ilevel))//'/ic_deltab'
     INQUIRE(file=filename,exist=ok_file2)
  endif
  ok_file=ok_file1.or.ok_file2
  if(ok_file)then
     call init_grafic(ilevel)
  else
  !----------------------------------------------------
  ! Compute initial conditions from subroutine condinit
  !----------------------------------------------------
  ! Mesh size at level ilevel in code units
  dx=boxlen/2**ilevel
  ! Loop over grids by vector sweeps
  do igrid=head(ilevel),tail(ilevel),nvector
     ngrid=MIN(nvector,tail(ilevel)-igrid+1)
     ! Loop over cells
     do ind=1,twotondim
        ! Compute cell centre position in code units
        do idim=1,ndim
           nstride=2**(idim-1)
           do i=1,ngrid
              xx(i,idim)=(2*grid(igrid+i-1)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx
           end do
        end do
        ! Call initial condition routine
        call condinit(xx,uu,dx,ngrid)
#ifdef HYDRO
        ! Scatter variables to main memory
        do ivar=1,nvar
           do i=1,ngrid
              grid(igrid+i-1)%uold(ind,ivar)=uu(i,ivar)
           end do
        end do
#endif
     end do
     ! End loop over cells
  end do
  ! End loop over grids
  endif

111 format('   Entering init_flow_fine for level ',I2)

end subroutine init_flow_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flow_fine_2(r,g,m,ilevel)
  use amr_parameters, only: ndim,twotondim,dp,nvector
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  
  ! Local variables
  integer::igrid,ngrid,ind,idim,nstride,i,ivar
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp)::dx
  logical::ok_file1,ok_file2,ok_file
  character(LEN=80)::filename

  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,111)ilevel

  !--------------------------------------
  ! Compute initial conditions from files
  !--------------------------------------
  filename=TRIM(r%initfile(ilevel))//'/ic_d'
  INQUIRE(file=filename,exist=ok_file1)
  if(r%multiple)then
     filename=TRIM(r%initfile(ilevel))//'/dir_deltab/ic_deltab.00001'
     INQUIRE(file=filename,exist=ok_file2)
  else
     filename=TRIM(r%initfile(ilevel))//'/ic_deltab'
     INQUIRE(file=filename,exist=ok_file2)
  endif
  ok_file=ok_file1.or.ok_file2
  if(ok_file)then
     call init_grafic_2(r,g,m,ilevel)
  else
  !----------------------------------------------------
  ! Compute initial conditions from subroutine condinit
  !----------------------------------------------------
  ! Mesh size at level ilevel in code units
  dx=r%boxlen/2**ilevel
  ! Loop over grids by vector sweeps
  do igrid=m%head(ilevel),m%tail(ilevel),nvector
     ngrid=MIN(nvector,m%tail(ilevel)-igrid+1)
     ! Loop over cells
     do ind=1,twotondim
        ! Compute cell centre position in code units
        do idim=1,ndim
           nstride=2**(idim-1)
           do i=1,ngrid
              xx(i,idim)=(2*m%grid(igrid+i-1)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx
           end do
        end do
        ! Call initial condition routine
        call condinit(xx,uu,dx,ngrid)
#ifdef HYDRO
        ! Scatter variables to main memory
        do ivar=1,nvar
           do i=1,ngrid
              m%grid(igrid+i-1)%uold(ind,ivar)=uu(i,ivar)
           end do
        end do
#endif
     end do
     ! End loop over cells
  end do
  ! End loop over grids
  endif

111 format('   Entering init_flow_fine for level ',I2)

end subroutine init_flow_fine_2
!################################################################
!################################################################
!################################################################
!################################################################
subroutine region_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn
  real(dp)::dx
  real(dp),dimension(1:nvector,1:nvar)::q
  real(dp),dimension(1:nvector,1:ndim)::x

  integer::i,ivar,k
  real(dp)::vol,r,xn,yn,zn,en

  ! Set some (tiny) default values in case n_region=0
  q(1:nn,1)=smallr
  q(1:nn,2)=0.0d0
#if NDIM>1
  q(1:nn,3)=0.0d0
#endif
#if NDIM>2
  q(1:nn,4)=0.0d0
#endif
  q(1:nn,ndim+2)=smallr*smallc**2/gamma
#if NVAR > NDIM + 2
  do ivar=ndim+3,nvar
     q(1:nn,ivar)=0.0d0
  end do
#endif

  ! Loop over initial conditions regions
  do k=1,nregion
     
     ! For "square" regions only:
     if(region_type(k) .eq. 'square')then
        ! Exponent of choosen norm
        en=exp_region(k)
        do i=1,nn
           ! Compute position in normalized coordinates
           xn=0.0d0; yn=0.0d0; zn=0.0d0
           xn=2.0d0*abs(x(i,1)-x_center(k))/length_x(k)
#if NDIM>1
           yn=2.0d0*abs(x(i,2)-y_center(k))/length_y(k)
#endif
#if NDIM>2
           zn=2.0d0*abs(x(i,3)-z_center(k))/length_z(k)
#endif
           ! Compute cell "radius" relative to region center
           if(exp_region(k)<10)then
              r=(xn**en+yn**en+zn**en)**(1.0/en)
           else
              r=max(xn,yn,zn)
           end if
           ! If cell lies within region,
           ! REPLACE primitive variables by region values
           if(r<1.0)then
              q(i,1)=d_region(k)
              q(i,2)=u_region(k)
#if NDIM>1
              q(i,3)=v_region(k)
#endif
#if NDIM>2
              q(i,4)=w_region(k)
#endif
              q(i,ndim+2)=p_region(k)
#if NENER>0
              do ivar=1,nener
                 q(i,ndim+2+ivar)=prad_region(k,ivar)
              enddo
#endif
#if NVAR>NDIM+2+NENER
              do ivar=ndim+3+nener,nvar
                 q(i,ivar)=var_region(k,ivar-ndim-2-nener)
              end do
#endif
           end if
        end do
     end if
     
     ! For "point" regions only:
     if(region_type(k) .eq. 'point')then
        ! Volume elements
        vol=dx**ndim
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=1.0; yn=1.0; zn=1.0
           xn=max(1.0-abs(x(i,1)-x_center(k))/dx,0.0_dp)
#if NDIM>1
           yn=max(1.0-abs(x(i,2)-y_center(k))/dx,0.0_dp)
#endif
#if NDIM>2
           zn=max(1.0-abs(x(i,3)-z_center(k))/dx,0.0_dp)
#endif
           r=xn*yn*zn
           ! If cell lies within CIC cloud, 
           ! ADD to primitive variables the region values
           q(i,1)=q(i,1)+d_region(k)*r/vol
           q(i,2)=q(i,2)+u_region(k)*r
#if NDIM>1
           q(i,3)=q(i,3)+v_region(k)*r
#endif
#if NDIM>2
           q(i,4)=q(i,4)+w_region(k)*r
#endif
           q(i,ndim+2)=q(i,ndim+2)+p_region(k)*r/vol
#if NENER>0
           do ivar=1,nener
              q(i,ndim+2+ivar)=q(i,ndim+2+ivar)+prad_region(k,ivar)*r/vol
           enddo
#endif
#if NVAR>NDIM+2+NENER
           do ivar=ndim+3+nener,nvar
              q(i,ivar)=var_region(k,ivar-ndim-2-nener)
           end do
#endif
        end do
     end if
  end do

  return
end subroutine region_condinit
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_grafic(ilevel)
  use amr_commons
  use hydro_commons
!  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !--------------------------------------
  ! Compute initial conditions from files
  ! with the grafic format.
  !--------------------------------------
  integer::i,igrid,ilun
  integer::ind,idim,ivar,ix,iy,iz,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count,info,nvar_in

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,rr,vx,vy=0,vz=0,ek,ei,pp,xx1,xx2,xx3,dx_loc,xval

  real(dp),allocatable,dimension(:,:,:)::init_array
  real(kind=4),allocatable,dimension(:,:)  ::init_plane

  logical::error,ok_file3
  character(LEN=80)::filename
  character(LEN=5)::nchar,ncharvar

  integer,parameter::tag=1107
  integer::dummy_io,info2

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh size at level ilevel in normalised units
  dx=0.5D0**ilevel
    
  ! Mesh size at level ilevel in code units
  dx_loc=boxlen*dx

  !-------------------------------------------------------------------------
  ! First step: compute level boundaries in terms of initial condition array
  !-------------------------------------------------------------------------
  i1_min=n1(ilevel)+1; i1_max=0
  i2_min=n2(ilevel)+1; i2_max=0
  i3_min=n3(ilevel)+1; i3_max=0
  do igrid=head(ilevel),tail(ilevel)
     do ind=1,twotondim
        ! Coordinates in normalised units (between 0 and 1)
        xx1=(2*grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx
        xx2=(2*grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx
        xx3=(2*grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx
        ! Scale to integer coordinates in the frame of the file
        xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
        xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
        xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
        ! Compute min and max
        i1_min=MIN(i1_min,int(xx1)+1); i1_max=MAX(i1_max,int(xx1)+1)
        i2_min=MIN(i2_min,int(xx2)+1); i2_max=MAX(i2_max,int(xx2)+1)
        i3_min=MIN(i3_min,int(xx3)+1); i3_max=MAX(i3_max,int(xx3)+1)
     end do
  end do
  if(noct(ilevel)>0)then
     error=.false.
     if(i1_min<1.or.i1_max>n1(ilevel))error=.true.
     if(i2_min<1.or.i2_max>n2(ilevel))error=.true.
     if(i3_min<1.or.i3_max>n3(ilevel))error=.true.
     if(error) then
        write(*,*)'Some grid are outside initial conditions sub-volume'
        write(*,*)'for ilevel=',ilevel
        write(*,*)i1_min,i1_max
        write(*,*)i2_min,i2_max
        write(*,*)i3_min,i3_max
        write(*,*)n1(ilevel),n2(ilevel),n3(ilevel)
        call clean_stop
     end if
  endif

  !------------------------------------------
  ! Second step: read initial condition files
  !------------------------------------------
  ! Allocate initial conditions array
  if(noct(ilevel)>0)allocate(init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
  allocate(init_plane(1:n1(ilevel),1:n2(ilevel)))
  ! Loop over input variables
  do ivar=1,nvar
     if(cosmo)then
        ! Read baryons initial overdensity and displacement at a=aexp
        if(multiple)then
           call title(myid,nchar)
           if(ivar==1)filename=TRIM(initfile(ilevel))//'/dir_deltab/ic_deltab.'//TRIM(nchar)
           if(ivar==2)filename=TRIM(initfile(ilevel))//'/dir_velcx/ic_velcx.'//TRIM(nchar)
           if(ivar==3)filename=TRIM(initfile(ilevel))//'/dir_velcy/ic_velcy.'//TRIM(nchar)
           if(ivar==4)filename=TRIM(initfile(ilevel))//'/dir_velcz/ic_velcz.'//TRIM(nchar)
           if(ivar==5)filename=TRIM(initfile(ilevel))//'/dir_tempb/ic_tempb.'//TRIM(nchar)
        else
           if(ivar==1)filename=TRIM(initfile(ilevel))//'/ic_deltab'
           if(ivar==2)filename=TRIM(initfile(ilevel))//'/ic_velcx'
           if(ivar==3)filename=TRIM(initfile(ilevel))//'/ic_velcy'
           if(ivar==4)filename=TRIM(initfile(ilevel))//'/ic_velcz'
           if(ivar==5)filename=TRIM(initfile(ilevel))//'/ic_tempb'
        endif
     else
        ! Read primitive variables
        if(ivar==1)filename=TRIM(initfile(ilevel))//'/ic_d'
        if(ivar==2)filename=TRIM(initfile(ilevel))//'/ic_u'
        if(ivar==3)filename=TRIM(initfile(ilevel))//'/ic_v'
        if(ivar==4)filename=TRIM(initfile(ilevel))//'/ic_w'
        if(ivar==5)filename=TRIM(initfile(ilevel))//'/ic_p'
     endif
     call title(ivar,ncharvar)
     if(ivar>5)then
        call title(ivar-5,ncharvar)
        filename=TRIM(initfile(ilevel))//'/ic_pvar_'//TRIM(ncharvar)
     endif
     
     INQUIRE(file=filename,exist=ok_file3)
     if(ok_file3)then
        ! Reading the existing file   
        if(myid==1)write(*,*)'Reading file '//TRIM(filename)
        if(multiple)then
           ilun=ncpu+myid+10
           
           ! Wait for the token
#ifndef WITHOUTMPI
           if(IOGROUPSIZE>0) then
              if (mod(myid-1,IOGROUPSIZE)/=0) then
                 call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                      & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
              end if
           endif
#endif
           
           open(ilun,file=filename,form='unformatted')
           rewind ilun
           read(ilun) ! skip first line
           do i3=1,n3(ilevel)
              read(ilun) ((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
              if(noct(ilevel)>0)then
                 if(i3.ge.i3_min.and.i3.le.i3_max)then
                    init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                         & init_plane(i1_min:i1_max,i2_min:i2_max)
                 end if
              endif
           end do
           close(ilun)
           ! Send the token
#ifndef WITHOUTMPI
           if(IOGROUPSIZE>0) then
              if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
                 dummy_io=1
                 call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
                      & MPI_COMM_WORLD,info2)
              end if
           endif
#endif
        else
           if(myid==1)then
              open(10,file=filename,form='unformatted')
              rewind 10
              read(10) ! skip first line
           endif
           do i3=1,n3(ilevel)
              if(myid==1)then
                 read(10) ((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
              else
                 init_plane=0.0
              endif
              buf_count=n1(ilevel)*n2(ilevel)
#ifndef WITHOUTMPI
              call MPI_BCAST(init_plane,buf_count,MPI_REAL,0,MPI_COMM_WORLD,info)
#endif
              if(noct(ilevel)>0)then
                 if(i3.ge.i3_min.and.i3.le.i3_max)then
                    init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                         & init_plane(i1_min:i1_max,i2_min:i2_max)
                 end if
              endif
           end do
           if(myid==1)close(10)
        endif
     else
        ! If file doesn't exist, initialize variable to default value 
        ! In most cases, this is zero (you can change that if necessary)
        if(myid==1)write(*,*)'File '//TRIM(filename)//' not found'
        if(myid==1)write(*,*)'Initialize corresponding variable to default value'
        if(noct(ilevel)>0)then
           init_array=0d0
!!$           ! Default value for metals
!!$           if(cosmo.and.ivar==imetal.and.metal)init_array=z_ave*0.02 ! from solar units
!!$           ! Default value for ionization fraction
!!$           if(cosmo)xval=sqrt(omega_m)/(h0/100.*omega_b) ! From the book of Peebles p. 173
!!$           if(cosmo.and.ivar==ixion.and.aton)init_array=1.2d-5*xval
        endif
     endif
     
     if(noct(ilevel)>0)then
        
        ! For cosmo runs, rescale initial conditions to code units
        if(cosmo)then
           ! Compute approximate average temperature in K
           T2_start=1.356d-2/aexp**2
           if(ivar==1)init_array=(1.0+dfact(ilevel)*init_array)*omega_b/omega_m
           if(ivar==2)init_array=dfact(ilevel)*vfact(1)*dx_loc/dxini(ilevel)*init_array/vfact(ilevel)
           if(ivar==3)init_array=dfact(ilevel)*vfact(1)*dx_loc/dxini(ilevel)*init_array/vfact(ilevel)
           if(ivar==4)init_array=dfact(ilevel)*vfact(1)*dx_loc/dxini(ilevel)*init_array/vfact(ilevel)
           if(ivar==ndim+2)init_array=(1.0+init_array)*T2_start/scale_T2
        endif
        
        ! Loop over cells
        do igrid=head(ilevel),tail(ilevel)
           do ind=1,twotondim
              ! Coordinates in normalised units (between 0 and 1)
              xx1=(2*grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx
              xx2=(2*grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx
              xx3=(2*grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx
              ! Scale to integer coordinates in the frame of the file
              xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
              xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
              xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
              ! Compute integer coordinates in the frame of the file
              i1=int(xx1)+1
              i2=int(xx2)+1
              i3=int(xx3)+1
#ifdef HYDRO
              ! Scatter to corresponding primitive variable
              grid(igrid)%uold(ind,ivar)=init_array(i1,i2,i3)
#endif
           end do
        end do
        ! End loop over cells

     endif
  end do
  ! End loop over input variables

  ! Deallocate initial conditions array
  if(noct(ilevel)>0)deallocate(init_array)
  deallocate(init_plane) 

  !----------------------------------------------------------------
  ! For cosmology runs: compute pressure, prevent negative density
  !----------------------------------------------------------------
  if(cosmo)then
     ! Loop over grids
     do igrid=head(ilevel),tail(ilevel)
        ! Loop over cells
        do ind=1,twotondim
#ifdef HYDRO
           ! Prevent negative densities
           rr=max(grid(igrid)%uold(ind,1),0.1*omega_b/omega_m)
           grid(igrid)%uold(ind,1)=rr
           ! Compute pressure from temperature and density
           grid(igrid)%uold(ind,ndim+2)=grid(igrid)%uold(ind,1)*grid(igrid)%uold(ind,ndim+2)
#endif
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end if

  !---------------------------------------------------
  ! Third step: compute initial conservative variables
  !---------------------------------------------------
  ! Loop over grids
  do igrid=head(ilevel),tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
#ifdef HYDRO
        ! Compute total energy density
        rr=grid(igrid)%uold(ind,1)
        vx=grid(igrid)%uold(ind,2)
#if NDIM>1
        vy=grid(igrid)%uold(ind,3)
#endif
#if NDIM>2
        vz=grid(igrid)%uold(ind,4)
#endif
        pp=grid(igrid)%uold(ind,ndim+2)
        ek=0.5d0*rr*(vx**2+vy**2+vz**2)
        ei=pp/(gamma-1.0)
        grid(igrid)%uold(ind,ndim+2)=ei+ek
        ! Compute momentum density
        do idim=1,ndim
           rr=grid(igrid)%uold(ind,1)
           grid(igrid)%uold(ind,idim+1)=rr*grid(igrid)%uold(ind,idim+1)
        end do
#if NVAR>NDIM+2
        ! Compute passive variable density
        do ivar=ndim+3,nvar
           rr=grid(igrid)%uold(ind,1)
           grid(igrid)%uold(ind,ivar)=rr*grid(igrid)%uold(ind,ivar)
        enddo
#endif
#endif
     end do
     ! End loop over cells
  end do
  ! End loop over grids

end subroutine init_grafic
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_grafic_2(r,g,m,ilevel)
  use amr_parameters, only: ndim,twotondim,dp,nvector
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !--------------------------------------
  ! Compute initial conditions from files
  ! with the grafic format.
  !--------------------------------------
  integer::i,igrid,ilun
  integer::ind,idim,ivar,ix,iy,iz,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count,info,nvar_in

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,rr,vx,vy=0,vz=0,ek,ei,pp,xx1,xx2,xx3,dx_loc,xval

  real(dp),allocatable,dimension(:,:,:)::init_array
  real(kind=4),allocatable,dimension(:,:)  ::init_plane

  logical::error,ok_file3
  character(LEN=80)::filename
  character(LEN=5)::nchar,ncharvar

  integer,parameter::tag=1107
  integer::dummy_io,info2

  ! Conversion factor from user units to cgs units
  call units_2(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh size at level ilevel in normalised units
  dx=0.5D0**ilevel
    
  ! Mesh size at level ilevel in code units
  dx_loc=r%boxlen*dx

  !-------------------------------------------------------------------------
  ! First step: compute level boundaries in terms of initial condition array
  !-------------------------------------------------------------------------
  i1_min=g%n1(ilevel)+1; i1_max=0
  i2_min=g%n2(ilevel)+1; i2_max=0
  i3_min=g%n3(ilevel)+1; i3_max=0
  do igrid=m%head(ilevel),m%tail(ilevel)
     do ind=1,twotondim
        ! Coordinates in normalised units (between 0 and 1)
        xx1=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx
        xx2=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx
        xx3=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx
        ! Scale to integer coordinates in the frame of the file
        xx1=(xx1*(g%dxini(ilevel)/dx)-g%xoff1(ilevel))/g%dxini(ilevel)
        xx2=(xx2*(g%dxini(ilevel)/dx)-g%xoff2(ilevel))/g%dxini(ilevel)
        xx3=(xx3*(g%dxini(ilevel)/dx)-g%xoff3(ilevel))/g%dxini(ilevel)
        ! Compute min and max
        i1_min=MIN(i1_min,int(xx1)+1); i1_max=MAX(i1_max,int(xx1)+1)
        i2_min=MIN(i2_min,int(xx2)+1); i2_max=MAX(i2_max,int(xx2)+1)
        i3_min=MIN(i3_min,int(xx3)+1); i3_max=MAX(i3_max,int(xx3)+1)
     end do
  end do
  if(m%noct(ilevel)>0)then
     error=.false.
     if(i1_min<1.or.i1_max>g%n1(ilevel))error=.true.
     if(i2_min<1.or.i2_max>g%n2(ilevel))error=.true.
     if(i3_min<1.or.i3_max>g%n3(ilevel))error=.true.
     if(error) then
        write(*,*)'Some grid are outside initial conditions sub-volume'
        write(*,*)'for ilevel=',ilevel
        write(*,*)i1_min,i1_max
        write(*,*)i2_min,i2_max
        write(*,*)i3_min,i3_max
        write(*,*)g%n1(ilevel),g%n2(ilevel),g%n3(ilevel)
        call clean_stop
     end if
  endif

  !------------------------------------------
  ! Second step: read initial condition files
  !------------------------------------------
  ! Allocate initial conditions array
  if(m%noct(ilevel)>0)allocate(init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
  allocate(init_plane(1:g%n1(ilevel),1:g%n2(ilevel)))
  ! Loop over input variables
  do ivar=1,nvar
     if(r%cosmo)then
        ! Read baryons initial overdensity and displacement at a=aexp
        if(r%multiple)then
           call title(g%myid,nchar)
           if(ivar==1)filename=TRIM(r%initfile(ilevel))//'/dir_deltab/ic_deltab.'//TRIM(nchar)
           if(ivar==2)filename=TRIM(r%initfile(ilevel))//'/dir_velcx/ic_velcx.'//TRIM(nchar)
           if(ivar==3)filename=TRIM(r%initfile(ilevel))//'/dir_velcy/ic_velcy.'//TRIM(nchar)
           if(ivar==4)filename=TRIM(r%initfile(ilevel))//'/dir_velcz/ic_velcz.'//TRIM(nchar)
           if(ivar==5)filename=TRIM(r%initfile(ilevel))//'/dir_tempb/ic_tempb.'//TRIM(nchar)
        else
           if(ivar==1)filename=TRIM(r%initfile(ilevel))//'/ic_deltab'
           if(ivar==2)filename=TRIM(r%initfile(ilevel))//'/ic_velcx'
           if(ivar==3)filename=TRIM(r%initfile(ilevel))//'/ic_velcy'
           if(ivar==4)filename=TRIM(r%initfile(ilevel))//'/ic_velcz'
           if(ivar==5)filename=TRIM(r%initfile(ilevel))//'/ic_tempb'
        endif
     else
        ! Read primitive variables
        if(ivar==1)filename=TRIM(r%initfile(ilevel))//'/ic_d'
        if(ivar==2)filename=TRIM(r%initfile(ilevel))//'/ic_u'
        if(ivar==3)filename=TRIM(r%initfile(ilevel))//'/ic_v'
        if(ivar==4)filename=TRIM(r%initfile(ilevel))//'/ic_w'
        if(ivar==5)filename=TRIM(r%initfile(ilevel))//'/ic_p'
     endif
     call title(ivar,ncharvar)
     if(ivar>5)then
        call title(ivar-5,ncharvar)
        filename=TRIM(r%initfile(ilevel))//'/ic_pvar_'//TRIM(ncharvar)
     endif
     
     INQUIRE(file=filename,exist=ok_file3)
     if(ok_file3)then
        ! Reading the existing file   
        if(g%myid==1)write(*,*)'Reading file '//TRIM(filename)
        if(r%multiple)then
           ilun=g%ncpu+g%myid+10
           
           ! Wait for the token
#ifndef WITHOUTMPI
           if(g%IOGROUPSIZE>0) then
              if (mod(g%myid-1,g%IOGROUPSIZE)/=0) then
                 call MPI_RECV(dummy_io,1,MPI_INTEGER,g%myid-1-1,tag,&
                      & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
              end if
           endif
#endif
           
           open(ilun,file=filename,form='unformatted')
           rewind ilun
           read(ilun) ! skip first line
           do i3=1,g%n3(ilevel)
              read(ilun) ((init_plane(i1,i2),i1=1,g%n1(ilevel)),i2=1,g%n2(ilevel))
              if(m%noct(ilevel)>0)then
                 if(i3.ge.i3_min.and.i3.le.i3_max)then
                    init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                         & init_plane(i1_min:i1_max,i2_min:i2_max)
                 end if
              endif
           end do
           close(ilun)
           ! Send the token
#ifndef WITHOUTMPI
           if(g%IOGROUPSIZE>0) then
              if(mod(g%myid,g%IOGROUPSIZE)/=0 .and.(g%myid.lt.g%ncpu))then
                 dummy_io=1
                 call MPI_SEND(dummy_io,1,MPI_INTEGER,g%myid-1+1,tag, &
                      & MPI_COMM_WORLD,info2)
              end if
           endif
#endif
        else
           if(g%myid==1)then
              open(10,file=filename,form='unformatted')
              rewind 10
              read(10) ! skip first line
           endif
           do i3=1,g%n3(ilevel)
              if(g%myid==1)then
                 read(10) ((init_plane(i1,i2),i1=1,g%n1(ilevel)),i2=1,g%n2(ilevel))
              else
                 init_plane=0.0
              endif
              buf_count=g%n1(ilevel)*g%n2(ilevel)
#ifndef WITHOUTMPI
              call MPI_BCAST(init_plane,buf_count,MPI_REAL,0,MPI_COMM_WORLD,info)
#endif
              if(m%noct(ilevel)>0)then
                 if(i3.ge.i3_min.and.i3.le.i3_max)then
                    init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                         & init_plane(i1_min:i1_max,i2_min:i2_max)
                 end if
              endif
           end do
           if(g%myid==1)close(10)
        endif
     else
        ! If file doesn't exist, initialize variable to default value 
        ! In most cases, this is zero (you can change that if necessary)
        if(g%myid==1)write(*,*)'File '//TRIM(filename)//' not found'
        if(g%myid==1)write(*,*)'Initialize corresponding variable to default value'
        if(m%noct(ilevel)>0)then
           init_array=0d0
!!$           ! Default value for metals
!!$           if(r%cosmo.and.ivar==g%imetal.and.r%metal)init_array=r%z_ave*0.02 ! from solar units
!!$           ! Default value for ionization fraction
!!$           if(r%cosmo)xval=sqrt(g%omega_m)/(g%h0/100.*g%omega_b) ! From the book of Peebles p. 173
!!$           if(r%cosmo.and.ivar==g%ixion.and.r%aton)init_array=1.2d-5*g%xval
        endif
     endif
     
     if(m%noct(ilevel)>0)then
        
        ! For cosmo runs, rescale initial conditions to code units
        if(r%cosmo)then
           ! Compute approximate average temperature in K
           g%T2_start=1.356d-2/g%aexp**2
           if(ivar==1)init_array=(1.0+g%dfact(ilevel)*init_array)*g%omega_b/g%omega_m
           if(ivar==2)init_array=g%dfact(ilevel)*g%vfact(1)*dx_loc/g%dxini(ilevel)*init_array/g%vfact(ilevel)
           if(ivar==3)init_array=g%dfact(ilevel)*g%vfact(1)*dx_loc/g%dxini(ilevel)*init_array/g%vfact(ilevel)
           if(ivar==4)init_array=g%dfact(ilevel)*g%vfact(1)*dx_loc/g%dxini(ilevel)*init_array/g%vfact(ilevel)
           if(ivar==ndim+2)init_array=(1.0+init_array)*g%T2_start/scale_T2
        endif
        
        ! Loop over cells
        do igrid=m%head(ilevel),m%tail(ilevel)
           do ind=1,twotondim
              ! Coordinates in normalised units (between 0 and 1)
              xx1=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx
              xx2=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx
              xx3=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx
              ! Scale to integer coordinates in the frame of the file
              xx1=(xx1*(g%dxini(ilevel)/dx)-g%xoff1(ilevel))/g%dxini(ilevel)
              xx2=(xx2*(g%dxini(ilevel)/dx)-g%xoff2(ilevel))/g%dxini(ilevel)
              xx3=(xx3*(g%dxini(ilevel)/dx)-g%xoff3(ilevel))/g%dxini(ilevel)
              ! Compute integer coordinates in the frame of the file
              i1=int(xx1)+1
              i2=int(xx2)+1
              i3=int(xx3)+1
#ifdef HYDRO
              ! Scatter to corresponding primitive variable
              m%grid(igrid)%uold(ind,ivar)=init_array(i1,i2,i3)
#endif
           end do
        end do
        ! End loop over cells

     endif
  end do
  ! End loop over input variables

  ! Deallocate initial conditions array
  if(m%noct(ilevel)>0)deallocate(init_array)
  deallocate(init_plane) 

  !----------------------------------------------------------------
  ! For cosmology runs: compute pressure, prevent negative density
  !----------------------------------------------------------------
  if(r%cosmo)then
     ! Loop over grids
     do igrid=m%head(ilevel),m%tail(ilevel)
        ! Loop over cells
        do ind=1,twotondim
#ifdef HYDRO
           ! Prevent negative densities
           rr=max(m%grid(igrid)%uold(ind,1),0.1*g%omega_b/g%omega_m)
           m%grid(igrid)%uold(ind,1)=rr
           ! Compute pressure from temperature and density
           m%grid(igrid)%uold(ind,ndim+2)=m%grid(igrid)%uold(ind,1)*m%grid(igrid)%uold(ind,ndim+2)
#endif
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end if

  !---------------------------------------------------
  ! Third step: compute initial conservative variables
  !---------------------------------------------------
  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
#ifdef HYDRO
        ! Compute total energy density
        rr=m%grid(igrid)%uold(ind,1)
        vx=m%grid(igrid)%uold(ind,2)
#if NDIM>1
        vy=m%grid(igrid)%uold(ind,3)
#endif
#if NDIM>2
        vz=m%grid(igrid)%uold(ind,4)
#endif
        pp=m%grid(igrid)%uold(ind,ndim+2)
        ek=0.5d0*rr*(vx**2+vy**2+vz**2)
        ei=pp/(r%gamma-1.0)
        m%grid(igrid)%uold(ind,ndim+2)=ei+ek
        ! Compute momentum density
        do idim=1,ndim
           rr=m%grid(igrid)%uold(ind,1)
           m%grid(igrid)%uold(ind,idim+1)=rr*m%grid(igrid)%uold(ind,idim+1)
        end do
#if NVAR>NDIM+2
        ! Compute passive variable density
        do ivar=ndim+3,nvar
           rr=m%grid(igrid)%uold(ind,1)
           m%grid(igrid)%uold(ind,ivar)=rr*m%grid(igrid)%uold(ind,ivar)
        enddo
#endif
#endif
     end do
     ! End loop over cells
  end do
  ! End loop over grids

end subroutine init_grafic_2
!################################################################
!################################################################
!################################################################
!################################################################
