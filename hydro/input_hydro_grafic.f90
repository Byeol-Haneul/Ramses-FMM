!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_input_hydro_grafic(pst,cpu_range,input_size,output_size,ilevel)
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::cpu_range,input_size,output_size
  integer::ilevel

  integer::next_range,next_cpu

  next_range=cpu_range/2
  next_cpu=pst%s%g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(pst%s%mdl,MDL_INPUT_HYDRO_GRAFIC,next_cpu,next_range,input_size,output_size,ilevel)
     call r_input_hydro_grafic(pst,next_range,input_size,output_size,ilevel)
  else
     call input_hydro_grafic(pst%s%r,pst%s%g,pst%s%m,ilevel)
  endif

end subroutine r_input_hydro_grafic
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_hydro_grafic(r,g,m,ilevel)
  use amr_parameters, only: ndim,twotondim,dp,nvector
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !--------------------------------------
  ! Compute initial conditions from files
  ! with the grafic format.
  !--------------------------------------
  integer::igrid,ilun
  integer::ind,idim,ivar
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,rr,vx,vy=0,vz=0,ek,ei,pp,xx1,xx2,xx3,dx_loc

  real(dp),allocatable,dimension(:,:,:)::init_array
  real(kind=4),allocatable,dimension(:,:)  ::init_plane

  logical::error,ok_file3
  character(LEN=80)::filename
  character(LEN=5)::nchar,ncharvar

  if(m%noct(ilevel)==0)return

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

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
  error=.false.
  if(i1_min<1.or.i1_max>g%n1(ilevel))error=.true.
  if(i2_min<1.or.i2_max>g%n2(ilevel))error=.true.
  if(i3_min<1.or.i3_max>g%n3(ilevel))error=.true.
  if(error) then
     write(*,*)'Some grid are outside initial conditions sub-volume'
     write(*,*)'for ilevel=',ilevel
     write(*,*)'and processor=',g%myid
     write(*,*)i1_min,i1_max
     write(*,*)i2_min,i2_max
     write(*,*)i3_min,i3_max
     write(*,*)g%n1(ilevel),g%n2(ilevel),g%n3(ilevel)
     call mdl_abort
  end if
  
  !------------------------------------------
  ! Second step: read initial condition files
  !------------------------------------------
  ! Allocate initial conditions array
  allocate(init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
  allocate(init_plane(1:g%n1(ilevel),1:g%n2(ilevel)))

  ! Loop over input variables
  do ivar=1,nvar
     if(r%cosmo)then
        ! Read baryons initial overdensity and displacement at a=aexp
        if(ivar==1)filename=TRIM(r%initfile(ilevel))//'/ic_deltab'
        if(ivar==2)filename=TRIM(r%initfile(ilevel))//'/ic_velcx'
        if(ivar==3)filename=TRIM(r%initfile(ilevel))//'/ic_velcy'
        if(ivar==4)filename=TRIM(r%initfile(ilevel))//'/ic_velcz'
        if(ivar==5)filename=TRIM(r%initfile(ilevel))//'/ic_tempb'
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
        open(10,file=filename,form='unformatted')
        rewind 10
        read(10) ! skip first line
        do i3=1,i3_min-1
           read(10)
        end do
        do i3=i3_min,i3_max
           read(10) ((init_plane(i1,i2),i1=1,g%n1(ilevel)),i2=1,g%n2(ilevel))
           init_array(i1_min:i1_max,i2_min:i2_max,i3) = init_plane(i1_min:i1_max,i2_min:i2_max)
        end do
        close(10)
     else
        ! If file doesn't exist, initialize variable to default value 
        ! In most cases, this is zero (you can change that if necessary)
        init_array=0d0
     endif
     
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
     
  end do
  ! End loop over input variables
  
  ! Deallocate initial conditions array
  deallocate(init_array)
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

end subroutine input_hydro_grafic
!################################################################
!################################################################
!################################################################
!################################################################
