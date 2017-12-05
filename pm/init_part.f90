!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_init_part(r,g,m,p,mdl,cpu_range,input_size,output_size)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size
  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to allocate
  ! particle-based arrays.
  !--------------------------------------------------------------------
  integer::next_range,next_cpu
  
  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_INIT_PART,next_cpu,next_range,input_size,output_size)
     call r_init_part(r,g,m,p,mdl,next_range,input_size,output_size)
  else
     call init_part(r,g,p)
  endif

end subroutine r_init_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_part(r,g,p)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  !------------------------------
  ! Allocate particle variables
  !------------------------------
  allocate(p%xp    (r%npartmax,ndim))
  allocate(p%vp    (r%npartmax,ndim))
  allocate(p%mp    (r%npartmax))
  allocate(p%levelp(r%npartmax))
  allocate(p%idp   (r%npartmax))
  allocate(p%sortp (r%npartmax))
  allocate(p%workp (r%npartmax))
#ifdef OUTPUT_PARTICLE_POTENTIAL
  allocate(p%phip  (r%npartmax))
#endif
  ! Allocate pointers to particle levels
  allocate(p%headp(r%levelmin:r%nlevelmax))
  allocate(p%tailp(r%levelmin:r%nlevelmax))
  ! No particle just yet
  p%headp=1
  p%tailp=0
end subroutine init_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_part_grid(r,g,m,p)
  use amr_parameters, only: dp,ndim,twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  real(kind=8)::mp_min_all
  integer::dummy_io,info,info2
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  !------------------------------------------------------------
  ! Initialize particle position from grid cell positions.
  ! Read particles displacement and velocities from grafic files.
  !------------------------------------------------------------
  integer::ipart,ipart_old,ilevel,idim
  integer::igrid
  integer::ind,ilun
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count
  integer,parameter::tagg=1109,tagg2=1110,tagg3=1111
  real(dp)::dx,dx_loc,xx1,xx2,xx3
  real(kind=8)::dispmax=0.0
  real(kind=4),allocatable,dimension(:,:)::init_plane,init_plane_x
  real(dp),allocatable,dimension(:,:,:)::init_array,init_array_x
  character(LEN=80)::filename,filename_x
  character(LEN=5)::nchar
  logical::ok,error,keep_part,read_pos=.false.

  if(r%verbose)write(*,*)'Entering init_part_from_grid'
  if(TRIM(r%filetype).NE.'grafic')return
  if(r%nrestart>0)return

  !----------------------------------------------------
  ! Reading initial conditions GRAFIC2 multigrid arrays
  !----------------------------------------------------
  ipart=0
  ! Loop over initial condition levels
  do ilevel=r%levelmin,r%nlevelmax
     
     if(r%initfile(ilevel)==' ')cycle
     
     ! Mesh size at level ilevel in normalised units
     dx=1.0/2**ilevel

     ! Mesh size at level ilevel in code units
     dx_loc=r%boxlen*dx
     
     !--------------------------------------------------------------
     ! First step: compute level boundaries and particle positions
     !--------------------------------------------------------------
     i1_min=g%n1(ilevel)+1; i1_max=0
     i2_min=g%n2(ilevel)+1; i2_max=0
     i3_min=g%n3(ilevel)+1; i3_max=0
     ipart_old=ipart

     ! Loop over grids
     do igrid=m%head(ilevel),m%tail(ilevel)
        do ind=1,twotondim
           ! Coordinates in normalised units (between 0 and 1)
           xx1=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx
           xx2=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx
           xx3=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx
           ! Compute integer coordinates in the frame of the file
           i1=int((xx1*(g%dxini(ilevel)/dx)-g%xoff1(ilevel))/g%dxini(ilevel))+1
           i2=int((xx2*(g%dxini(ilevel)/dx)-g%xoff2(ilevel))/g%dxini(ilevel))+1
           i3=int((xx3*(g%dxini(ilevel)/dx)-g%xoff3(ilevel))/g%dxini(ilevel))+1
           ! Compute min and max
           i1_min=MIN(i1_min,i1); i1_max=MAX(i1_max,i1)
           i2_min=MIN(i2_min,i2); i2_max=MAX(i2_max,i2)
           i3_min=MIN(i3_min,i3); i3_max=MAX(i3_max,i3)
           ! Create particle
           keep_part=.NOT.m%grid(igrid)%refined(ind)
           if(keep_part)then
              ipart=ipart+1
              if(ipart>r%npartmax)then
                 write(*,*)'Maximum number of particles incorrect'
                 write(*,*)'npartmax should be greater than',ipart
                 call clean_stop(g)
              endif
              if(ndim>0)p%xp(ipart,1)=xx1
              if(ndim>1)p%xp(ipart,2)=xx2
              if(ndim>2)p%xp(ipart,3)=xx3
              p%mp(ipart)=0.5d0**(3*ilevel)*(1.0d0-g%omega_b/g%omega_m)
              p%idp(ipart)=1+(i1-1)+(i2-1)*g%n1(ilevel)+(i3-1)*g%n1(ilevel)*g%n2(ilevel)
           endif
        end do
     end do

     ! Check that all grids are within initial condition region
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
           call clean_stop(g)
        end if
     endif
     
     !---------------------------------------------------------------------
     ! Second step: read initial condition file and set particle velocities
     !---------------------------------------------------------------------
     ! Allocate initial conditions array
     if(m%noct(ilevel)>0)then
        allocate(init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
        allocate(init_array_x(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
        init_array=0d0
        init_array_x=0d0
     end if
     allocate(init_plane(1:g%n1(ilevel),1:g%n2(ilevel)))
     allocate(init_plane_x(1:g%n1(ilevel),1:g%n2(ilevel)))

     ! Loop over input variables
     do idim=1,ndim
        
        ! Read dark matter initial displacement field
        if(idim==1)filename=TRIM(r%initfile(ilevel))//'/ic_velcx'
        if(idim==2)filename=TRIM(r%initfile(ilevel))//'/ic_velcy'
        if(idim==3)filename=TRIM(r%initfile(ilevel))//'/ic_velcz'
        
        if(idim==1)filename_x=TRIM(r%initfile(ilevel))//'/ic_poscx'
        if(idim==2)filename_x=TRIM(r%initfile(ilevel))//'/ic_poscy'
        if(idim==3)filename_x=TRIM(r%initfile(ilevel))//'/ic_poscz'
        
        INQUIRE(file=filename_x,exist=ok)
        if(.not.ok)then
           read_pos = .false.
        else
           read_pos = .true.
           if(g%myid==1)write(*,*)'Reading file '//TRIM(filename_x)
        end if
        
        if(g%myid==1)write(*,*)'Reading file '//TRIM(filename)

        if(g%myid==1)then
           open(10,file=filename,form='unformatted')
           rewind 10
           read(10) ! skip first line
        end if
        do i3=1,g%n3(ilevel)
           if(g%myid==1)then
              if(r%debug.and.mod(i3,10)==0)write(*,*)'Reading plane ',i3
              read(10)((init_plane(i1,i2),i1=1,g%n1(ilevel)),i2=1,g%n2(ilevel))
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

        if(read_pos) then
           if(g%myid==1)then
              open(10,file=filename_x,form='unformatted')
              rewind 10
              read(10) ! skip first line
           end if
           do i3=1,g%n3(ilevel)
              if(g%myid==1)then
                 if(r%debug.and.mod(i3,10)==0)write(*,*)'Reading plane ',i3
                 read(10)((init_plane_x(i1,i2),i1=1,g%n1(ilevel)),i2=1,g%n2(ilevel))
              else
                 init_plane_x=0.0
              endif
              buf_count=g%n1(ilevel)*g%n2(ilevel)
#ifndef WITHOUTMPI
              call MPI_BCAST(init_plane_x,buf_count,MPI_REAL,0,MPI_COMM_WORLD,info)
#endif
              if(m%noct(ilevel)>0)then
                 if(i3.ge.i3_min.and.i3.le.i3_max)then
                    init_array_x(i1_min:i1_max,i2_min:i2_max,i3) = &
                         & init_plane_x(i1_min:i1_max,i2_min:i2_max)
                 end if
              endif
           end do
           if(g%myid==1)close(10)
        end if

        if(m%noct(ilevel)>0)then
           ! Rescale initial displacement field to code units
           init_array=g%dfact(ilevel)*dx/g%dxini(ilevel)*init_array/g%vfact(ilevel)
           if(read_pos)then
              init_array_x = init_array_x/g%boxlen_ini
           endif
           ipart=ipart_old
           ! Loop over grids
           do igrid=m%head(ilevel),m%tail(ilevel)
              ! Loop over cells
              do ind=1,twotondim
                 ! Coordinates in normalised units (between 0 and 1)
                 xx1=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx
                 xx2=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx
                 xx3=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx
                 ! Compute integer coordinates in the frame of the file
                 i1=int((xx1*(g%dxini(ilevel)/dx)-g%xoff1(ilevel))/g%dxini(ilevel))+1
                 i2=int((xx2*(g%dxini(ilevel)/dx)-g%xoff2(ilevel))/g%dxini(ilevel))+1
                 i3=int((xx3*(g%dxini(ilevel)/dx)-g%xoff3(ilevel))/g%dxini(ilevel))+1
                 ! Compute particle velocity and identity
                 keep_part=.NOT.m%grid(igrid)%refined(ind)
                 if(keep_part)then
                    ipart=ipart+1
                    p%vp(ipart,idim)=init_array(i1,i2,i3)
                    if(.not. read_pos)then
                       dispmax=max(dispmax,abs(init_array(i1,i2,i3)/dx))
                    else
                       if(idim==1)p%xp(ipart,idim)=xx1+init_array_x(i1,i2,i3)
                       if(idim==2)p%xp(ipart,idim)=xx2+init_array_x(i1,i2,i3)
                       if(idim==3)p%xp(ipart,idim)=xx3+init_array_x(i1,i2,i3)
                       dispmax=max(dispmax,abs(init_array_x(i1,i2,i3)/dx))
                    endif
                 end if
              end do
              ! End loop over cells
           end do
           ! End loop over grids
        endif

     end do
     ! End loop over input variables

     ! Deallocate initial conditions array
     if(m%noct(ilevel)>0)then
        deallocate(init_array,init_array_x)
     end if
     deallocate(init_plane,init_plane_x)
     
     if(r%debug)write(*,*)'npart=',ipart,'/',r%npartmax,' for PE=',g%myid

  end do
  ! End loop over levels

  ! Initial particle number
  p%npart=ipart

  ! Move particle according to Zeldovich approximation
  if(.not. read_pos)then
     do ipart=1,p%npart
        p%xp(ipart,1:ndim)=p%xp(ipart,1:ndim)+p%vp(ipart,1:ndim)
     enddo
  endif

  ! Scale displacement to velocity
  do ipart=1,p%npart
     p%vp(ipart,1:ndim)=g%vfact(1)*p%vp(ipart,1:ndim)
  end do

  ! Periodic box
  do ipart=1,p%npart
#if NDIM>0
     if(p%xp(ipart,1)<   0.0d0 )p%xp(ipart,1)=p%xp(ipart,1)+r%boxlen
     if(p%xp(ipart,1)>=r%boxlen)p%xp(ipart,1)=p%xp(ipart,1)-r%boxlen
#endif
#if NDIM>1
     if(p%xp(ipart,2)<   0.0d0 )p%xp(ipart,2)=p%xp(ipart,2)+r%boxlen
     if(p%xp(ipart,2)>=r%boxlen)p%xp(ipart,2)=p%xp(ipart,2)-r%boxlen
#endif
#if NDIM>2
     if(p%xp(ipart,3)<   0.0d0 )p%xp(ipart,3)=p%xp(ipart,3)+r%boxlen
     if(p%xp(ipart,3)>=r%boxlen)p%xp(ipart,3)=p%xp(ipart,3)-r%boxlen
#endif
  end do

  ! Compute particle initial level
  do ipart=1,p%npart
     p%levelp(ipart)=r%levelmin
  end do

  ! Put all particles in levelmin
  p%headp=p%npart+1
  p%tailp=p%npart
  p%headp(r%levelmin)=1
  p%tailp(r%levelmin)=p%npart

end subroutine init_part_grid
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

