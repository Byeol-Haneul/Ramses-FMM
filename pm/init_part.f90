!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_part_file_2(r,g,p)
  use amr_parameters, only: ndim,dp,i8b
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  real(dp)::mp_min_all
  integer::info
#endif
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  !------------------------------------------------------------
  ! Allocate particle-based arrays.
  ! Read particles positions and velocities from various files
  ! including gadget, ascii or restart files.
  ! grafic initial conditions are performed after the AMR grid 
  ! has been constructed.
  !------------------------------------------------------------
  integer::ipart,jpart,jpart_loc=0,ipart_old,idim
  integer::i,ilun,icpu
  integer::indglob
  real(dp)::xx1,xx2,xx3,vv1,vv2,vv3,mm1
  integer::ncpu_file
  integer::ileft,iright,nrest,ipos
  integer(i8b)::npart_tot_file
  integer(i8b)::istart,iend,nleft,nright

  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::isp
  integer(i8b),allocatable,dimension(:)::isp8
  integer,allocatable,dimension(:)::npart_file
  integer(i8b),allocatable,dimension(:)::ncum_file

  integer,dimension(1:g%ncpu+1)::start_ind
  character(LEN=80)::filename
  character(LEN=80)::fileloc,file_part
  character(LEN=20)::filetype_loc
  character(LEN=5)::nchar,ncharcpu

  if(r%verbose)write(*,*)'Entering init_part_from_file'

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
  
  if(r%nrestart>0)then

     !----------------------------------------
     ! In case of restart, read RAMSES file
     !----------------------------------------

     ! Read particle files header
     call title(r%nrestart,nchar)
     fileloc='output_'//TRIM(nchar)//'/part_header.txt'
     call input_header_2(r,g,fileloc,npart_tot_file,ncpu_file)
     if(g%myid==1)write(*,'(" Restart snapshot has ",I8," particles")')npart_tot_file

     ! Compute new local particle number
     p%npart_tot=npart_tot_file
     p%npart=npart_tot_file/g%ncpu
     nleft=(g%myid-1)*p%npart
     nright=g%myid*p%npart
     nrest=npart_tot_file-p%npart*g%ncpu
     if(g%myid.LE.nrest)then
        p%npart=p%npart+1
        nright=nright+g%myid
        nleft=nleft+(g%myid-1)
     endif

     ! Compute number of particles in each file
     allocate(npart_file(1:ncpu_file))
     allocate(ncum_file(0:ncpu_file))
     if(g%myid==1)then
        ncum_file(0)=0
        do icpu=1,ncpu_file
           call title(icpu,ncharcpu)
           file_part='output_'//TRIM(nchar)//'/part.out'//TRIM(ncharcpu)
           ilun=10
           open(unit=ilun,file=TRIM(file_part),access="stream"&
                & ,action="read",form='unformatted')
           read(ilun,POS=5)npart_file(icpu)
           ncum_file(icpu)=ncum_file(icpu-1)+npart_file(icpu)
           close(ilun)
        end do
     endif
#ifndef WITHOUTMPI
     call MPI_BCAST(npart_file,ncpu_file,MPI_INTEGER,0,MPI_COMM_WORLD,info)
     call MPI_BCAST(ncum_file,ncpu_file+1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
#endif

     ! Compute interval of file to open for current process
     ileft=0
     iright=-1
     if(nright.GT.nleft)then
        do icpu=1,ncpu_file
           if(icpu>1)then
              if(ncum_file(icpu).GT.nleft.AND.ncum_file(icpu-1).LT.nright)then
                 if(ileft==0)ileft=icpu
                 iright=MAX(icpu,iright)
              endif
           else
              if(ncum_file(icpu).GT.nleft)then
                 if(ileft==0)ileft=icpu
                 iright=MAX(icpu,iright)
              endif
           endif
        end do
     endif

     ! Loop over relevant files
     ipart=0
     ipart_old=0
     do icpu=ileft,iright
        if(icpu>1)then
           istart=MAX(nleft-ncum_file(icpu-1),0)+1
           iend=MIN(nright-ncum_file(icpu-1),npart_file(icpu))
        else
           istart=nleft+1
           iend=MIN(nright,npart_file(icpu))
        endif

        ! Open the PART file
        call title(icpu,ncharcpu)
        file_part='output_'//TRIM(nchar)//'/part.out'//TRIM(ncharcpu)
        open(unit=10,file=TRIM(file_part),access="stream"&
             & ,action="read",form='unformatted')

        ! Read positions
        allocate(xdp(istart:iend))
        do idim=1,ndim
           ipos=9+8*(idim-1)*npart_file(icpu)+8*(istart-1)
           read(10,POS=ipos)xdp
           ipart=ipart_old
           do i=istart,iend
              ipart=ipart+1
              p%xp(ipart,idim)=xdp(i)
           end do
        end do

        ! Read velocities
        do idim=1,ndim
           ipos=9+8*(ndim+idim-1)*npart_file(icpu)+8*(istart-1)
           read(10,POS=ipos)xdp
           ipart=ipart_old
           do i=istart,iend
              ipart=ipart+1
              p%vp(ipart,idim)=xdp(i)
           end do
        end do

        ! Read masses
        ipos=9+8*(ndim+ndim)*npart_file(icpu)+8*(istart-1)
        read(10,POS=ipos)xdp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%mp(ipart)=xdp(i)
        end do
        deallocate(xdp)

        ! Read identity
        allocate(isp8(istart:iend))
#ifndef LONGINT
        ipos=9+8*(ndim+ndim+1)*npart_file(icpu)+4*(istart-1)
#else
        ipos=9+8*(ndim+ndim+1)*npart_file(icpu)+8*(istart-1)
#endif
        read(10,POS=ipos)isp8
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%idp(ipart)=isp8(i)
        end do
        deallocate(isp8)

        ! Read level
        allocate(isp(istart:iend))
#ifndef LONGINT
        ipos=9+8*(ndim+ndim+1)*npart_file(icpu)+4*npart_file(icpu)+4*(istart-1)
#else
        ipos=9+8*(ndim+ndim+1)*npart_file(icpu)+8*npart_file(icpu)+4*(istart-1)
#endif
        read(10,POS=ipos)isp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%levelp(ipart)=isp(i)
        end do
        deallocate(isp)

        ! Close the PART file
        close(10)
        ipart_old=ipart

     end do
     ! End loop over files

     ! Deallocate local arrays
     deallocate(npart_file,ncum_file)
     
  else     

     !--------------------------------------
     ! Read ASCII initial conditions file
     !--------------------------------------
     
     filetype_loc=r%filetype
     if(.not. r%cosmo)filetype_loc='ascii'
     
     select case (filetype_loc)
        
     case ('grafic')
        ! Particle data read later, nothing is done here
        return

     case ('ascii')

        if(TRIM(r%initfile(r%levelmin)).NE.' ')then
           
           filename=TRIM(r%initfile(r%levelmin))//'/ic_part'
           if(g%myid==1)write(*,*)'Opening file '//TRIM(filename)
           open(10,file=filename,form='formatted')

           ! Figure out starting index for each cpu
           ! as well as total number of particle in file
           jpart=0
           do while (1==1)
              read(10,*,end=101)xx1,xx2,xx3,vv1,vv2,vv3,mm1
              if(ABS(xx1)<r%boxlen/2..AND.ABS(xx2)<r%boxlen/2..AND.ABS(xx3)<r%boxlen/2.)then
                 jpart=jpart+1
              endif
           end do
101        continue
           p%npart_tot=jpart
           if(g%myid==1)write(*,*)'Found npart_tot=',p%npart_tot
           do icpu=1,g%ncpu+1
              start_ind(icpu)=1+((icpu-1)*jpart)/g%ncpu
           end do
           
           rewind(10)
           
           jpart=0
           jpart_loc=0
           indglob=0
           do 
              read(10,*,end=100)xx1,xx2,xx3,vv1,vv2,vv3,mm1
              if(ABS(xx1)<r%boxlen/2..AND.ABS(xx2)<r%boxlen/2..AND.ABS(xx3)<r%boxlen/2.)then
                 jpart=jpart+1
                 indglob=indglob+1
                 if(jpart >= start_ind(g%myid) .and. jpart < start_ind(g%myid+1))then
                    jpart_loc=jpart_loc+1
                    if(jpart_loc>r%npartmax)then
                       write(*,*)'Maximum number of particles incorrect'
                       write(*,*)'npartmax should be greater than',start_ind(2)
                       call clean_stop(g)
                    endif
                    p%xp(jpart_loc,1)=xx1+r%boxlen/2.0
                    p%xp(jpart_loc,2)=xx2+r%boxlen/2.0
                    p%xp(jpart_loc,3)=xx3+r%boxlen/2.0
                    p%vp(jpart_loc,1)=vv1
                    p%vp(jpart_loc,2)=vv2
                    p%vp(jpart_loc,3)=vv3
                    p%mp(jpart_loc  )=mm1
                    p%idp(jpart_loc )=indglob
                    p%levelp(jpart_loc)=r%levelmin
                 end if
              endif
           end do
100        continue
           close(10)
           
        end if
        p%npart=jpart_loc
        
     case ('gadget')
        write(*,*) 'Gadget format not supported '
        call clean_stop(g)
        
     case DEFAULT
        write(*,*) 'Unsupported format file ' // r%filetype
        call clean_stop(g)
        
     end select
  end if

  if(p%npart==0)then
     g%mp_min=-1.0
     return
  endif

  ! Put all particles in levelmin
  p%headp=p%npart+1
  p%tailp=p%npart
  p%headp(r%levelmin)=1
  p%tailp(r%levelmin)=p%npart
        
  ! Compute minimum dark matter particle mass
  g%mp_min=MINVAL(p%mp(1:p%npart))
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(g%mp_min,mp_min_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)  
  g%mp_min=mp_min_all
#endif
  if(g%myid==1)write(*,*)'Minimum particle mass=',g%mp_min
  
  ! Compute maximum number of particles in all processors
  p%npart_max=p%npart
  p%npart_tot=p%npart
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(p%npart,p%npart_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)     
  call MPI_ALLREDUCE(p%npart,p%npart_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif

end subroutine init_part_file_2
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_part_grid_2(r,g,m,p)
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

  ! Allocate default IC grid
  allocate(init_array(1:1,1:1,1:1))
  allocate(init_array_x(1:1,1:1,1:1))

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
        deallocate(init_array,init_array_x)
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
        if(r%multiple)then
           call title(g%myid,nchar)
           if(idim==1)filename=TRIM(r%initfile(ilevel))//'/dir_velcx/ic_velcx.'//TRIM(nchar)
           if(idim==2)filename=TRIM(r%initfile(ilevel))//'/dir_velcy/ic_velcy.'//TRIM(nchar)
           if(idim==3)filename=TRIM(r%initfile(ilevel))//'/dir_velcz/ic_velcz.'//TRIM(nchar)
        else
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
           
        endif
        
        if(g%myid==1)write(*,*)'Reading file '//TRIM(filename)

        if(r%multiple)then

           ! Wait for the token
#ifndef WITHOUTMPI
           if(g%IOGROUPSIZE>0) then
              if (mod(g%myid-1,g%IOGROUPSIZE)/=0) then
                 call MPI_RECV(dummy_io,1,MPI_INTEGER,g%myid-1-1,tagg2,&
                      & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
              end if
           endif
#endif
           ilun=g%myid+10
           open(ilun,file=filename,form='unformatted')
           rewind(ilun)
           read(ilun) ! skip first line
           do i3=1,g%n3(ilevel)
              read(ilun)((init_plane(i1,i2),i1=1,g%n1(ilevel)),i2=1,g%n2(ilevel))
              
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
                 call MPI_SEND(dummy_io,1,MPI_INTEGER,g%myid-1+1,tagg2, &
                      & MPI_COMM_WORLD,info2)
              end if
           endif
#endif
        else
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

        endif

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
        
  ! Compute minimum dark matter particle mass
  g%mp_min=MINVAL(p%mp(1:p%npart))
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(g%mp_min,mp_min_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  g%mp_min=mp_min_all
#endif
  if(g%myid==1)write(*,*)'Minimum particle mass=',g%mp_min

  ! Compute maximum and total number of particles
  p%npart_max=p%npart
  p%npart_tot=p%npart
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(p%npart,p%npart_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(p%npart,p%npart_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif

end subroutine init_part_grid_2
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

