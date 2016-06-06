!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_part_file
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------
  ! Allocate particle-based arrays.
  ! Read particles positions and velocities from various files
  ! including gadget, ascii or restart files.
  ! grafic initial conditions are performed after the AMR grid 
  ! has been constructed.
  !------------------------------------------------------------
  integer::dummyint
  integer::npart2,ndim2,ncpu2
  integer::ipart,jpart,jpart_loc,ipart_old,ilevel,idim
  integer::i,igrid,ngrid,iskip,nsink
  integer::ind,ix,iy,iz,ilun,info,icpu,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count,indglob,npart_new
  real(dp)::dx,xx1,xx2,xx3,vv1,vv2,vv3,mm1,ll1,ll2,ll3
  real(dp)::scale,dx_loc,rr,rmax,dx_min,mp_min_all
  integer::ncode,bit_length,temp,ncpu_file
  integer::ileft,iright,nrest,ipos
  integer(i8b)::npart_tot_file
  integer(i8b)::istart,iend,nleft,nright

  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::isp
  integer(i8b),allocatable,dimension(:)::isp8
  integer,allocatable,dimension(:)::npart_file
  integer(i8b),allocatable,dimension(:)::ncum_file

  integer,dimension(1:ncpu+1)::start_ind
  logical::error,keep_part,eof,jumped,read_pos=.false.,ok
  character(LEN=80)::filename,filename_x
  character(LEN=80)::fileloc,file_part
  character(LEN=20)::filetype_loc
  character(LEN=5)::nchar,ncharcpu

  if(verbose)write(*,*)'Entering init_part_from_file'

  !------------------------------
  ! Allocate particle variables
  !------------------------------
  allocate(xp    (npartmax,ndim))
  allocate(vp    (npartmax,ndim))
  allocate(mp    (npartmax))
  allocate(levelp(npartmax))
  allocate(idp   (npartmax))
  allocate(sortp (npartmax))
  allocate(workp (npartmax))
#ifdef OUTPUT_PARTICLE_POTENTIAL
  allocate(phip  (npartmax))
#endif
  
  ! Allocate pointers to particle levels
  allocate(headp(levelmin:nlevelmax))
  allocate(tailp(levelmin:nlevelmax))

  ! No particle just yet
  headp=1
  tailp=0
  
  if(nrestart>0)then

     !----------------------------------------
     ! In case of restart, read RAMSES file
     !----------------------------------------

     ! Read particle files header
     call title(nrestart,nchar)
     fileloc='output_'//TRIM(nchar)//'/part_header.txt'
     call input_header(fileloc,npart_tot_file,ncpu_file)
     if(myid==1)then
     if(myid==1)write(*,'(" Restart snapshot has ",I8," particles")')npart_tot_file
     endif

     ! Compute new local particle number
     npart_tot=npart_tot_file
     npart=npart_tot_file/ncpu
     nleft=(myid-1)*npart
     nright=myid*npart
     nrest=npart_tot_file-npart*ncpu
     if(myid.LE.nrest)then
        npart=npart+1
        nright=nright+myid
        nleft=nleft+(myid-1)
     endif

     ! Compute number of particles in each file
     allocate(npart_file(1:ncpu_file))
     allocate(ncum_file(0:ncpu_file))
     if(myid==1)then
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
              xp(ipart,idim)=xdp(i)
           end do
        end do

        ! Read velocities
        do idim=1,ndim
           ipos=9+8*(ndim+idim-1)*npart_file(icpu)+8*(istart-1)
           read(10,POS=ipos)xdp
           ipart=ipart_old
           do i=istart,iend
              ipart=ipart+1
              vp(ipart,idim)=xdp(i)
           end do
        end do

        ! Read masses
        ipos=9+8*(ndim+ndim)*npart_file(icpu)+8*(istart-1)
        read(10,POS=ipos)xdp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           mp(ipart)=xdp(i)
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
           idp(ipart)=isp8(i)
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
           levelp(ipart)=isp(i)
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
     
     filetype_loc=filetype
     if(.not. cosmo)filetype_loc='ascii'
     
     select case (filetype_loc)
        
     case ('grafic')
        ! Particle data read later, nothing is done here
        return

     case ('ascii')

        if(TRIM(initfile(levelmin)).NE.' ')then
           
           filename=TRIM(initfile(levelmin))//'/ic_part'
           if(myid==1)write(*,*)'Opening file '//TRIM(filename)
           open(10,file=filename,form='formatted')

           ! Figure out starting index for each cpu
           ! as well as total number of particle in file
           jpart=0
           do while (1==1)
              read(10,*,end=101)xx1,xx2,xx3,vv1,vv2,vv3,mm1
              if(ABS(xx1)<boxlen/2..AND.ABS(xx2)<boxlen/2..AND.ABS(xx3)<boxlen/2.)then
                 jpart=jpart+1
              endif
           end do
101        continue
           npart_tot=jpart
           if(myid==1)write(*,*)'Found npart_tot=',npart_tot
           do icpu=1,ncpu+1
              start_ind(icpu)=1+((icpu-1)*jpart)/ncpu
           end do
           
           rewind(10)
           
           jpart=0
           jpart_loc=0
           indglob=0
           do 
              read(10,*,end=100)xx1,xx2,xx3,vv1,vv2,vv3,mm1
              if(ABS(xx1)<boxlen/2..AND.ABS(xx2)<boxlen/2..AND.ABS(xx3)<boxlen/2.)then
                 jpart=jpart+1
                 indglob=indglob+1
                 if(jpart >= start_ind(myid) .and. jpart < start_ind(myid+1))then
                    jpart_loc=jpart_loc+1
                    if(jpart_loc>npartmax)then
                       write(*,*)'Maximum number of particles incorrect'
                       write(*,*)'npartmax should be greater than',start_ind(2)
                       call clean_stop
                    endif
                    xp(jpart_loc,1)=xx1+boxlen/2.0
                    xp(jpart_loc,2)=xx2+boxlen/2.0
                    xp(jpart_loc,3)=xx3+boxlen/2.0
                    vp(jpart_loc,1)=vv1
                    vp(jpart_loc,2)=vv2
                    vp(jpart_loc,3)=vv3
                    mp(jpart_loc  )=mm1
                    idp(jpart_loc )=indglob
                    levelp(jpart_loc)=levelmin
                 end if
              endif
           end do
100        continue
           close(10)
           
        end if
        npart=jpart_loc
        
     case ('gadget')
        write(*,*) 'Gadget format not supported '
        call clean_stop
        
     case DEFAULT
        write(*,*) 'Unsupported format file ' // filetype
        call clean_stop
        
     end select
  end if

  if(npart==0)return

  ! Put all particles in levelmin
  headp=npart+1
  tailp=npart
  headp(levelmin)=1
  tailp(levelmin)=npart
        
  ! Compute minimum dark matter particle mass
  mp_min=MINVAL(mp(1:npart))
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(mp_min,mp_min_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)  
  mp_min=mp_min_all
#endif
  if(myid==1)write(*,*)'Minimum particle mass=',mp_min
  
  ! Compute maximum number of particles in all processors
  npart_max=npart
  npart_tot=npart
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npart,npart_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)     
  call MPI_ALLREDUCE(npart,npart_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif

end subroutine init_part_file

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

subroutine init_part_grid
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------
  ! Initialize particle position from grid cell positions.
  ! Read particles displacement and velocities from grafic files.
  !------------------------------------------------------------
  integer::dummyint
  integer::npart2,ndim2,ncpu2
  integer::ipart,jpart,jpart_loc,ipart_old,ilevel,idim
  integer::i,igrid,ngrid,iskip,nsink
  integer::ind,ix,iy,iz,ilun,info,icpu,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count,indglob,npart_new
  integer,parameter::tagg=1109,tagg2=1110,tagg3=1111
  integer::dummy_io,info2
  real(dp)::dx,dx_loc,xx1,xx2,xx3,xn1,xn2,xn3,vv1,vv2,vv3
  real(kind=8)::dispmax=0.0,dispall,mp_min_all
  real(kind=4),allocatable,dimension(:,:)::init_plane,init_plane_x
  real(dp),allocatable,dimension(:,:,:)::init_array,init_array_x
  character(LEN=80)::filename,filename_x
  character(LEN=5)::nchar,ncharcpu
  logical::ok,error,keep_part,read_pos=.false.

  if(verbose)write(*,*)'Entering init_part_from_grid'
  if(TRIM(filetype).NE.'grafic')return
  if(nrestart>0)return

  !----------------------------------------------------
  ! Reading initial conditions GRAFIC2 multigrid arrays
  !----------------------------------------------------
  ipart=0
  ! Loop over initial condition levels
  do ilevel=levelmin,nlevelmax
     
     if(initfile(ilevel)==' ')cycle
     
     ! Mesh size at level ilevel in normalised units
     dx=1.0/2**ilevel

     ! Mesh size at level ilevel in code units
     dx_loc=boxlen*dx
     
     !--------------------------------------------------------------
     ! First step: compute level boundaries and particle positions
     !--------------------------------------------------------------
     i1_min=n1(ilevel)+1; i1_max=0
     i2_min=n2(ilevel)+1; i2_max=0
     i3_min=n3(ilevel)+1; i3_max=0
     ipart_old=ipart

     ! Loop over grids
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           ! Coordinates in normalised units (between 0 and 1)
           xx1=(2*grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx
           xx2=(2*grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx
           xx3=(2*grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx
           ! Compute integer coordinates in the frame of the file
           i1=int((xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel))+1
           i2=int((xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel))+1
           i3=int((xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel))+1
           ! Compute min and max
           i1_min=MIN(i1_min,i1); i1_max=MAX(i1_max,i1)
           i2_min=MIN(i2_min,i2); i2_max=MAX(i2_max,i2)
           i3_min=MIN(i3_min,i3); i3_max=MAX(i3_max,i3)
           ! Create particle
           keep_part=.NOT.grid(igrid)%refined(ind)
           if(keep_part)then
              ipart=ipart+1
              if(ipart>npartmax)then
                 write(*,*)'Maximum number of particles incorrect'
                 write(*,*)'npartmax should be greater than',ipart
                 call clean_stop
              endif
              if(ndim>0)xp(ipart,1)=xx1
              if(ndim>1)xp(ipart,2)=xx2
              if(ndim>2)xp(ipart,3)=xx3
              mp(ipart)=0.5d0**(3*ilevel)*(1.0d0-omega_b/omega_m)
              idp(ipart)=1+(i1-1)+(i2-1)*n1(ilevel)+(i3-1)*n1(ilevel)*n2(ilevel)
           endif
        end do
     end do

     ! Check that all grids are within initial condition region
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
     
     !---------------------------------------------------------------------
     ! Second step: read initial condition file and set particle velocities
     !---------------------------------------------------------------------
     ! Allocate initial conditions array
     if(noct(ilevel)>0)then
        allocate(init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
        allocate(init_array_x(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
        init_array=0d0
        init_array_x=0d0
     end if
     allocate(init_plane(1:n1(ilevel),1:n2(ilevel)))
     allocate(init_plane_x(1:n1(ilevel),1:n2(ilevel)))

     ! Loop over input variables
     do idim=1,ndim
        
        ! Read dark matter initial displacement field
        if(multiple)then
           call title(myid,nchar)
           if(idim==1)filename=TRIM(initfile(ilevel))//'/dir_velcx/ic_velcx.'//TRIM(nchar)
           if(idim==2)filename=TRIM(initfile(ilevel))//'/dir_velcy/ic_velcy.'//TRIM(nchar)
           if(idim==3)filename=TRIM(initfile(ilevel))//'/dir_velcz/ic_velcz.'//TRIM(nchar)
        else
           if(idim==1)filename=TRIM(initfile(ilevel))//'/ic_velcx'
           if(idim==2)filename=TRIM(initfile(ilevel))//'/ic_velcy'
           if(idim==3)filename=TRIM(initfile(ilevel))//'/ic_velcz'

           if(idim==1)filename_x=TRIM(initfile(ilevel))//'/ic_poscx'
           if(idim==2)filename_x=TRIM(initfile(ilevel))//'/ic_poscy'
           if(idim==3)filename_x=TRIM(initfile(ilevel))//'/ic_poscz'
           
           INQUIRE(file=filename_x,exist=ok)
           if(.not.ok)then
              read_pos = .false.
           else
              read_pos = .true.
              if(myid==1)write(*,*)'Reading file '//TRIM(filename_x)
           end if
           
        endif
        
        if(myid==1)write(*,*)'Reading file '//TRIM(filename)

        if(multiple)then

           ! Wait for the token
#ifndef WITHOUTMPI
           if(IOGROUPSIZE>0) then
              if (mod(myid-1,IOGROUPSIZE)/=0) then
                 call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tagg2,&
                      & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
              end if
           endif
#endif
           ilun=myid+10
           open(ilun,file=filename,form='unformatted')
           rewind(ilun)
           read(ilun) ! skip first line
           do i3=1,n3(ilevel)
              read(ilun)((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
              
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
                 call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tagg2, &
                      & MPI_COMM_WORLD,info2)
              end if
           endif
#endif
        else
           if(myid==1)then
              open(10,file=filename,form='unformatted')
              rewind 10
              read(10) ! skip first line
           end if
           do i3=1,n3(ilevel)
              if(myid==1)then
                 if(debug.and.mod(i3,10)==0)write(*,*)'Reading plane ',i3
                 read(10)((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
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

           if(read_pos) then
              if(myid==1)then
                 open(10,file=filename_x,form='unformatted')
                 rewind 10
                 read(10) ! skip first line
              end if
              do i3=1,n3(ilevel)
                 if(myid==1)then
                    if(debug.and.mod(i3,10)==0)write(*,*)'Reading plane ',i3
                    read(10)((init_plane_x(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                 else
                    init_plane_x=0.0
                 endif
                 buf_count=n1(ilevel)*n2(ilevel)
#ifndef WITHOUTMPI
                 call MPI_BCAST(init_plane_x,buf_count,MPI_REAL,0,MPI_COMM_WORLD,info)
#endif
                 if(noct(ilevel)>0)then
                    if(i3.ge.i3_min.and.i3.le.i3_max)then
                       init_array_x(i1_min:i1_max,i2_min:i2_max,i3) = &
                            & init_plane_x(i1_min:i1_max,i2_min:i2_max)
                    end if
                 endif
              end do
              if(myid==1)close(10)
           end if

        endif

        if(noct(ilevel)>0)then
           ! Rescale initial displacement field to code units
           init_array=dfact(ilevel)*dx/dxini(ilevel)*init_array/vfact(ilevel)
           if(read_pos)then
              init_array_x = init_array_x/boxlen_ini
           endif
           ipart=ipart_old
           ! Loop over grids
           do igrid=head(ilevel),tail(ilevel)
              ! Loop over cells
              do ind=1,twotondim
                 ! Coordinates in normalised units (between 0 and 1)
                 xx1=(2*grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx
                 xx2=(2*grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx
                 xx3=(2*grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx
                 ! Compute integer coordinates in the frame of the file
                 i1=int((xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel))+1
                 i2=int((xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel))+1
                 i3=int((xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel))+1
                 ! Compute particle velocity and identity
                 keep_part=.NOT.grid(igrid)%refined(ind)
                 if(keep_part)then
                    ipart=ipart+1
                    vp(ipart,idim)=init_array(i1,i2,i3)
                    if(.not. read_pos)then
                       dispmax=max(dispmax,abs(init_array(i1,i2,i3)/dx))
                    else
                       if(idim==1)xp(ipart,idim)=xx1+init_array_x(i1,i2,i3)
                       if(idim==2)xp(ipart,idim)=xx2+init_array_x(i1,i2,i3)
                       if(idim==3)xp(ipart,idim)=xx3+init_array_x(i1,i2,i3)
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
     if(noct(ilevel)>0)then
        deallocate(init_array,init_array_x)
     end if
     deallocate(init_plane,init_plane_x)
     
     if(debug)write(*,*)'npart=',ipart,'/',npartmax,' for PE=',myid

  end do
  ! End loop over levels

  ! Initial particle number
  npart=ipart

  ! Move particle according to Zeldovich approximation
  if(.not. read_pos)then
     do ipart=1,npart
        xp(ipart,1:ndim)=xp(ipart,1:ndim)+vp(ipart,1:ndim)
     enddo
  endif

  ! Scale displacement to velocity
  do ipart=1,npart
     vp(ipart,1:ndim)=vfact(1)*vp(ipart,1:ndim)
  end do

  ! Periodic box
  do ipart=1,npart
#if NDIM>0
     if(xp(ipart,1)< 0.0d0 )xp(ipart,1)=xp(ipart,1)+boxlen
     if(xp(ipart,1)>=boxlen)xp(ipart,1)=xp(ipart,1)-boxlen
#endif
#if NDIM>1
     if(xp(ipart,2)< 0.0d0 )xp(ipart,2)=xp(ipart,2)+boxlen
     if(xp(ipart,2)>=boxlen)xp(ipart,2)=xp(ipart,2)-boxlen
#endif
#if NDIM>2
     if(xp(ipart,3)< 0.0d0 )xp(ipart,3)=xp(ipart,3)+boxlen
     if(xp(ipart,3)>=boxlen)xp(ipart,3)=xp(ipart,3)-boxlen
#endif
  end do

  ! Compute particle initial level
  do ipart=1,npart
     levelp(ipart)=levelmin
  end do

  ! Put all particles in levelmin
  headp=npart+1
  tailp=npart
  headp(levelmin)=1
  tailp(levelmin)=npart
        
  ! Compute minimum dark matter particle mass
  mp_min=MINVAL(mp(1:npart))
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(mp_min,mp_min_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  mp_min=mp_min_all
#endif
  if(myid==1)write(*,*)'Minimum particle mass=',mp_min

  ! Compute maximum and total number of particles
  npart_max=npart
  npart_tot=npart
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npart,npart_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(npart,npart_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif

end subroutine init_part_grid

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

