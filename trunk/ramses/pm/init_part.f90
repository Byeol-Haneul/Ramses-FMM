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
  integer::ncode,bit_length,temp

  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::isp
  integer(i8b),allocatable,dimension(:)::isp8

  integer,dimension(1:ncpu+1)::start_ind
  logical::error,keep_part,eof,jumped,read_pos=.false.,ok
  character(LEN=80)::filename,filename_x
  character(LEN=80)::fileloc
  character(LEN=20)::filetype_loc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering init_part'

  ! Allocate particle variables
  allocate(xp    (npartmax,ndim))
  allocate(vp    (npartmax,ndim))
  allocate(mp    (npartmax))
  allocate(levelp(npartmax))
  allocate(idp   (npartmax))
  allocate(sortp (npartmax))
  allocate(workp (npartmax))
#ifdef OUTPUT_PARTICLE_POTENTIAL
  allocate(ptcl_phi(npartmax))
#endif
  
  ! Allocate pointers to particle levels
  allocate(headp(levelmin:nlevelmax))
  allocate(tailp(levelmin:nlevelmax))

  ! No particle just yet
  headp=1
  tailp=0
  
  !--------------------
  ! Read part.tmp file
  !--------------------

  if(nrestart>0)then

     ilun=2*ncpu+myid+10
     call title(nrestart,nchar)
     fileloc='output_'//TRIM(nchar)//'/part_'//TRIM(nchar)//'.out'
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)

     open(unit=ilun,file=fileloc,form='unformatted')
     rewind(ilun)
     read(ilun)ndim2
     read(ilun)npart2
     if(ndim2.ne.ndim.or.npart2.gt.npartmax)then
        write(*,*)'File part.tmp not compatible'
        write(*,*)'Found   =',ndim2,npart2
        write(*,*)'Expected=',ndim,npartmax
        call clean_stop
     end if
     ! Read position
     allocate(xdp(1:npart2))
     do idim=1,ndim
        read(ilun)xdp
        xp(1:npart2,idim)=xdp
     end do
     ! Read velocity
     do idim=1,ndim
        read(ilun)xdp
        vp(1:npart2,idim)=xdp
     end do
     ! Read mass
     read(ilun)xdp
     mp(1:npart2)=xdp
     deallocate(xdp)
     ! Read identity
     allocate(isp8(1:npart2))
     read(ilun)isp8
     idp(1:npart2)=isp8
     deallocate(isp8)
     ! Read level
     allocate(isp(1:npart2))
     read(ilun)isp
     levelp(1:npart2)=isp
     deallocate(isp)
     close(ilun)
     if(debug)write(*,*)'part.tmp read for processor ',myid
     npart=npart2

     ! Put all particles in levelmin
     headp=npart+1
     tailp=npart
     headp(levelmin)=1
     tailp(levelmin)=npart

  else     
  !--------------------------------------
  ! Read particle initial conditions file
  !--------------------------------------

     filetype_loc=filetype
     if(.not. cosmo)filetype_loc='ascii'
     
     select case (filetype_loc)
        
     case ('grafic')
        write(*,*) 'grafic IC: particle will be generated from the grid later.'

     case ('ascii')

        if(TRIM(initfile(levelmin)).NE.' ')then
           
           filename=TRIM(initfile(levelmin))//'/ic_part'
           if(myid==1)write(*,*)' Opening file '//TRIM(filename)
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
           if(myid==1)write(*,*)' Found npart_tot=',npart_tot
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
        if(myid==1)write(*,*)' Read npart=',npart
        
        ! Put all particles in levelmin
        headp=npart+1
        tailp=npart
        headp(levelmin)=1
        tailp(levelmin)=npart
        
        if(debug)write(*,*)'npart=',npart,'/',npart_tot
        
     case ('gadget')
        write(*,*) 'Gadget format not supported '
        call clean_stop
        
     case DEFAULT
        write(*,*) 'Unsupported format file ' // filetype
        call clean_stop
        
     end select
  end if

  ! Compute minimum dark matter particle mass
  mp_min=MINVAL(mp(1:npart))
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(mp_min,mp_min_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)  
  mp_min=mp_min_all
#endif
  if(myid==1)write(*,*)'mass minimum=',mp_min

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
  real(dp)::dx,xx1,xx2,xx3,vv1,vv2,vv3,mm1,ll1,ll2,ll3


end subroutine init_part_grid

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

