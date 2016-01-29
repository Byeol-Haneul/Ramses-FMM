subroutine init_part
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------
  ! Allocate particle-based arrays.
  ! Read particles positions and velocities from grafic files
  !------------------------------------------------------------
  integer::dummyint
  integer::npart2,ndim2,ncpu2
  integer::ipart,jpart,jpart_loc,ipart_old,ilevel,idim
  integer::i,igrid,ngrid,iskip,nsink
  integer::ind,ix,iy,iz,ilun,info,icpu,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count,indglob,npart_new
  real(dp)::dx,xx1,xx2,xx3,vv1,vv2,vv3,mm1,ll1,ll2,ll3
  real(dp)::scale,dx_loc,rr,rmax,dx_min
  integer::ncode,bit_length,temp
  real(kind=8)::bscale
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nvector)::ind_grid,ind_cell,cc,ii
  integer(i8b),dimension(1:ncpu)::npart_cpu,npart_all
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::isp
  integer(i8b),allocatable,dimension(:)::isp8
  logical,allocatable,dimension(:)::nb
  real(kind=4),allocatable,dimension(:,:)::init_plane,init_plane_x
  real(dp),allocatable,dimension(:,:,:)::init_array,init_array_x
  real(kind=8),dimension(1:nvector,1:3)::xx,vv,xs
  real(dp),dimension(1:nvector,1:3)::xx_dp
  integer,dimension(1:nvector)::ixx,iyy,izz
  real(qdp),dimension(1:nvector)::order
  real(kind=8),dimension(1:nvector)::mm
  real(kind=8)::dispmax=0.0,dispall
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:3)::centerofmass

  integer::ibuf,tag=101,tagf=102,tagu=102
  integer::countsend,countrecv
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE,2*ncpu)::statuses
  integer,dimension(2*ncpu)::reqsend,reqrecv
  integer,dimension(ncpu)::sendbuf,recvbuf
#endif

  integer,dimension(1:ncpu+1)::start_ind
  logical::error,keep_part,eof,jumped,read_pos=.false.,ok
  character(LEN=80)::filename,filename_x
  character(LEN=80)::fileloc
  character(LEN=20)::filetype_loc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering init_part'

  if(allocated(xp))then
     if(verbose)write(*,*)'Initial conditions already set'
     return
  end if

  ! Allocate particle variables
  allocate(xp    (npartmax,ndim))
  allocate(vp    (npartmax,ndim))
  allocate(ap    (npartmax,ndim))
  allocate(mp    (npartmax))
  allocate(levelp(npartmax))
  allocate(idp   (npartmax))
  allocate(part_ref_mask (npartmax))
  allocate(part_hkey(npartmax,0:2))
  allocate(current_state(npartmax))
  allocate(part_ind_permutation(npartmax))
  allocate(part_ind_permutation2(npartmax))
  allocate(part_level_offset(levelmin:nlevelmax + 1))
#ifdef OUTPUT_PARTICLE_POTENTIAL
  stop
  allocate(ptcl_phi(npartmax))
#endif

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
     read(ilun)ncpu2
     read(ilun)ndim2
     read(ilun)npart2
     read(ilun)dummyint
     read(ilun)dummyint
     read(ilun)dummyint
     read(ilun)dummyint
     read(ilun)dummyint
     if(ncpu2.ne.ncpu.or.ndim2.ne.ndim.or.npart2.gt.npartmax)then
        write(*,*)'File part.tmp not compatible'
        write(*,*)'Found   =',ncpu2,ndim2,npart2
        write(*,*)'Expected=',ncpu,ndim,npartmax
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

     ! put all particles to levelmin
     part_level_offset = npart
     part_level_offset(levelmin) = 0
  else     

     filetype_loc=filetype
     if(.not. cosmo)filetype_loc='ascii'
     
     select case (filetype_loc)
        
     case ('grafic')
        write(*,*) 'grafic ICs currently not supported by experimental particle implementation'
        stop

     case ('ascii')

        ! Local particle count
        ipart=0

        if(TRIM(initfile(levelmin)).NE.' ')then
           
           filename=TRIM(initfile(levelmin))//'/ic_part'
           open(10,file=filename,form='formatted')
           indglob=0

           !figure out starting indices for domains
           jpart=0
           do while (1==1)
              read(10,*,end=101)xx1,xx2,xx3,vv1,vv2,vv3,mm1
              jpart=jpart+1
           end do
101        continue

           do icpu=1,ncpu+1
              start_ind(icpu)=1+((icpu-1)*jpart)/ncpu
           end do

           rewind(10)

           jpart=0
           jpart_loc=0
           do 
              read(10,*,end=100)xx1,xx2,xx3,vv1,vv2,vv3,mm1
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
           end do
100        continue
           close(10)

        end if
        npart=jpart_loc

        ! put all particles to levelmin
        part_level_offset = npart
        part_level_offset(levelmin) = 0

        ! Compute total number of particle
        npart_cpu=0; npart_all=0
        npart_cpu(myid)=npart
#ifndef WITHOUTMPI
#ifndef LONGINT
        call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
        call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
        npart_cpu(1)=npart_all(1)
#endif
        do icpu=2,ncpu
           npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
        end do
        if(debug)write(*,*)'npart=',npart,'/',npart_cpu(ncpu),jpart
        
        ! don't support gadget for now...
        !     case ('gadget')
        !     call load_gadget

     case DEFAULT
        write(*,*) 'Unsupported format file ' // filetype
        call clean_stop

     end select
  end if

end subroutine init_part
#define TIME_START(cs) call SYSTEM_CLOCK(COUNT=cs)
#define TIME_END(ce) call SYSTEM_CLOCK(COUNT=ce)
#define TIME_SPENT(cs,ce,cr) REAL((ce-cs)/cr)
