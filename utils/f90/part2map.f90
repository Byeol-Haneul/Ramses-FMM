program part2map
  !--------------------------------------------------------------------------
  ! This code projects RAMSES data onto a map.
  ! This is based in the MINI-RAMSES prototype.
  ! R. Teyssier, Princeton, February 2nd 2023
  !--------------------------------------------------------------------------
  implicit none

  integer,parameter::flen=200
  
  integer::npart,nfile
  integer::nx_sample=128,ny_sample=128,nx,ny
  integer::i,j,idim,jdim,icpu
  integer::ix,iy,ixp1,iyp1
  integer::npart_file, ndim_file, ipos, npart_actual

  real(KIND=8)::xmin=-1,xmax=-1,ymin=-1,ymax=-1,zmin=-1,zmax=-1
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax
  real(KIND=8)::dx,dy,dz,xx,yy,zz,mtot
  real(KIND=8)::jxin=0, jyin=0, jzin=0
  real(KIND=8)::xcenter, ycenter, zcenter, z_coord, y_coord, x_coord
  real(kind=8)::jx, jy, jz, kx, ky, kz, lx, ly, lz, ddx, ddy, dex, dey

  character(LEN=1)::proj='z'
  character(LEN=5)::nchar,ncharcpu
  character(LEN=flen)::nomfich,repository,prefix='part',outfich
  character(LEN=flen)::file_amr,file_hydro
  character(LEN=flen)::filetype='bin'
  character(LEN=flen)::file_part

  logical::ok,ok_part,ok_cell,do_max,do_grav,do_peak,do_rt
  logical::backup_file=.false.,check_ramses_exist
  logical::rotation=.false., periodic=.false., sideon=.false.

  real(KIND=8),dimension(:,:),allocatable::map
  real(KIND=4),dimension(:,:),allocatable::toto
  real(KIND=4),dimension(:,:),allocatable::x
  real(KIND=4),dimension(:),allocatable::m

  type params
     integer::ndim
     integer::ncpu
     integer::nvar
     integer::levelmin
     integer::nlevelmax
     real(kind=8)::boxlen
     real(kind=8)::t
     real(kind=8)::texp
     real(kind=8)::aexp
     real(kind=8)::gamma
     real(kind=8)::unit_d
     real(kind=8)::unit_l
     real(kind=8)::unit_t
  end type params
  
  type(params)::p
  
  write(*,*)'Starting part2map'

  !------------------------------
  ! Read parameter and info files
  !------------------------------

  ! Read part2map parameters
  call read_params

  ! Check that all files exist
  if(.NOT. check_ramses_exist(repository,prefix))then
     write(*,*)'Repository '//TRIM(repository)//' incomplete.'
     write(*,*)'Stopping.'
     stop
  endif
  if(index(repository,'output')==0)backup_file=.true.

  ! Read RAMSES params
  call read_ramses_params
  write(*,*)'time =',real(p%t,kind=4)
  
  ! Read RAMSES info
  call read_info
  
  ! Read PART header file
  call read_part_header
  write(*,*)'npart=',npart
  write(*,*)'nfile=',nfile
  
  !-----------------
  ! Set up geometry 
  !-----------------
  if(xmin<0)xmin=0
  if(xmax<0)xmax=p%boxlen
  if(ymin<0)ymin=0
  if(ymax<0)ymax=p%boxlen
  if(zmin<0)zmin=0
  if(zmax<0)zmax=p%boxlen
  if(abs(jxin)+abs(jyin)+abs(jzin)==0)then
     if (proj=='x')then
        idim=2
        jdim=3
        xxmin=ymin ; xxmax=ymax
        yymin=zmin ; yymax=zmax
        zzmin=xmin ; zzmax=xmax
     else if (proj=='y') then
        idim=1
        jdim=3
        xxmin=xmin ; xxmax=xmax
        yymin=zmin ; yymax=zmax
        zzmin=ymin ; zzmax=ymax
     else
        idim=1
        jdim=2
        xxmin=xmin ; xxmax=xmax
        yymin=ymin ; yymax=ymax
        zzmin=zmin ; zzmax=zmax
     end if
     rotation=.false.
  else
     xcenter=0.5*(xmin+xmax)
     ycenter=0.5*(ymin+ymax)
     zcenter=0.5*(zmin+zmax)
     jx=jxin/sqrt(jxin**2+jyin**2+jzin**2)
     jy=jyin/sqrt(jxin**2+jyin**2+jzin**2)
     jz=jzin/sqrt(jxin**2+jyin**2+jzin**2)
     kx=0
     ky=-jz
     kz=jy
     lx=jy*kz-jz*ky
     ly=jz*kx-jx*kz
     lz=jx*ky-jy*kx
     rotation=.true.
  endif
  dx=(xxmax-xxmin)/dble(nx)
  dy=(yymax-yymin)/dble(ny)
  
  !-----------------------------------------------
  ! Compute projected variables
  !----------------------------------------------
  ! Loop over processor files
  npart_actual=0
  allocate(map(0:nx,0:ny))
  do icpu=1,nfile
     call title(icpu,ncharcpu)

     ! Prepare reading the AMR file
     file_part=TRIM(repository)//'/'//TRIM(prefix)//'.'//TRIM(ncharcpu)

     open(unit=10,file=file_part,access="stream",action="read",form='unformatted')
     ipos=1
     read(10,POS=ipos)ndim_file
     ipos=5
     read(10,POS=ipos)npart_file

     allocate(x(1:npart_file,1:ndim_file))
     allocate(m(1:npart_file))
     ipos=9
     read(10,POS=ipos)x
     ipos=9+4*npart_file*2*ndim_file
     read(10,POS=ipos)m

     do i=1,npart_file
        ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
             &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
             &   x(i,3)>=zmin.and.x(i,3)<=zmax)

        if(ok_part)then
           npart_actual=npart_actual+1
           if(rotation)then
              ! Perform rotation
              xx=x(i,1)-xcenter
              yy=x(i,2)-ycenter
              zz=x(i,3)-zcenter
              z_coord=xx*jx+yy*jy+zz*jz+zcenter
              y_coord=xx*kx+yy*ky+zz*kz+ycenter
              x_coord=xx*lx+yy*ly+zz*lz+xcenter
              ddx=(x_coord-xmin)/dx
              ddy=(y_coord-ymin)/dy
              if(sideon)then
                 ddx=(x_coord-xmin)/dx
                 ddy=(z_coord-zmin)/dy
              endif
           else
              ddx=(x(i,idim)-xxmin)/dx
              ddy=(x(i,jdim)-yymin)/dy
           endif
           ! CIC mass deposition 
           ix=ddx
           iy=ddy
           ddx=ddx-ix
           ddy=ddy-iy
           dex=1.0-ddx
           dey=1.0-ddy
           if(periodic)then
              if(ix<0)ix=ix+nx
              if(ix>=nx)ix=ix-nx
              if(iy<0)iy=iy+ny
              if(iy>=ny)iy=iy-ny
           endif
           ixp1=ix+1
           iyp1=iy+1
           if(periodic)then
              if(ixp1<0)ixp1=ixp1+nx
              if(ixp1>=nx)ixp1=ixp1-nx
              if(iyp1<0)iyp1=iyp1+ny
              if(iyp1>=ny)iyp1=iyp1-ny
           endif
           if(ix>=0.and.ix<nx.and.iy>=0.and.iy<ny.and.ddx>=0.and.ddy>=0)then
              map(ix  ,iy  )=map(ix  ,iy  )+m(i)*dex*dey
              map(ix  ,iyp1)=map(ix  ,iyp1)+m(i)*dex*ddy
              map(ixp1,iy  )=map(ixp1,iy  )+m(i)*ddx*dey
              map(ixp1,iyp1)=map(ixp1,iyp1)+m(i)*ddx*ddy
              mtot=mtot+m(i)
           endif
        end if

     end do

     deallocate(x,m)
     ! Close file
     close(10)
     close(11)

  end do
  ! End loop over cpu

  write(*,*)'Data read and projected.'
  write(*,*)'Total number of part=',npart_actual
  write(*,*)'Total deposited mass=',mtot

  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Writing data to '//TRIM(nomfich)
  if(periodic)then
     allocate(toto(0:nx-1,0:ny-1))
     toto=map(0:nx-1,0:ny-1)
  else
     allocate(toto(0:nx,0:ny))
     toto=map(0:nx,0:ny)
  endif
  ! Binary format (to be read by ramses utilities)
  if(TRIM(filetype).eq.'bin')then
     open(unit=10,file=nomfich,form='unformatted')
     if(periodic)then
        write(10)p%t,xmax-xmin,ymax-ymin,zmax-zmin
        write(10)nx,ny
        write(10)toto
        write(10)xmin,xmax
        write(10)ymin,ymax
     else
        write(10)p%t,xmax-xmin,ymax-ymin,zmax-zmin
        write(10)nx+1,ny+1
        write(10)toto
        write(10)xmin,xmax
        write(10)ymin,ymax
     endif
     close(10)
  endif
  ! Ascii format (to be read by gnuplot)
  if(TRIM(filetype).eq.'ascii')then
     open(unit=10,file=nomfich,form='formatted')
     if(periodic)then
        do j=0,ny-1
           do i=0,nx-1
              xx=xmin+(dble(i)+0.5)/dble(nx)*(xmax-xmin)
              yy=ymin+(dble(j)+0.5)/dble(ny)*(ymax-ymin)
              write(10,*)xx,yy,toto(i,j)
           end do
           write(10,*) " "
        end do
     else
        do j=0,ny
           do i=0,nx
              xx=xmin+dble(i)/dble(nx)*(xmax-xmin)
              yy=ymin+dble(j)/dble(ny)*(ymax-ymin)
              write(10,*)xx,yy,toto(i,j)
           end do
           write(10,*) " "
        end do
     endif
     close(10)
  endif

contains
  
  subroutine read_params
    
    implicit none
    
    integer       :: i,n
    integer       :: iargc
    character(len=4)   :: opt
    character(len=128) :: arg
    LOGICAL       :: bad, ok
    real(kind=8)  :: xcen=0.5,ycen=0.5,zcen=0.5,rad=0
    
    n = iargc()
    if (n < 4) then
       print *, 'usage: part2map -inp  input_dir'
       print *, '                -out  output_file'
       print *, '               [-dir axis] '
       print *, '               [-jx  jxin] '
       print *, '               [-jy  jyin] '
       print *, '               [-jz  jzin] '
       print *, '               [-xce xcen] '
       print *, '               [-yce ycen] '
       print *, '               [-zce zcen] '
       print *, '               [-rad rad] '
       print *, '               [-xmi xmin] '
       print *, '               [-xma xmax] '
       print *, '               [-ymi ymin] '
       print *, '               [-yma ymax] '
       print *, '               [-zmi zmin] '
       print *, '               [-zma zmax] '
       print *, '               [-nx  nx] '
       print *, '               [-ny  ny] '
       print *, '               [-fil filetype] '
       print *, '               [-pre prefix] '
       print *, 'ex: part2map -inp output_00001 -out map.dat'// &
            &   ' -dir z -xmi 0.1 -xma 0.7'
       print *, ' '
       stop
    end if
    
    do i = 1,n,2
       call getarg(i,opt)
       if (i == n) then
          print '("option ",a2," has no argument")', opt
          stop 2
       end if
       call getarg(i+1,arg)
       select case (opt)
       case ('-inp')
          repository = trim(arg)
       case ('-out')
          outfich = trim(arg)
       case ('-dir')
          proj = trim(arg) 
       case ('-xce')
          read (arg,*) xcen
       case ('-yce')
          read (arg,*) ycen
       case ('-zce')
          read (arg,*) zcen
       case ('-rad')
          read (arg,*) rad
       case ('-xmi')
          read (arg,*) xmin
       case ('-xma')
          read (arg,*) xmax
       case ('-ymi')
          read (arg,*) ymin
       case ('-yma')
          read (arg,*) ymax
       case ('-zmi')
          read (arg,*) zmin
       case ('-zma')
          read (arg,*) zmax
       case ('-nx')
          read (arg,*) nx_sample
       case ('-ny')
          read (arg,*) ny_sample
       case ('-fil')
          read (arg,*) filetype
       case ('-pre')
          read (arg,*) prefix
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
    
    if(rad>0)then
       xmin=xcen-rad
       xmax=xcen+rad
       ymin=ycen-rad
       ymax=ycen+rad
       zmin=zcen-rad
       zmax=zcen+rad
    endif

    if((abs(jxin)+abs(jyin)+abs(jzin))>0.and.rad==0)then
       write(*,*)'rotations only allowed with spherical selection'
       stop
    endif

    nx=nx_sample
    ny=ny_sample

    return
    
  end subroutine read_params  
  !================================================================
  !================================================================
  !================================================================
  !================================================================
  subroutine read_ramses_params
    !-----------------------------------------------
    ! Read RAMSES parameters file
    !-----------------------------------------------
    character(LEN=128)::nomfich
    integer::ilun,ilevel,noutput,skip
    integer::nfile1,ncpu1
    
    nomfich=TRIM(repository)//'/params.bin'
    
    ilun=10
    open(unit=ilun,file=nomfich,access="stream",action="read",form='unformatted')
    read(ilun,POS=1)nfile1
    read(ilun,POS=5)ncpu1
    read(ilun,POS=9)p%ndim
    read(ilun,POS=13)p%levelmin
    read(ilun,POS=17)p%nlevelmax
    read(ilun,POS=21)p%boxlen
    read(ilun,POS=29)noutput
    skip=4*(11+4*noutput)+1
    read(ilun,POS=skip)p%t
    skip=skip+4*(2+4*p%nlevelmax+2+2*17)
    read(ilun,POS=skip)p%gamma
    close(ilun)

    if(backup_file)then
       p%ncpu=ncpu1
    else
       p%ncpu=nfile1
    endif

  end subroutine read_ramses_params
  !================================================================
  !================================================================
  !================================================================
  !================================================================
  subroutine read_part_header
    !-----------------------------------------------
    ! Read PART header file
    !-----------------------------------------------
    character(LEN=128)::nomfich,fields
    integer::ilun
    
    nomfich=TRIM(repository)//'/'//TRIM(prefix)//'_header.txt'
    
    ilun=10
    open(unit=ilun,file=nomfich,form='formatted')
    read(ilun,*)
    read(ilun,*)npart
    read(ilun,*)
    read(ilun,*)nfile
    read(ilun,*)
    read(ilun,*)fields

  end subroutine read_part_header
  !================================================================
  !================================================================
  !================================================================
  !================================================================
  subroutine read_info
    !-----------------------------------------------
    ! Read ramses info file
    !-----------------------------------------------
    character(LEN=128)::nomfich
    character(LEN=80)::GMGM
    integer::ilun,nfile1,ncpu1,ndim1,levelmin1,levelmax1
    real(kind=8)::boxlen,t,omega_m,omega_b,omega_l,omega_k,gamma,h0

    nomfich=TRIM(repository)//'/info.txt'

    ilun=10
    open(unit=ilun,file=nomfich,form='formatted')
    read(ilun,'(A13,I11)')GMGM,nfile1
    read(ilun,'(A13,I11)')GMGM,ncpu1
    read(ilun,'(A13,I11)')GMGM,ndim1
    read(ilun,'(A13,I11)')GMGM,levelmin1
    read(ilun,'(A13,I11)')GMGM,levelmax1
    read(ilun,*) ! ngridmax
    read(ilun,*) ! nstep_coarse
    read(ilun,*)

    read(ilun,'(A13,E23.15)')GMGM,boxlen
    read(ilun,'(A13,E23.15)')GMGM,t
    read(ilun,'(A13,E23.15)')GMGM,p%texp
    read(ilun,'(A13,E23.15)')GMGM,p%aexp
    read(ilun,'(A13,E23.15)')GMGM,h0
    read(ilun,'(A13,E23.15)')GMGM,omega_m
    read(ilun,'(A13,E23.15)')GMGM,omega_l
    read(ilun,'(A13,E23.15)')GMGM,omega_k
    read(ilun,'(A13,E23.15)')GMGM,omega_b
    read(ilun,'(A13,E23.15)')GMGM,gamma
    read(ilun,'(A13,E23.15)')GMGM,p%unit_l
    read(ilun,'(A13,E23.15)')GMGM,p%unit_d
    read(ilun,'(A13,E23.15)')GMGM,p%unit_t
    read(ilun,*)

    close(ilun)

  end subroutine read_info
  !================================================================
  !================================================================
  !================================================================
  !================================================================
end program part2map

function check_ramses_exist(repository,prefix)
  !-----------------------------------------------
  ! Check that RAMSES files are there
  !-----------------------------------------------
  logical::check_ramses_exist
  character(len=80)::repository,prefix
  character(LEN=128)::nomfich_part  
  check_ramses_exist=.true.
  nomfich_part=TRIM(repository)//'/'//TRIM(prefix)//'.00001'
  inquire(file=nomfich_part, exist=check_ramses_exist) ! verify input file 
end function check_ramses_exist
!================================================================
!================================================================
!================================================================
!================================================================

!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif

end subroutine title
