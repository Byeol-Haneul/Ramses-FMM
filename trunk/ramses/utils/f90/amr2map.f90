program amr2map
  use amr_commons
  use hydro_commons
  use hilbert
  !--------------------------------------------------------------------------
  ! This code projects RAMSES data onto a map.
  ! This is based in the MINI-RAMSES prototype.
  ! R. Teyssier, Zurich, 15/04/16.
  !--------------------------------------------------------------------------
  implicit none
  integer::n,i,j,k,type=0,domax=0
  integer::ivar,lmax=0
  integer::ilevel,idim,jdim,kdim,ldim,icell
  integer::nlevelmaxs,nlevel
  integer::ind,icpu,ncpu_read

  integer::nx_sample=0,ny_sample=0
  integer::imin,imax,jmin,jmax,kmin,kmax
  integer::nstride,nx_full,ny_full,lmin
  integer::ix,iy,iz,ndom,impi,bit_length,maxdom
  integer::iskip_amr,iskip_hydro,noct_skip,noct_file,noct_tmp
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max

  real(KIND=8)::dxline,weight
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax
  real(KIND=8)::dx,dy,dz,xx,yy,zz
  real(KIND=8)::rho,map
  real(KIND=8)::metmax=0d0

  character(LEN=1)::proj='z'
  character(LEN=5)::nchar,ncharcpu
  character(LEN=128)::nomfich,repository,outfich,mapfiletype='bin'
  character(LEN=128)::file_amr,file_hydro
  logical::ok,ok_part,ok_cell,do_max
  logical::check_ramses_exist

  real(KIND=4),dimension(:,:),allocatable::toto
  real(KIND=8),dimension(:),allocatable::bound_key
  integer,dimension(:),allocatable::cpu_list

  integer,dimension(1:ndim)::ckey,cart_key
  logical,dimension(1:twotondim)::refined
  real(dp),dimension(1:twotondim,1:nvar)::uold

  type level  
     integer::ilevel
     real(KIND=8),dimension(:,:),pointer::map
     real(KIND=8),dimension(:,:),pointer::rho
     integer::imin
     integer::imax
     integer::jmin
     integer::jmax
  end type level

  type(level),dimension(1:100)::mapgrid

  write(*,*)'Starting amr2map'

  ! Read amr2map parameters
  call read_params

  ! Check that all files exist
  if(.NOT. check_ramses_exist(repository))then
     write(*,*)'Repository '//TRIM(repository)//' incomplete.'
     write(*,*)'Stopping.'
     stop
  endif

  ! Read RAMSES params
  call read_ramses_params(repository)
  write(*,*)'time=',t

  !------------------------------
  ! Set up geometry parameters
  !------------------------------
  if(lmax==0)then
     lmax=nlevelmax
  endif
  write(*,*)'Working resolution =',2**lmax
  zzmax=1.0
  zzmin=0.0
  if(ndim>2)then
     if (proj=='x')then
        idim=2
        jdim=3
        kdim=1
        xxmin=ymin ; xxmax=ymax
        yymin=zmin ; yymax=zmax
        zzmin=xmin ; zzmax=xmax
     else if (proj=='y') then
        idim=1
        jdim=3
        kdim=2
        xxmin=xmin ; xxmax=xmax
        yymin=zmin ; yymax=zmax
        zzmin=ymin ; zzmax=ymax
     else
        idim=1
        jdim=2
        kdim=3
        xxmin=xmin ; xxmax=xmax
        yymin=ymin ; yymax=ymax
        zzmin=zmin ; zzmax=zmax
     end if
  else
     idim=1
     jdim=2
     xxmin=xmin ; xxmax=xmax
     yymin=ymin ; yymax=ymax
  end if

  !-----------------------------
  ! Allocate map hierarchy
  !-----------------------------
  do ilevel=levelmin,lmax
     nx_full=2**ilevel
     ny_full=2**ilevel
     imin=int(xxmin*dble(nx_full))+1
     imax=int(xxmax*dble(nx_full))+1
     jmin=int(yymin*dble(ny_full))+1
     jmax=int(yymax*dble(ny_full))+1
     allocate(mapgrid(ilevel)%map(imin:imax,jmin:jmax))
     allocate(mapgrid(ilevel)%rho(imin:imax,jmin:jmax))
     mapgrid(ilevel)%map(:,:)=0.0
     mapgrid(ilevel)%rho(:,:)=0.0
     mapgrid(ilevel)%imin=imin
     mapgrid(ilevel)%imax=imax
     mapgrid(ilevel)%jmin=jmin
     mapgrid(ilevel)%jmax=jmax
  end do

  !-----------------------------------------------
  ! Compute projected variables
  !----------------------------------------------
  allocate(cpu_list(1:ncpu))

  ! Loop over levels 
  do ilevel=levelmin,lmax

     dx=0.5**ilevel
     zmax=1.0/dx
     dxline=1
     if(ndim==3)dxline=dx

     write(*,*)'ilevel=',ilevel,dx,zmax

     ! Get cpu list within bounding box
     call cmp_cpu_list(xmin,xmax,ymin,ymax,zmin,zmax,ilevel,ncpu_read,cpu_list)
     
     ! Loop over processor files
     do k=1,ncpu_read
        icpu=cpu_list(k)
        call title(icpu,ncharcpu)

        ! Prepare reading the AMR file
        file_amr=TRIM(repository)//'/amr.out'//TRIM(ncharcpu)
        noct_skip=0
        open(unit=10,file=file_amr,access="stream",action="read",form='unformatted')
        do i=levelmin,ilevel-1
           read(10,POS=13+4*(i-levelmin))noct_tmp
           noct_skip=noct_skip+noct_tmp
        end do
        read(10,POS=13+4*(ilevel-levelmin))noct_file
        iskip_amr=13+4*(nlevelmax-levelmin+1)+(4*ndim+4*twotondim)*noct_skip

        ! Prepare reading the HYDRO file
        file_hydro=TRIM(repository)//'/hydro.out'//TRIM(ncharcpu)
        open(unit=11,file=file_hydro,access="stream",action="read",form='unformatted')
        iskip_hydro=17+4*(nlevelmax-levelmin+1)+(8*twotondim*nvar)*noct_skip

        ! Loop over useful octs in file
        do i=1,noct_file

           ! Read values from files
           read(10,POS=iskip_amr+(4*ndim+4*twotondim)*(i-1))ckey
           read(10,POS=iskip_amr+(4*ndim+4*twotondim)*(i-1)+4*ndim)refined
           read(11,POS=iskip_hydro+(8*twotondim*nvar)*(i-1))uold

           ! Loop over 2**ndim cells
           do ind=1,twotondim

              ! Compute Cartesian key
              do ldim=1,ndim
                 nstride=2**(ldim-1)
                 cart_key(ldim)=2*ckey(ldim)+MOD((ind-1)/nstride,2)
              end do

              ok_cell=.not.refined(ind).OR.ilevel.eq.lmax

              ! Get variable to map out
              rho=uold(ind,1)
              select case (type)
              case (-1) ! Cpu map
                 map = icpu
              case (0) ! Refinement map
                 map = ilevel
              case (1) ! Density map
                 if(do_max)then
                    map = uold(ind,1)
                 else
                    map = uold(ind,1)**2
                 endif
              case (2) ! Mass weighted x-velocity
                 map = uold(ind,2)
              case (3) ! Mass weighted y-velocity
                 map = uold(ind,3)
              case (4) ! Mass weighted z-velocity
                 map = uold(ind,4)
              case (5) ! Temperature
                 if(do_max)then
                    map = uold(ind,5)/uold(ind,1)
                 else
                    map = uold(ind,5)
                 endif
              case default ! Passive scalar
                 if(do_max)then
                    map = uold(ind,type)/uold(ind,1)
                 else
                    map = uold(ind,type)
                 endif
                 metmax=max(metmax,uold(ind,type))
              end select

              if(ok_cell)then
                 ix=cart_key(idim)+1
                 iy=cart_key(jdim)+1
                 zz=(cart_key(kdim)+0.5)*dx

                 if(ndim==3)then
                    weight=(min(zz+dx/2.,zzmax)-max(zz-dx/2.,zzmin))/dx
                    weight=min(1.0d0,max(weight,0.0d0))
                 else
                    weight=1.0
                 endif

                 if(    ix>=mapgrid(ilevel)%imin.and.&
                      & iy>=mapgrid(ilevel)%jmin.and.&
                      & ix<=mapgrid(ilevel)%imax.and.&
                      & iy<=mapgrid(ilevel)%jmax)then

                    if(do_max)then
                       mapgrid(ilevel)%map(ix,iy)=max(mapgrid(ilevel)%map(ix,iy),map)
                       mapgrid(ilevel)%rho(ix,iy)=max(mapgrid(ilevel)%rho(ix,iy),rho)
                    else
                       mapgrid(ilevel)%map(ix,iy)=mapgrid(ilevel)%map(ix,iy)+map*dxline*weight/(zzmax-zzmin)
                       mapgrid(ilevel)%rho(ix,iy)=mapgrid(ilevel)%rho(ix,iy)+rho*dxline*weight/(zzmax-zzmin)
                    endif
                 endif
              end if

           end do
           ! End loop over cells

        end do
        ! End loop over octs

        ! Close file
        close(10)
        close(11)

     end do
     ! End loop over cpu
     
  end do
  ! End loop over levels

  write(*,*)'Data read and projected.'
  if(type==6.or.type==7)then
     write(*,*)metmax
  endif
  
  nx_full=2**lmax
  ny_full=2**lmax
  imin=int(xxmin*dble(nx_full))+1
  imax=int(xxmax*dble(nx_full))
  jmin=int(yymin*dble(ny_full))+1
  jmax=int(yymax*dble(ny_full))

  do ix=imin,imax
     xmin=((ix-0.5)/2**lmax)
     do iy=jmin,jmax
        ymin=((iy-0.5)/2**lmax)
        do ilevel=levelmin,lmax-1
           ndom=2**ilevel
           i=int(xmin*ndom)+1
           j=int(ymin*ndom)+1
           if(do_max) then
              mapgrid(lmax)%map(ix,iy)=max(mapgrid(lmax)%map(ix,iy), &
                   & mapgrid(ilevel)%map(i,j))
              mapgrid(lmax)%rho(ix,iy)=max(mapgrid(lmax)%rho(ix,iy), &
                   & mapgrid(ilevel)%rho(i,j))
           else
              mapgrid(lmax)%map(ix,iy)=mapgrid(lmax)%map(ix,iy) + &
                   & mapgrid(ilevel)%map(i,j)
              mapgrid(lmax)%rho(ix,iy)=mapgrid(lmax)%rho(ix,iy) + &
                   & mapgrid(ilevel)%rho(i,j)
           endif
        end do
     end do
  end do

  write(*,*)'Norm=',sum(mapgrid(lmax)%rho(imin:imax,jmin:jmax))/(imax-imin+1)/(jmax-jmin+1)

  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Writing map in '//TRIM(nomfich)

  if (mapfiletype=='bin')then
     write(*,*)'unformatted binary file'
     open(unit=20,file=nomfich,form='unformatted')
     if(nx_sample==0)then
        write(20)t, xxmax-xxmin, yymax-yymin, zzmax-zzmin
        write(20)imax-imin+1,jmax-jmin+1
        allocate(toto(imax-imin+1,jmax-jmin+1))
        if(do_max)then
           toto=mapgrid(lmax)%map(imin:imax,jmin:jmax)
        else
           toto=mapgrid(lmax)%map(imin:imax,jmin:jmax)/mapgrid(lmax)%rho(imin:imax,jmin:jmax)
        endif
        write(20)toto
     else
        if(ny_sample==0)ny_sample=nx_sample
        write(20)t, xxmax-xxmin, yymax-yymin, zzmax-zzmin
        write(20)nx_sample+1,ny_sample+1
        allocate(toto(0:nx_sample,0:ny_sample))
        do i=0,nx_sample
           ix=int(dble(i)/dble(nx_sample)*dble(imax-imin+1))+imin
           ix=min(ix,imax)
           do j=0,ny_sample
              iy=int(dble(j)/dble(ny_sample)*dble(jmax-jmin+1))+jmin
              iy=min(iy,jmax)
              if(do_max)then
                 toto(i,j)=mapgrid(lmax)%map(ix,iy)
              else
                 toto(i,j)=mapgrid(lmax)%map(ix,iy)/mapgrid(lmax)%rho(ix,iy)
              endif
           end do
        end do
        write(20)toto
     endif
     close(20)
  endif
  if (mapfiletype=='ascii')then
     write(*,*)'formatted text file'
     open(unit=20,file=nomfich,form='formatted')
     if(nx_sample==0)then
        do j=jmin,jmax
           do i=imin,imax
              xx=xxmin+(dble(i-imin)+0.5)/dble(imax-imin+1)*(xxmax-xxmin)
              yy=yymin+(dble(j-jmin)+0.5)/dble(jmax-jmin+1)*(yymax-yymin)
              if(do_max)then
                 write(20,*)xx,yy,mapgrid(lmax)%map(i,j)
              else
                 write(20,*)xx,yy,mapgrid(lmax)%map(i,j)/mapgrid(lmax)%rho(i,j)
              endif
           end do
           write(20,*) " "
        end do
     else  
        if(ny_sample==0)ny_sample=nx_sample
        allocate(toto(0:nx_sample,0:ny_sample))
        do i=0,nx_sample
           ix=int(dble(i)/dble(nx_sample)*dble(imax-imin+1))+imin
           ix=min(ix,imax)
           do j=0,ny_sample
              iy=int(dble(j)/dble(ny_sample)*dble(jmax-jmin+1))+jmin
              iy=min(iy,jmax)
              if(do_max)then
                 toto(i,j)=mapgrid(lmax)%map(ix,iy)
              else
                 toto(i,j)=mapgrid(lmax)%map(ix,iy)/mapgrid(lmax)%rho(ix,iy)
           endif
           end do
        end do
        do j=0,ny_sample
           do i=0,nx_sample
              xx=xxmin+dble(i)/dble(nx_sample)*(xxmax-xxmin)
              yy=yymin+dble(j)/dble(ny_sample)*(yymax-yymin)
              write(20,*)xx,yy,toto(i,j)
           end do
           write(20,*) " "
        end do
     endif
     close(20)
  endif

contains
  
  subroutine read_params
    
    implicit none
    
    integer       :: i,n
    integer       :: iargc
    character(len=4)   :: opt
    character(len=128) :: arg
    LOGICAL       :: bad, ok
    
    n = iargc()
    if (n < 4) then
       print *, 'usage: amr2map   -inp  input_dir'
       print *, '                 -out  output_file'
       print *, '                 [-dir axis] '
       print *, '                 [-xmi xmin] '
       print *, '                 [-xma xmax] '
       print *, '                 [-ymi ymin] '
       print *, '                 [-yma ymax] '
       print *, '                 [-zmi zmin] '
       print *, '                 [-zma zmax] '
       print *, '                 [-lma lmax] '
       print *, '                 [-typ type] '
       print *, '                 [-fil filetype] '
       print *, '                 [-max maxi] '
       print *, 'ex: amr2map -inp output_00001 -out map.dat'// &
            &   ' -dir z -xmi 0.1 -xma 0.7 -lma 12'
       print *, ' '
       print *, ' type :-1 = cpu number'
       print *, ' type : 0 = ref. level (default)'
       print *, ' type : 1 = gas density'
       print *, ' type : 2 = X velocity'
       print *, ' type : 3 = Y velocity'
       print *, ' type : 4 = Z velocity'
       print *, ' type : 5 = gas pressure'
       print *, ' type : 6 = gas metallicity'
       print *, ' '
       print *, ' maxi : 0 = average along line of sight (default)'
       print *, ' maxi : 1 = maximum along line of sight'
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
       case ('-lma')
          read (arg,*) lmax
       case ('-nx')
          read (arg,*) nx_sample
       case ('-ny')
          read (arg,*) ny_sample
       case ('-typ')
          read (arg,*) type
       case ('-fil')
          read (arg,*) mapfiletype
       case ('-max')
          read (arg,*) domax
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
    
    do_max=.false.
    if(domax==1)do_max=.true.

    return
    
  end subroutine read_params
  
end program amr2map

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
!================================================================
!================================================================
!================================================================
!================================================================
function check_ramses_exist(repository)
  implicit none
  character(LEN=128)::repository
  logical::check_ramses_exist
  !-----------------------------------------------
  ! Check that RAMSES files are there
  !-----------------------------------------------
  integer::ipos
  character(LEN=5)::char
  character(LEN=128)::nomfich
  logical::ok

  check_ramses_exist=.true.
  ipos=INDEX(repository,'output_')
  nomfich=TRIM(repository)//'/hydro.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     check_ramses_exist=.false.
  endif
  nomfich=TRIM(repository)//'/amr.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     check_ramses_exist=.false.
  endif
  nomfich=TRIM(repository)//'/params.out'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     check_ramses_exist=.false.
  endif

end function check_ramses_exist
!================================================================
!================================================================
!================================================================
!================================================================
subroutine read_ramses_params(repository)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  character(LEN=128)::repository
  !-----------------------------------------------
  ! Read RAMSES parameters file
  !-----------------------------------------------
  character(LEN=128)::nomfich
  integer::ilun,ilevel,ndim2,nhilbert2

  nomfich=TRIM(repository)//'/params.out'

  ilun=10
  open(unit=ilun,file=nomfich,access="stream",action="read",form='unformatted')
  ! Read grid variables
  read(ilun)ncpu
  read(ilun)ndim2
  if(ndim.NE.ndim2)then
     write(*,*)'Recompile'
     write(*,*)'ndim=',ndim
     write(*,*)'ndim in file=',ndim2
     stop
  endif
  read(ilun)levelmin
  read(ilun)nlevelmax
  ! Overwrite boxlen with value from file
  read(ilun)boxlen
  ! Read time variables
  read(ilun)noutput,iout,ifout
  read(ilun)tout(1:noutput)
  read(ilun)aout(1:noutput)
  read(ilun)t
  read(ilun)dtold(1:nlevelmax)
  read(ilun)dtnew(1:nlevelmax)
  read(ilun)nstep,nstep_coarse
  ! Read various constants
  read(ilun)const,mass_tot_0,rho_tot
  read(ilun)omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini
  read(ilun)aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
  read(ilun)mass_sph
  read(ilun)nhilbert2
  if(nhilbert.NE.nhilbert2)then
     write(*,*)'Recompile'
     write(*,*)'nhilbert=',nhilbert
     write(*,*)'nhilbert in file=',nhilbert2
     stop
  endif
  allocate(bound_key_level(1:nhilbert,0:ncpu,1:nlevelmax+1))
  do ilevel=levelmin,nlevelmax
     read(ilun)bound_key_level(1:nhilbert,0:ncpu,ilevel)
  end do
  close(ilun)

end subroutine read_ramses_params
!================================================================
!================================================================
!================================================================
!================================================================
subroutine cmp_cpu_list(xmin,xmax,ymin,ymax,zmin,zmax,ilevel,ncpu_read,cpu_list)
  use amr_commons
  use hilbert
  implicit none
  real(dp)::xmin,xmax,ymin,ymax,zmin,zmax
  integer::ilevel,ncpu_read
  integer,dimension(1:ncpu)::cpu_list
  !----------------------------------------------
  ! Set up Hilbert key for optimal file reading
  !----------------------------------------------
  logical,dimension(1:ncpu)::cpu_read
  integer(kind=4),dimension(1:nvector),save::dummy_state
  integer(kind=8),dimension(1:nvector,1:ndim)::ix
  integer(kind=8),dimension(1:nvector,1:nhilbert)::hk
  integer(kind=8),dimension(1:nhilbert)::one_key,order_min,order_max
  integer(kind=8),dimension(1:nhilbert,1:8)::bounding_min,bounding_max
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(dp)::dmax,dx,maxdom
  integer::i,j,ilev,lmin,imin,imax,jmin,jmax,kmin,kmax,bit_length,ndom,impi
  logical::ok1,ok2

  ncpu_read=ncpu
  do i=1,ncpu
     cpu_list(i)=i
  end do
  return

  dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
  do ilev=1,nlevelmax
     dx=0.5d0**ilev
     if(dx.lt.dmax)exit
  end do
  lmin=MIN(ilev,ilevel)
  bit_length=lmin-1
  maxdom=2**bit_length
  imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
  if(bit_length>0)then
     imin=int(xmin*dble(maxdom))
     imax=imin+1
     jmin=int(ymin*dble(maxdom))
     jmax=jmin+1
     kmin=int(zmin*dble(maxdom))
     kmax=kmin+1
  endif
  
  ndom=1
  if(bit_length>0)ndom=8
  idom(1)=imin; idom(2)=imax
  idom(3)=imin; idom(4)=imax
  idom(5)=imin; idom(6)=imax
  idom(7)=imin; idom(8)=imax
  jdom(1)=jmin; jdom(2)=jmin
  jdom(3)=jmax; jdom(4)=jmax
  jdom(5)=jmin; jdom(6)=jmin
  jdom(7)=jmax; jdom(8)=jmax
  kdom(1)=kmin; kdom(2)=kmin
  kdom(3)=kmin; kdom(4)=kmin
  kdom(5)=kmax; kdom(6)=kmax
  kdom(7)=kmax; kdom(8)=kmax
  
  one_key=0
  one_key(1)=1
  do i=1,ndom
     if(bit_length>0)then
        ix(1,1)=idom(i)
        if(ndim>1)ix(1,2)=jdom(i)
        if(ndim>2)ix(1,3)=kdom(i)
        call hilbert_key(ix,hk,dummy_state,0,bit_length,1)
        order_min(1:nhilbert)=hk(1,1:nhilbert)
     else
        order_min(1:nhilbert)=0
     endif
     bounding_min(1:nhilbert,i)=order_min
     bounding_max(1:nhilbert,i)=order_min+one_key
  end do

  if(ilevel>bit_length+1)then
     do ilev=bit_length+1,ilevel
        do i=1,ndom
           bounding_min(1:nhilbert,i)=refine_key(bounding_min(1:nhilbert,i),ilevel)
        end do
     end do
  endif

  cpu_min=0; cpu_max=0
  do impi=1,ncpu
     do i=1,ndom
        ok1=ge_keys(bounding_min(1:nhilbert,i),bound_key_level(1:nhilbert,impi-1,ilevel))
        ok2=gt_keys(bound_key_level(1:nhilbert,impi,ilevel),bounding_min(1:nhilbert,i))
        if(ok1.and.ok2)then
           cpu_min(i)=impi
        endif
        ok1=gt_keys(bounding_max(1:nhilbert,i),bound_key_level(1:nhilbert,impi-1,ilevel))
        ok2=ge_keys(bound_key_level(1:nhilbert,impi,ilevel),bounding_max(1:nhilbert,i))
        if(ok1.and.ok2)then
           cpu_max(i)=impi
        endif
     end do
  end do
  
  ncpu_read=0
  do i=1,ndom
     do j=cpu_min(i),cpu_max(i)
        if(.not. cpu_read(j))then
           ncpu_read=ncpu_read+1
           cpu_list(ncpu_read)=j
           cpu_read(j)=.true.
        endif
     enddo
  enddo

end subroutine cmp_cpu_list
!================================================================
!================================================================
!================================================================
!================================================================
