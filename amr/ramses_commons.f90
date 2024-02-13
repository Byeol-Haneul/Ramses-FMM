module ramses_commons
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_module, only: mdl_t
  use clfind_commons, only: clump_t
  use cooling_module, only: cooling_t

  type ramses_t

     type(run_t)::r
     type(global_t)::g
     type(mesh_t)::m
     type(part_t)::p
     type(part_t)::s
     type(clump_t)::c
     type(cooling_t)::cool
     type(mdl_t),pointer::mdl => null()

  end type ramses_t

  type pst_t
     
     type(ramses_t),pointer::s => null()
     type(pst_t),pointer::pLower => null()
     integer::iUpper = -1
     integer::nLower = 0
     integer::nUpper = 0
     
  end type pst_t

contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine open_file(s,filename,nskip,ilun)
  use hydro_parameters, only: nvar
  use amr_parameters, only: ndim,flen,twotondim
#ifndef WITHOUTMPI
  use mpi
#endif
  type(ramses_t)::s
  character(len=flen),intent(in)::filename
  integer,dimension(s%r%levelmin:s%r%nlevelmax),intent(out)::nskip
  integer,intent(out)::ilun
#ifndef WITHOUTMPI
  integer::info,reqsend,reqrecv,tag1=100,tag2=200
  integer,dimension(MPI_STATUS_SIZE)::reqstatus
#endif
  character(LEN=flen)::fileloc
  character(LEN=5)::nchar
  integer,dimension(1:s%r%nfile+1)::istart
  integer::i,ifile,ncpufile,nremain,ilevel,ierr,iskip=0
  integer,dimension(s%r%levelmin:s%r%nlevelmax)::noct
  logical::file_exist
  
  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ncpufile=g%ncpu/r%nfile
  nremain=g%ncpu-r%nfile*ncpufile
  istart(1)=1
  do i=1,r%nfile
     if(i.LE.nremain)then
        istart(i+1)=istart(i)+ncpufile+1
     else
        istart(i+1)=istart(i)+ncpufile
     endif
  end do

  ifile=1
  do i=1,r%nfile
     if(g%myid.GE.istart(i).AND.g%myid.LT.istart(i+1))then
        ifile=i
        exit
     end if
  end do

#ifndef WITHOUTMPI
  if(g%myid.EQ.istart(ifile+1)-1)then
     noct=m%noct
     if(g%myid.GT.istart(ifile))then
     call MPI_ISEND(noct,r%nlevelmax-r%levelmin+1,MPI_INTEGER,g%myid-2,tag1,MPI_COMM_WORLD,reqsend,info)
     call MPI_WAIT(reqsend,reqstatus,info)
     endif
  elseif(g%myid.GT.istart(ifile))then
     call MPI_IRECV(noct,r%nlevelmax-r%levelmin+1,MPI_INTEGER,g%myid,tag1,MPI_COMM_WORLD,reqrecv,info)
     call MPI_WAIT(reqrecv,reqstatus,info)
     noct=noct+m%noct
     call MPI_ISEND(noct,r%nlevelmax-r%levelmin+1,MPI_INTEGER,g%myid-2,tag1,MPI_COMM_WORLD,reqsend,info)
     call MPI_WAIT(reqsend,reqstatus,info)
  elseif(g%myid.EQ.istart(ifile))then
     call MPI_IRECV(noct,r%nlevelmax-r%levelmin+1,MPI_INTEGER,g%myid,tag1,MPI_COMM_WORLD,reqrecv,info)
     call MPI_WAIT(reqrecv,reqstatus,info)
     noct=noct+m%noct
  else
     write(*,*)'Problem in open_file'
  endif
#else
  noct=m%noct
#endif

  if(g%myid.EQ.istart(ifile))then

     ilun=10+istart(ifile)
     call title(ifile,nchar)
     fileloc=TRIM(filename)//TRIM(nchar)
     inquire(file=fileloc, exist=file_exist)
     if (file_exist) then
        open(unit=ilun,file=fileloc,iostat=ierr)
        close(ilun,status="delete")
     end if

     if(index(filename,'clump').NE.0)then
        open(unit=ilun,file=fileloc,form='formatted')
        write(ilun,'(144A)')'     index       halo  lev   parent      ncell    peak_x             peak_y             peak_z     '//&
             '        rho-               rho+               rho_av             mass_cl            relevance   '
     elseif(index(filename,'halo').NE.0)then
        open(unit=ilun,file=fileloc,form='formatted')
        write(ilun,'(135A)')'     index      ncell    peak_x             peak_y             peak_z     '//&
             '        rho+               mass      '
     else
        open(unit=ilun,file=fileloc,access="stream",action="write",form='unformatted')
        write(ilun)ndim
        if(index(filename,'hydro').NE.0)write(ilun)nvar
        if(index(filename,'grav').NE.0)write(ilun)ndim+1
        if(index(filename,'peak').NE.0)write(ilun)2
        write(ilun)r%levelmin
        write(ilun)r%nlevelmax
        do ilevel=r%levelmin,r%nlevelmax
           write(ilun)noct(ilevel)
        end do
     endif

     if(index(filename,'amr').NE.0)iskip=13+(r%nlevelmax-r%levelmin+1)*4
     if(index(filename,'hydro').NE.0)iskip=17+(r%nlevelmax-r%levelmin+1)*4
     if(index(filename,'grav').NE.0)iskip=17+(r%nlevelmax-r%levelmin+1)*4
     if(index(filename,'peak').NE.0)iskip=17+(r%nlevelmax-r%levelmin+1)*4

     do ilevel=r%levelmin,r%nlevelmax
        nskip(ilevel)=iskip
        if(index(filename,'amr').NE.0)iskip=iskip+(4*ndim+4*twotondim)*noct(ilevel)
        if(index(filename,'hydro').NE.0)iskip=iskip+(4*twotondim*nvar)*noct(ilevel)
        if(index(filename,'grav').NE.0)iskip=iskip+(4*twotondim*(ndim+1))*noct(ilevel)
        if(index(filename,'peak').NE.0)iskip=iskip+(4*twotondim*2)*noct(ilevel)
     end do

  elseif(g%myid.GT.istart(ifile))then

#ifndef WITHOUTMPI
     call MPI_IRECV(nskip,r%nlevelmax-r%levelmin+1,MPI_INTEGER,g%myid-2,tag2,MPI_COMM_WORLD,reqrecv,info)
     call MPI_WAIT(reqrecv,reqstatus,info)
#endif
     ilun=10+istart(ifile)
     call title(ifile,nchar)
     fileloc=TRIM(filename)//TRIM(nchar)
     inquire(file=fileloc, exist=file_exist)
     if (.not. file_exist) then
        write(*,*)'Problem in open_file'
        stop
     end if

     if(index(filename,'clump').NE.0)then
        open(unit=ilun,file=fileloc,form='formatted',access='append')
     elseif(index(filename,'halo').NE.0)then
        open(unit=ilun,file=fileloc,form='formatted',access='append')
     else
        open(unit=ilun,file=fileloc,access="stream",action="write",form='unformatted')
     endif

  end if

  end associate

end subroutine open_file
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine close_file(s,filename,nskip,ilun)
  use hydro_parameters, only: nvar
  use amr_parameters, only: ndim,flen,twotondim
#ifndef WITHOUTMPI
  use mpi
#endif
  type(ramses_t)::s
  character(len=flen),intent(in)::filename
  integer,dimension(s%r%levelmin:s%r%nlevelmax),intent(inout)::nskip
  integer,intent(in)::ilun
#ifndef WITHOUTMPI
  integer::info,reqsend,tag2=200
  integer,dimension(MPI_STATUS_SIZE)::reqstatus
#endif
  integer,dimension(1:s%r%nfile+1)::istart
  integer::ncpufile,nremain
  integer::i,iskip,ifile,ilevel

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ncpufile=g%ncpu/r%nfile
  nremain=g%ncpu-r%nfile*ncpufile
  istart(1)=1
  do i=1,r%nfile
     if(i.LE.nremain)then
        istart(i+1)=istart(i)+ncpufile+1
     else
        istart(i+1)=istart(i)+ncpufile
     endif
  end do

  ifile=1
  do i=1,r%nfile
     if(g%myid.GE.istart(i).AND.g%myid.LT.istart(i+1))then
        ifile=i
        exit
     end if
  end do

  close(ilun)

  if(g%myid.LT.istart(ifile+1)-1)then

     do ilevel=r%levelmin,r%nlevelmax
        iskip=nskip(ilevel)
        if(index(filename,'amr').NE.0)iskip=iskip+(4*ndim+4*twotondim)*m%noct(ilevel)
        if(index(filename,'hydro').NE.0)iskip=iskip+(4*twotondim*nvar)*m%noct(ilevel)
        if(index(filename,'grav').NE.0)iskip=iskip+(4*twotondim*(ndim+1))*m%noct(ilevel)
        if(index(filename,'peak').NE.0)iskip=iskip+(4*twotondim*2)*m%noct(ilevel)
        nskip(ilevel)=iskip
     end do

#ifndef WITHOUTMPI
     call MPI_ISEND(nskip,r%nlevelmax-r%levelmin+1,MPI_INTEGER,g%myid,tag2,MPI_COMM_WORLD,reqsend,info)
     call MPI_WAIT(reqsend,reqstatus,info)
#endif
  end if

  end associate

end subroutine close_file
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine open_part_file(s,p,filename,nskip,ilun)
  use amr_parameters, only: ndim,flen,twotondim,i8b
  use pm_commons, only: part_t
#ifndef WITHOUTMPI
  use mpi
#endif
  type(ramses_t)::s
  type(part_t)::p
  character(len=flen),intent(in)::filename
  integer, dimension(1:p%nvaralloc+1),intent(out)::nskip
  integer, intent(out)::ilun
#ifndef WITHOUTMPI
  integer::info,reqsend,reqrecv,tag1=100,tag2=200
  integer,dimension(MPI_STATUS_SIZE)::reqstatus
#endif
  character(LEN=flen)::fileloc
  character(LEN=5)::nchar
  integer,dimension(1:s%r%nfile+1)::istart
  integer::i,ivar,ifile,ncpufile,nremain,ilevel,ierr
  integer::npart
  logical::file_exist

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  ncpufile=g%ncpu/r%nfile
  nremain=g%ncpu-r%nfile*ncpufile
  istart(1)=1
  do i=1,r%nfile
     if(i.LE.nremain)then
        istart(i+1)=istart(i)+ncpufile+1
     else
        istart(i+1)=istart(i)+ncpufile
     endif
  end do

  ifile=1
  do i=1,r%nfile
     if(g%myid.GE.istart(i).AND.g%myid.LT.istart(i+1))then
        ifile=i
        exit
     end if
  end do

#ifndef WITHOUTMPI
  if(g%myid.EQ.istart(ifile+1)-1)then
     npart=p%npart
     if(g%myid.GT.istart(ifile))then
     call MPI_ISEND(npart,1,MPI_INTEGER,g%myid-2,tag1,MPI_COMM_WORLD,reqsend,info)
     call MPI_WAIT(reqsend,reqstatus,info)
     endif
  elseif(g%myid.GT.istart(ifile))then
     call MPI_IRECV(npart,1,MPI_INTEGER,g%myid,tag1,MPI_COMM_WORLD,reqrecv,info)
     call MPI_WAIT(reqrecv,reqstatus,info)
     npart=npart+p%npart
     call MPI_ISEND(npart,1,MPI_INTEGER,g%myid-2,tag1,MPI_COMM_WORLD,reqsend,info)
     call MPI_WAIT(reqsend,reqstatus,info)
  elseif(g%myid.EQ.istart(ifile))then
     call MPI_IRECV(npart,1,MPI_INTEGER,g%myid,tag1,MPI_COMM_WORLD,reqrecv,info)
     call MPI_WAIT(reqrecv,reqstatus,info)
     npart=npart+p%npart
  else
     write(*,*)'Problem in open_part_file'
  endif
#else
  npart=p%npart
#endif

  if(g%myid.EQ.istart(ifile))then

     ilun=10+istart(ifile)
     call title(ifile,nchar)
     fileloc=TRIM(filename)//TRIM(nchar)
     inquire(file=fileloc, exist=file_exist)
     if (file_exist) then
        open(unit=ilun,file=fileloc,iostat=ierr)
        close(ilun,status="delete")
     end if
     open(unit=ilun,file=fileloc,access="stream",action="write",form='unformatted')

     write(ilun)ndim
     write(ilun)npart
     nskip(1)=9

     ! Positions
     do ivar=1,ndim
        nskip(ivar+1)=nskip(ivar)+4*npart
     end do
     ! Velocities
     do ivar=ndim+1,2*ndim
        nskip(ivar+1)=nskip(ivar)+4*npart
     end do
     ! Masses
     ivar=2*ndim+1
     nskip(ivar+1)=nskip(ivar)+4*npart
     ! Metallicities
     if(allocated(p%zp))then
        ivar=ivar+1
        nskip(ivar+1)=nskip(ivar)+4*npart
     endif
     ! Birth times
     if(allocated(p%tp))then
        ivar=ivar+1
        nskip(ivar+1)=nskip(ivar)+4*npart
     endif
     ! Identities
     ivar=ivar+1
     nskip(ivar+1)=nskip(ivar)+i8b*npart
     ! Levels
     ivar=ivar+1
     nskip(ivar+1)=nskip(ivar)+4*npart
#ifdef OUTPUT_PARTICLE_POTENTIAL
     ! Potentials
     ivar=ivar+1
     nskip(ivar+1)=nskip(ivar)+4*npart
#endif

  elseif(g%myid.GT.istart(ifile))then

#ifndef WITHOUTMPI
     call MPI_IRECV(nskip,p%nvaralloc+1,MPI_INTEGER,g%myid-2,tag2,MPI_COMM_WORLD,reqrecv,info)
     call MPI_WAIT(reqrecv,reqstatus,info)
#endif
     ilun=10+istart(ifile)
     call title(ifile,nchar)
     fileloc=TRIM(filename)//TRIM(nchar)
     inquire(file=fileloc, exist=file_exist)
     if (.not. file_exist) then
        write(*,*)'Problem in open_part_file'
        stop
     end if
     open(unit=ilun,file=fileloc,access="stream",action="write",form='unformatted')

  end if

  end associate

end subroutine open_part_file
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine close_part_file(s,p,filename,nskip,ilun)
  use amr_parameters, only: ndim,flen,twotondim,i8b
  use pm_commons, only: part_t
#ifndef WITHOUTMPI
  use mpi
#endif
  type(ramses_t)::s
  type(part_t)::p
  character(len=flen),intent(in)::filename
  integer, dimension(1:p%nvaralloc+1),intent(inout)::nskip
  integer, intent(in)::ilun
#ifndef WITHOUTMPI
  integer::info,reqsend,tag2=200
  integer,dimension(MPI_STATUS_SIZE)::reqstatus
#endif
  integer,dimension(1:s%r%nfile+1)::istart
  integer::ncpufile,nremain
  integer::i,ivar,ifile,ilevel

  associate(r=>s%r,g=>s%g)

  ncpufile=g%ncpu/r%nfile
  nremain=g%ncpu-r%nfile*ncpufile
  istart(1)=1
  do i=1,r%nfile
     if(i.LE.nremain)then
        istart(i+1)=istart(i)+ncpufile+1
     else
        istart(i+1)=istart(i)+ncpufile
     endif
  end do

  ifile=1
  do i=1,r%nfile
     if(g%myid.GE.istart(i).AND.g%myid.LT.istart(i+1))then
        ifile=i
        exit
     end if
  end do

  close(ilun)

  if(g%myid.LT.istart(ifile+1)-1)then

     ! Positions
     do ivar=1,ndim
        nskip(ivar)=nskip(ivar)+4*p%npart
     end do
     ! Velocities
     do ivar=ndim+1,2*ndim
        nskip(ivar)=nskip(ivar)+4*p%npart
     end do
     ! Masses
     ivar=2*ndim+1
     nskip(ivar)=nskip(ivar)+4*p%npart
     ! Metallicities
     if(allocated(p%zp))then
        ivar=ivar+1
        nskip(ivar)=nskip(ivar)+4*p%npart
     endif
     ! Birth times
     if(allocated(p%tp))then
        ivar=ivar+1
        nskip(ivar)=nskip(ivar)+4*p%npart
     endif
     ! Identities
     ivar=ivar+1
     nskip(ivar)=nskip(ivar)+i8b*p%npart
     ! Levels
     ivar=ivar+1
     nskip(ivar)=nskip(ivar)+4*p%npart
#ifdef OUTPUT_PARTICLE_POTENTIAL
     ! Potentials
     ivar=ivar+1
     nskip(ivar)=nskip(ivar)+4*p%npart
#endif

#ifndef WITHOUTMPI
     call MPI_ISEND(nskip,p%nvaralloc+1,MPI_INTEGER,g%myid,tag2,MPI_COMM_WORLD,reqsend,info)
     call MPI_WAIT(reqsend,reqstatus,info)
#endif
  end if

  end associate

end subroutine close_part_file
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module ramses_commons

