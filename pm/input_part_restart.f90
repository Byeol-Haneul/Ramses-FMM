!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part_restart(pst)
  use amr_parameters, only: ndim,dp,i8b
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from a Ramses restart file.
  !--------------------------------------------------------------------
  integer::icpu,ilun,ncpu_file
  integer(i8b)::npart_tot_file,npart_tot_check
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::file_head,file_part
  integer,allocatable,dimension(:)::npart_file
  
  if(pst%s%r%verbose)write(*,*)'Entering input_part_restart'

  ! Read particle files header
  call title(pst%s%r%nrestart,nchar)
  file_head='output_'//TRIM(nchar)//'/part_header.txt'
  call input_header(pst%s%r,pst%s%g,file_head,npart_tot_file,ncpu_file)
  write(*,'(" Restart snapshot has ",I8," particles")')npart_tot_file

  ! Allocate local array
  allocate(npart_file(1:ncpu_file))

  ! Read number of particles in each file
  npart_tot_check=0
  do icpu=1,ncpu_file
     call title(icpu,ncharcpu)
     file_part='output_'//TRIM(nchar)//'/part.out'//TRIM(ncharcpu)
     ilun=10
     open(unit=ilun,file=TRIM(file_part),access="stream",action="read",form='unformatted')
     read(ilun,POS=5)npart_file(icpu)
     npart_tot_check=npart_tot_check+npart_file(icpu)
     close(ilun)
  end do
  if(npart_tot_check.NE.npart_tot_file)then
     write(*,*)' Input file corrupted'
     call mdl_abort
  endif

  ! Call recursive slave routine
  call r_input_part_restart(pst,pst%s%mdl%ncpu,ncpu_file,0,npart_file)

  ! Deallocate local array
  deallocate(npart_file)

end subroutine m_input_part_restart
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_input_part_restart(pst,cpu_range,input_size,output_size,input_array)
  use amr_parameters, only: dp
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::cpu_range,input_size,output_size
  integer,dimension(1:input_size)::input_array
  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to read and dispatch
  ! particles from a Ramses restart file.
  !--------------------------------------------------------------------
  integer::next_range,next_cpu
  
  next_range=cpu_range/2
  next_cpu=pst%s%g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(pst%s%mdl,MDL_INPUT_PART_RESTART,next_cpu,next_range,input_size,output_size,input_array)
     call r_input_part_restart(pst,next_range,input_size,output_size,input_array)
  else
     call input_part_restart(pst%s%r,pst%s%g,pst%s%p,input_size,input_array)
  endif

end subroutine r_input_part_restart
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_part_restart(r,g,p,ncpu_file,npart_file)
  use amr_parameters, only: ndim,dp,i8b
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  integer::ncpu_file
  integer,dimension(1:ncpu_file)::npart_file
  !------------------------------------------------------------
  ! Read particles positions and velocities from a Ramses 
  ! restart file and allocate particle-based arrays.
  !------------------------------------------------------------
  integer::ipart,ipart_old
  integer::i,idim,icpu,ileft,iright,nrest,ipos
  integer(i8b)::istart,iend,nleft,nright,npart_tot
  integer(i8b),dimension(0:ncpu_file)::ncum_file

  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::isp
  integer(i8b),allocatable,dimension(:)::isp8

  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::file_part
  
  !-------------------------------------
  ! Compute local particle number
  !-------------------------------------
  ncum_file(0)=0
  do icpu=1,ncpu_file
     ncum_file(icpu)=ncum_file(icpu-1)+npart_file(icpu)
  end do
  npart_tot=ncum_file(ncpu_file)
  p%npart_tot=npart_tot
  
  p%npart=npart_tot/g%ncpu
  nleft=(g%myid-1)*p%npart
  nright=g%myid*p%npart
  nrest=npart_tot-p%npart*g%ncpu
  if(g%myid.LE.nrest)then
     p%npart=p%npart+1
     nleft=nleft+(g%myid-1)
     nright=nright+g%myid
  endif
  
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
  call title(r%nrestart,nchar)
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
     open(unit=10,file=TRIM(file_part),access="stream",action="read",form='unformatted')
     
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

  ! Put all particles in levelmin
  p%headp=p%npart+1
  p%tailp=p%npart
  p%headp(r%levelmin)=1
  p%tailp(r%levelmin)=p%npart
        
end subroutine input_part_restart
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

