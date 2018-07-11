!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part_grafic(pst)
  use amr_parameters, only: dp,i8b
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from a Ramses restart file.
  !--------------------------------------------------------------------
  integer::dummy
  integer,allocatable,dimension(:)::input_array

  if(pst%s%r%verbose)write(*,*)'Entering input_part_grafic'
  if(TRIM(pst%s%r%filetype).NE.'grafic')return
  if(pst%s%r%nrestart>0)return

  ! Compute total number of particles in file
  if(TRIM(pst%s%r%initfile(pst%s%r%levelmin)).NE.' ')then
     pst%s%p%npart_tot=pst%s%g%n1(pst%s%r%levelmin)*pst%s%g%n2(pst%s%r%levelmin)*pst%s%g%n3(pst%s%r%levelmin)
     write(*,*)'Found npart_tot=',pst%s%p%npart_tot
  else
     pst%s%p%npart_tot=0
  endif

  ! If no particle found, no need to read
  if(pst%s%p%npart_tot==0)then
     return
  endif

  ! Call recursive slave routine
  allocate(input_array(1:storage_size(pst%s%p%npart_tot)/32))
  input_array=transfer(pst%s%p%npart_tot,input_array)
  call r_input_part_grafic(pst,input_array,storage_size(pst%s%p%npart_tot)/32,dummy,0)
  deallocate(input_array)

end subroutine m_input_part_grafic
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_input_part_grafic(pst,input_array,input_size,output_array,output_size)
  use amr_parameters, only: i8b
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array
  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to read and dispatch
  ! particles from a Ramses restart file.
  !--------------------------------------------------------------------
  integer(i8b)::npart_tot

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_INPUT_PART_GRAFIC,pst%iUpper+1,input_size,output_size,input_array)
     call r_input_part_grafic(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     npart_tot=transfer(input_array,npart_tot)
     call input_part_grafic(pst%s%r,pst%s%g,pst%s%p,npart_tot)
  endif

end subroutine r_input_part_grafic
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_part_grafic(r,g,p,npart_tot)
  use amr_parameters, only: i8b,dp,ndim,twotondim
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  integer(i8b)::npart_tot
  !----------------------------------------------------
  ! Reading initial conditions from single GRAFIC file
  !----------------------------------------------------
  integer::icpu,ipart,ipart_grafic,idim
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::plane_size
  real(dp)::dx,xx1,xx2,xx3
  real(kind=8)::dispmax=0.0
  integer,dimension(1:g%ncpu)::npart_loc
  integer(i8b),dimension(1:g%ncpu+1)::start_ind
  real(kind=4),dimension(:,:),allocatable::init_plane,init_plane_x
  character(LEN=80)::filename,filename_x
  character(LEN=5)::nchar
  logical::ok,error,keep_part,read_pos=.false.

  !-------------------------------------------
  ! Mesh size at levelmin in normalised units
  !-------------------------------------------
  dx=0.5d0**r%levelmin

  !--------------------------------------
  ! Compute starting index for each cpu
  !--------------------------------------
  p%npart_tot=npart_tot
  do icpu=1,g%ncpu+1
     start_ind(icpu)=((icpu-1)*npart_tot)/g%ncpu
     if(icpu>1)then
        npart_loc(icpu-1)=start_ind(icpu)-start_ind(icpu-1)
     endif
  end do
  p%npart=npart_loc(g%myid)

  !--------------------------------------
  ! Initialize particles in planes
  !--------------------------------------
  plane_size=g%n1(r%levelmin)*g%n2(r%levelmin)
  i3_min=start_ind(g%myid)/plane_size
  i3_max=(start_ind(g%myid+1)-1)/plane_size

  ipart=1
  ipart_grafic=i3_min*plane_size

  ! Initialize positions, masses and ids
  do i3=i3_min,i3_max
     do i2=0,g%n2(r%levelmin)-1
        do i1=0,g%n1(r%levelmin)-1
           keep_part=(ipart_grafic>=start_ind(g%myid).AND.ipart<=p%npart)
           if(keep_part)then
              xx1=(dble(i1)+0.5)/dble(g%n1(r%levelmin))
              xx2=(dble(i2)+0.5)/dble(g%n2(r%levelmin))
              xx3=(dble(i3)+0.5)/dble(g%n3(r%levelmin))
              if(ndim>0)p%xp(ipart,1)=xx1
              if(ndim>1)p%xp(ipart,2)=xx2
              if(ndim>2)p%xp(ipart,3)=xx3
              p%mp(ipart)=0.5d0**(3*r%levelmin)*(1.0d0-g%omega_b/g%omega_m)
              p%idp(ipart)=ipart_grafic+1
              ipart=ipart+1
           endif
           ipart_grafic=ipart_grafic+1
        end do
     end do
  end do
  
  !--------------------------------------
  ! Allocate temporary arrays
  !--------------------------------------
  allocate(init_plane(1:g%n1(r%levelmin),1:g%n2(r%levelmin)))
  filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscx'
  INQUIRE(file=filename_x,exist=ok)
  read_pos=.false.
  if(ok)then
     read_pos=.true.
     allocate(init_plane_x(1:g%n1(r%levelmin),1:g%n2(r%levelmin)))
  else
     if(g%myid==1)write(*,*)'File '//TRIM(filename_x)//' not found.'
  endif
  
  !--------------------------------------
  ! Loop over displacements dimensions
  !--------------------------------------
  do idim=1,ndim
     
     ! Read dark matter Zeldovich initial displacement field
     if(idim==1)filename=TRIM(r%initfile(r%levelmin))//'/ic_velcx'
     if(idim==2)filename=TRIM(r%initfile(r%levelmin))//'/ic_velcy'
     if(idim==3)filename=TRIM(r%initfile(r%levelmin))//'/ic_velcz'
     
     if(g%myid==1)write(*,*)'Reading file '//TRIM(filename)     
     open(10,file=filename,form='unformatted')
     rewind 10
     read(10) ! skip first line
     do i3=0,i3_min-1
        read(10) ! skip unnecessary planes
     end do

     ! If present, read higher order dark matter initial displacement field
     if(read_pos)then
        if(idim==1)filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscx'
        if(idim==2)filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscy'
        if(idim==3)filename_x=TRIM(r%initfile(r%levelmin))//'/ic_poscz'
        
        if(g%myid==1)write(*,*)'Reading file '//TRIM(filename_x)
        open(11,file=filename_x,form='unformatted')
        rewind 11
        read(11) ! skip first line
        do i3=0,i3_min-1
           read(11) ! skip unnecessary planes
        end do
     end if

     ! Read useful planes
     ipart=1
     ipart_grafic=i3_min*plane_size
     do i3=i3_min,i3_max
        read(10)((init_plane(i1,i2),i1=1,g%n1(r%levelmin)),i2=1,g%n2(r%levelmin))
        if(read_pos)then
           read(11)((init_plane_x(i1,i2),i1=1,g%n1(r%levelmin)),i2=1,g%n2(r%levelmin))
        endif
        ! Rescale initial displacement field to code units
        init_plane = init_plane*g%dfact(r%levelmin)*dx/g%dxini(r%levelmin)/g%vfact(r%levelmin)
        if(read_pos)then
           init_plane_x = init_plane_x/g%boxlen_ini
        endif
        do i2=1,g%n2(r%levelmin)
           do i1=1,g%n1(r%levelmin)
              keep_part=(ipart_grafic>=start_ind(g%myid).AND.ipart<=p%npart)
              if(keep_part)then
                 p%vp(ipart,idim)=init_plane(i1,i2)
                 if(.not. read_pos)then
                    dispmax=max(dispmax,abs(init_plane(i1,i2)/dx))
                 else
                    p%xp(ipart,idim)=p%xp(ipart,idim)+init_plane_x(i1,i2)
                    dispmax=max(dispmax,abs(init_plane_x(i1,i2)/dx))
                 endif
                 ipart=ipart+1
              endif
              ipart_grafic=ipart_grafic+1
           end do
        end do
     end do
     ! End loop over planes
     close(10)
     if(read_pos)close(11)
     
  end do
  ! End loop over dimensions

  ! Deallocate temporary array
  deallocate(init_plane)
  if(read_pos)deallocate(init_plane_x)
  
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

end subroutine input_part_grafic
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

