!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part_ascii(pst)
  use amr_parameters, only: dp,i8b
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from a Ramses restart file.
  !--------------------------------------------------------------------
  real(dp)::xx1,xx2,xx3,vv1,vv2,vv3,mm1
  integer(i8b)::npart_tot
  character(LEN=80)::filename
  integer::dummy
  integer,allocatable,dimension(:)::input_array

  associate(s=>pst%s)
  
  if(s%r%nrestart>0)return
  if(s%r%verbose)write(*,*)'Entering init_part_ascii'
  
  ! Compute total number of particles in file
  if(TRIM(s%r%initfile(s%r%levelmin)).NE.' ')then
     filename=TRIM(s%r%initfile(s%r%levelmin))//'/ic_part'
     write(*,*)'Opening file '//TRIM(filename)
     open(10,file=filename,form='formatted')
     npart_tot=0
     do
        read(10,*,end=101)xx1,xx2,xx3,vv1,vv2,vv3,mm1
        if(ABS(xx1)<s%r%boxlen/2.0d0.AND.ABS(xx2)<s%r%boxlen/2.0d0.AND.ABS(xx3)<s%r%boxlen/2.0d0)then
           npart_tot=npart_tot+1
        endif
     end do
101  continue
     s%p%npart_tot=npart_tot
     write(*,*)'Found npart_tot=',s%p%npart_tot
     close(10)
  else
     s%p%npart_tot=0
  endif

  ! If no particle found, no need to read
  if(s%p%npart_tot==0)then
     return
  endif

  ! Call recursive slave routine
  allocate(input_array(1:storage_size(npart_tot)/32))
  input_array=transfer(npart_tot,input_array)
  call r_input_part_ascii(pst,input_array,storage_size(npart_tot)/32,dummy,0)
  deallocate(input_array)

  end associate
  
end subroutine m_input_part_ascii
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_input_part_ascii(pst,input_array,input_size,output_array,output_size)
  use mdl_module
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
     call mdl_send_request(pst%s%mdl,MDL_INPUT_PART_ASCII,pst%iUpper+1,input_size,output_size,input_array)
     call r_input_part_ascii(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     npart_tot=transfer(input_array,npart_tot)
     call input_part_ascii(pst%s%r,pst%s%g,pst%s%p,npart_tot)
  endif

end subroutine r_input_part_ascii
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_part_ascii(r,g,p,npart_tot)
  use mdl_module
  use amr_parameters, only: dp,i8b
  use amr_commons, only: run_t,global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  integer(i8b)::npart_tot
  !------------------------------------------------------------
  ! Allocate particle-based arrays.
  ! Read particles positions and velocities from various files
  ! including gadget, ascii or restart files.
  ! grafic initial conditions are performed after the AMR grid 
  ! has been constructed.
  !------------------------------------------------------------
  integer::jpart_loc
  integer::i,ilun,icpu
  integer(i8b)::jpart,indglob
  integer(i8b),dimension(1:g%ncpu+1)::start_ind
  real(dp)::xx1,xx2,xx3,vv1,vv2,vv3,mm1
  character(LEN=80)::filename

  !--------------------------------------
  ! Compute starting index for each cpu
  !--------------------------------------
  p%npart_tot=npart_tot
  do icpu=1,g%ncpu+1
     start_ind(icpu)=1+((icpu-1)*npart_tot)/g%ncpu
  end do
  
  !--------------------------------------
  ! Read ASCII initial conditions file
  !--------------------------------------  
  filename=TRIM(r%initfile(r%levelmin))//'/ic_part'
  open(10+g%myid,file=filename,form='formatted')
  jpart=0
  indglob=0
  jpart_loc=0
  do 
     read(10,*,end=100)xx1,xx2,xx3,vv1,vv2,vv3,mm1
     if(ABS(xx1)<r%boxlen/2.0d0.AND.ABS(xx2)<r%boxlen/2.0d0.AND.ABS(xx3)<r%boxlen/2.0d0)then
        jpart=jpart+1
        indglob=indglob+1
        if(jpart >= start_ind(g%myid) .and. jpart < start_ind(g%myid+1))then
           jpart_loc=jpart_loc+1
           if(jpart_loc>r%npartmax)then
              write(*,*)'Maximum number of particles incorrect'
              write(*,*)'npartmax should be greater than',start_ind(2)
              call mdl_abort(g%mdl)
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
100 continue
  close(10)
     
  p%npart=jpart_loc
  
  ! Put all particles in levelmin
  p%headp=p%npart+1
  p%tailp=p%npart
  p%headp(r%levelmin)=1
  p%tailp(r%levelmin)=p%npart
        
end subroutine input_part_ascii
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

