!################################################################
!################################################################
!################################################################
!################################################################
!=======================================================================
real(kind=8) function wallclock()
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
#ifdef WITHOUTMPI
  integer,      save :: tstart
  integer            :: tcur
  integer            :: count_rate
#else
  real(kind=8), save :: tstart
  real(kind=8)       :: tcur
#endif
  logical,      save :: first_call=.true.
  real(kind=8), save :: norm, offset=0.
  !---------------------------------------------------------------------
  if (first_call) then
#ifdef WITHOUTMPI
     call system_clock(count=tstart, count_rate=count_rate)
     norm=1d0/count_rate
#else
     norm = 1d0
     tstart = MPI_Wtime()
#endif
     first_call=.false.
  end if
#ifdef WITHOUTMPI
  call system_clock(count=tcur)
#else
  tcur = MPI_Wtime()
#endif
  wallclock = (tcur-tstart)*norm + offset
  if (wallclock < 0.) then
     offset = offset + 24d0*3600d0
     wallclock = wallclock + 24d0*3600d0
  end if
end function wallclock
!################################################################
!################################################################
!################################################################
!################################################################
module timer_module
  implicit none
  integer,            parameter         :: mtimer=200       ! max nr of timers
  real(kind=8),       dimension(mtimer) :: start, time
  integer                               :: ntimer=0, itimer
  character(len=72), dimension(mtimer)  :: labels
contains
  subroutine findit (label)
    implicit none
    character(len=*) label
    do itimer=1,ntimer
       if (trim(label) == trim(labels(itimer))) return
    end do
    ntimer = ntimer+1
    itimer = ntimer
    labels(itimer) = label
    time(itimer) = 0.
  end subroutine findit
end module timer_module
!################################################################
!################################################################
!################################################################
!################################################################
subroutine m_timer(pst,label,cmd)
  use mdl_module
  use ramses_commons, only: pst_t
  use timer_module
  implicit none
  type(pst_t)::pst
  character(len=*) label, cmd
  real(kind=8) wallclock, current
  current = wallclock()                                                 ! current time
  if (itimer > 0) then                                                  ! if timer is active ..
     time(itimer) = time(itimer) + current - start(itimer)              ! add to it
  end if
  call findit (label)                                                   ! locate timer slot
  if (cmd == 'start') then                                              ! start command
     start(itimer) = current                                            ! register start time
  else if (cmd == 'stop') then                                          ! stop command
     itimer = 0                                                         ! turn off timer
  end if
end subroutine m_timer
!################################################################
!################################################################
!################################################################
!################################################################
subroutine m_output_timer(pst,write_file,filename)
  use amr_parameters, only: flen
  use mdl_module
  use ramses_commons, only: pst_t
  use timer_module
  implicit none
  type(pst_t)::pst
  real(kind=8)::total
  character(len=*)::filename
  character(len=flen)::fileloc
  logical::write_file
  integer::ilun
  if(write_file)then
     ilun=11
     fileloc=TRIM(filename) ! Open file for timing info
     open(unit=ilun,file=fileloc,form='formatted')
  else
     ilun=6
  endif
  write (ilun,'(/a,i7,a)') '     seconds         %    STEP'
  total = 1e-9
  do itimer = 1,ntimer
     total = total + time(itimer)
  end do
  do itimer = 1,ntimer
     if (time(itimer)/total > 0.0001) write (ilun,'(f12.3,4x,f6.1,4x,a24)') &
          time(itimer), 100.*time(itimer)/total,labels(itimer)
  end do
  write (ilun,'(f12.3,4x,f6.1,4x,a)') total, 100., 'TOTAL'
  write(ilun,*)
  if(write_file)close(ilun)
end subroutine m_output_timer
!################################################################
!################################################################
!################################################################
!################################################################


  









