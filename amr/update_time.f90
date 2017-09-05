!=======================================================================
real(kind=8) function wallclock()
  implicit none
#ifdef WITHOUTMPI
  integer,      save :: tstart
  integer            :: tcur
  integer            :: count_rate
#else
  integer            :: info
  real(kind=8), save :: tstart
  real(kind=8)       :: tcur
#endif
  logical,      save :: first_call=.true.
  real(kind=8), save :: norm, offset=0.
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
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
!=======================================================================
module timer_m
  implicit none
  integer,            parameter         :: mtimer=200                    ! max nr of timers
  real(kind=8),       dimension(mtimer) :: start, time
  integer                               :: ntimer=0, itimer
  character(len=72), dimension(mtimer)  :: labels
contains
!-----------------------------------------------------------------------
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
end subroutine
end module
!=======================================================================
subroutine timer (label, cmd)
  use timer_m
  implicit none
  character(len=*) label, cmd
  real(kind=8) wallclock, current
  integer ierr
!-----------------------------------------------------------------------
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
end subroutine
!=======================================================================
subroutine finalize_timer
  use amr_parameters
  use amr_commons
  use timer_m
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(kind=8) :: total, gtotal, avtime, rmstime
  real(kind=8), dimension(ncpu) :: vtime
  integer,      dimension(ncpu) :: all_ntimer
  logical,      dimension(ncpu) :: gprint_timer
  integer      :: imn, imx, mpi_err, icpu
  logical      :: o, print_timer
!-----------------------------------------------------------------------
  o = myid == 1
  total = 1e-9
  if (o .and. ncpu==1) write (*,'(/a,i7,a)') '     seconds         %    STEP (rank=',myid,')'
  do itimer = 1,ntimer
     total = total + time(itimer)
  end do
  if (ncpu==1) then
     do itimer = 1,ntimer
        if (o .and. time(itimer)/total > 0.001) write (*,'(f12.3,4x,f6.1,4x,a24)') &
          time(itimer), 100.*time(itimer)/total,labels(itimer)
     end do
     if (o) write (*,'(f12.3,4x,f6.1,4x,a)') total, 100., 'TOTAL'
  end if
#ifndef WITHOUTMPI
  if (ncpu > 1) then
     ! Check that timers are consistent across ranks
     call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call MPI_GATHER(ntimer,1,MPI_INTEGER,all_ntimer,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
     if (o) then
        if (maxval(all_ntimer) .ne. minval(all_ntimer)) then
           write (*,*)
           write (*,*) '--------------------------------------------------------------------'
           write (*,*) 'Error: Inconsistent number of timers on each rank. Min, max nr:', minval(all_ntimer), maxval(all_ntimer)
           write (*,*) 'Timing summary below can be misleading'
           write (*,*) 'Labels of timer on rank==1 :'
           write (*,*) '--------------------------------------------------------------------'
           do itimer=1,ntimer
              write(*,'(i3,1x,a)') itimer, labels(itimer)
           enddo
        endif
        ! Find first occurence of a rank with a different number of timers -- if it exists
        gprint_timer=.false.
        do icpu=1,ncpu
           if (all_ntimer(icpu) .ne. ntimer) then
              gprint_timer(icpu) = .true.
              exit
           endif
        enddo
        if (any(gprint_timer)) call sleep(1) ! Make sure that master rank finished, before we print from other rank.
     endif
     call MPI_SCATTER(gprint_timer,1,MPI_LOGICAL,print_timer,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_err)
     if (print_timer) then
        write (*,*)
        write (*,*) 'Labels of timer on rank==',myid
        write (*,*) '--------------------------------------------------------------------'
        do itimer=1,ntimer
           write(*,'(i3,1x,a)') itimer, labels(itimer)
        enddo
        write (*,*)
     endif

     call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call MPI_ALLREDUCE(total,gtotal,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
     gtotal = gtotal / ncpu

     if (o) write (*,*) '--------------------------------------------------------------------'
     if (o) write (*,'(/a)') '     minimum       average       maximum' // &
                  '  standard dev        std/av       %   rmn   rmx  TIMER'
     do itimer = 1,ntimer
        call MPI_GATHER(real(time(itimer),kind=8),1,MPI_REAL8,vtime,1,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
        if (o) then
           if (maxval(vtime)/gtotal > 0.001) then
              avtime  = sum(vtime) / ncpu ! average time used
              imn     = minloc(vtime,1)
              imx     = maxloc(vtime,1)
              rmstime = sqrt(sum((vtime - avtime)**2)/ncpu)
              write (*,'(5(f12.3,2x),f6.1,2x,2i4,4x,a24)') &
                 vtime(imn), avtime, vtime(imx), rmstime, rmstime/avtime, 100.*avtime/gtotal, imn, imx, labels(itimer)
           endif
        endif
     end do
     if (o) write (*,'(f12.3,4x,f6.1,4x,a)') total, 100., 'TOTAL'
  endif
#endif
end subroutine
!=======================================================================
subroutine reset_timer
   use timer_m
   implicit none
#ifndef WITHOUTMPI
   include 'mpif.h'
#endif
!-----------------------------------------------------------------------
   do itimer = 1,ntimer
      time(itimer)=0.0
   end do
end subroutine
!=======================================================================
subroutine update_time(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif  
  integer::ilevel

  real(dp)::dt,econs,mcons
  real(kind=8)::ttend
  real(kind=8),save::ttstart=0
  integer::i,itest,info

  ! Local constants
  dt=dtnew(ilevel)
  itest=0

#ifndef WITHOUTMPI
  if(myid==1)then
     if(ttstart.eq.0.0)ttstart=MPI_WTIME(info)
  endif
#endif

  !-------------------------------------------------------------
  ! At this point, IF nstep_coarse has JUST changed, all levels
  ! are synchronised, and all new refinements have been done.
  !-------------------------------------------------------------
  if(nstep_coarse .ne. nstep_coarse_old)then

     !--------------------------
     ! Check mass conservation
     !--------------------------
     if(mass_tot_0==0.0D0)then
        mass_tot_0=mass_tot
        mcons=0.0D0
     else
        mcons=(mass_tot-mass_tot_0)/mass_tot_0
     end if

     !----------------------------
     ! Check energy conservation
     !----------------------------
     if(epot_tot_old.ne.0)then
        epot_tot_int=epot_tot_int + &
             & 0.5D0*(epot_tot_old+epot_tot)*log(aexp/aexp_old)
     end if
     epot_tot_old=epot_tot
     aexp_old=aexp
     if(const==0.0D0)then
        const=epot_tot+ekin_tot  ! initial total energy
        econs=0.0D0
     else
        econs=(ekin_tot+epot_tot-epot_tot_int-const) / &
             &(-(epot_tot-epot_tot_int-const)+ekin_tot)
     end if

     if(mod(nstep_coarse,ncontrol)==0.or.output_done)then
        if(myid==1)then
           
           !-------------------------------
           ! Output AMR structure to screen
           !-------------------------------
           write(*,*)'Mesh structure'
           do i=levelmin,nlevelmax
              if(noct_tot(i)>0)write(*,999)i,noct_tot(i),noct_min(i),noct_max(i),noct_tot(i)/ncpu
           end do
           !----------------------------------------------
           ! Output mass and energy conservation to screen
           !----------------------------------------------
           if(cooling.or.pressure_fix)then
              write(*,778)nstep_coarse,mcons,econs,epot_tot,ekin_tot,eint_tot
           else
              write(*,777)nstep_coarse,mcons,econs,epot_tot,ekin_tot
           end if
           if(pic)then
              write(*,888)nstep,t,dt,aexp,real(100.0D0*dble(noct_used_max)/dble(ngridmax)),&
                   & real(100.0D0*dble(npart_max)/dble(npartmax+1))
           else
              write(*,888)nstep,t,dt,aexp,real(100.0D0*dble(noct_used_max)/dble(ngridmax))
           endif
           itest=1
        end if
        output_done=.false.
     end if

     !---------------
     ! Exit program
     !---------------
     if(t>=tout(noutput).or.aexp>=aout(noutput).or. &
          & nstep_coarse>=nstepmax)then
        if(myid==1)then
           write(*,*)'Run completed'
           do i=levelmin,nlevelmax
              call write_screen(i)
           end do           
#ifndef WITHOUTMPI
           ttend=MPI_WTIME(info)
           write(*,*)'Total elapsed time:',ttend-ttstart
#endif
        endif
        call clean_stop
     end if

  end if
  nstep_coarse_old=nstep_coarse

  !----------------------------
  ! Output controls to screen
  !----------------------------
  if(mod(nstep,ncontrol)==0)then
     if(myid==1.and.itest==0)then
        if(pic)then
           write(*,888)nstep,t,dt,aexp,real(100.0D0*dble(noct_used_max)/dble(ngridmax)),&
                & real(100.0D0*dble(npart_max)/dble(npartmax+1))
        else
           write(*,888)nstep,t,dt,aexp,real(100.0D0*dble(noct_used_max)/dble(ngridmax))
        endif
     end if
  end if

  !------------------------
  ! Update time variables
  !------------------------
  t=t+dt
  nstep=nstep+1
  if(cosmo)then
     ! Find neighboring times
     i=1
     do while(tau_frw(i)>t.and.i<n_frw)
        i=i+1
     end do
     ! Interpolate expansion factor
     aexp = aexp_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          & aexp_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
     hexp = hexp_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          & hexp_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
     texp =    t_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          &    t_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
  else
     aexp = 1.0
     hexp = 0.0
     texp = t
  end if

777 format(' Main step=',i6,' mcons=',1pe9.2,' econs=',1pe9.2, &
         & ' epot=',1pe9.2,' ekin=',1pe9.2)
778 format(' Main step=',i6,' mcons=',1pe9.2,' econs=',1pe9.2, &
         & ' epot=',1pe9.2,' ekin=',1pe9.2,' eint=',1pe9.2)
888 format(' Fine step=',i6,' t=',1pe12.5,' dt=',1pe10.3, &
         & ' a=',1pe10.3,' mem=',0pF4.1,'% ',0pF4.1,'%')
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')
 
end subroutine update_time
!=======================================================================
subroutine update_time_2(r,g,m,p,ilevel)
  use amr_parameters, only: dp,n_frw
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif  
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  integer::ilevel

  ! Local variables
  real(dp)::dt,econs,mcons
  real(kind=8)::ttend
  real(kind=8),save::ttstart=0
  integer::i,itest,info

  ! Local constants
  dt=g%dtnew(ilevel)
  itest=0

#ifndef WITHOUTMPI
  if(g%myid==1)then
     if(ttstart.eq.0.0)ttstart=MPI_WTIME(info)
  endif
#endif

  !-------------------------------------------------------------
  ! At this point, IF nstep_coarse has JUST changed, all levels
  ! are synchronised, and all new refinements have been done.
  !-------------------------------------------------------------
  if(g%nstep_coarse .ne. g%nstep_coarse_old)then

     !--------------------------
     ! Check mass conservation
     !--------------------------
     if(g%mass_tot_0==0.0D0)then
        g%mass_tot_0=g%mass_tot
        mcons=0.0D0
     else
        mcons=(g%mass_tot-g%mass_tot_0)/g%mass_tot_0
     end if

     !----------------------------
     ! Check energy conservation
     !----------------------------
     if(g%epot_tot_old.ne.0)then
        g%epot_tot_int=g%epot_tot_int + &
             & 0.5D0*(g%epot_tot_old+g%epot_tot)*log(g%aexp/g%aexp_old)
     end if
     g%epot_tot_old=g%epot_tot
     g%aexp_old=g%aexp
     if(g%const==0.0D0)then
        g%const=g%epot_tot+g%ekin_tot  ! initial total energy
        econs=0.0D0
     else
        econs=(g%ekin_tot+g%epot_tot-g%epot_tot_int-g%const) / &
             &(-(g%epot_tot-g%epot_tot_int-g%const)+g%ekin_tot)
     end if

     if(mod(g%nstep_coarse,r%ncontrol)==0.or.g%output_done)then
        if(g%myid==1)then
           
           !-------------------------------
           ! Output AMR structure to screen
           !-------------------------------
           write(*,*)'Mesh structure'
           do i=r%levelmin,r%nlevelmax
              if(m%noct_tot(i)>0)write(*,999)i,m%noct_tot(i),m%noct_min(i),m%noct_max(i),m%noct_tot(i)/g%ncpu
           end do
           !----------------------------------------------
           ! Output mass and energy conservation to screen
           !----------------------------------------------
           if(r%cooling.or.r%pressure_fix)then
              write(*,778)g%nstep_coarse,mcons,econs,g%epot_tot,g%ekin_tot,g%eint_tot
           else
              write(*,777)g%nstep_coarse,mcons,econs,g%epot_tot,g%ekin_tot
           end if
           if(r%pic)then
              write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(r%ngridmax)),&
                   & real(100.0D0*dble(p%npart_max)/dble(r%npartmax+1))
           else
              write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(r%ngridmax))
           endif
           itest=1
        end if
        g%output_done=.false.
     end if

     !---------------
     ! Exit program
     !---------------
     if(g%t>=r%tout(r%noutput).or.g%aexp>=r%aout(r%noutput).or. &
          & g%nstep_coarse>=r%nstepmax)then
        if(g%myid==1)then
           write(*,*)'Run completed'
           do i=r%levelmin,r%nlevelmax
              call write_screen_2(r,g,m,i)
           end do           
#ifndef WITHOUTMPI
           ttend=MPI_WTIME(info)
           write(*,*)'Total elapsed time:',ttend-ttstart
#endif
        endif
        call clean_stop
     end if

  end if
  g%nstep_coarse_old=g%nstep_coarse

  !----------------------------
  ! Output controls to screen
  !----------------------------
  if(mod(g%nstep,r%ncontrol)==0)then
     if(g%myid==1.and.itest==0)then
        if(r%pic)then
           write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(r%ngridmax)),&
                & real(100.0D0*dble(p%npart_max)/dble(r%npartmax+1))
        else
           write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(r%ngridmax))
        endif
     end if
  end if

  !------------------------
  ! Update time variables
  !------------------------
  g%t=g%t+dt
  g%nstep=g%nstep+1
  if(r%cosmo)then
     ! Find neighboring times
     i=1
     do while(g%tau_frw(i)>g%t.and.i<n_frw)
        i=i+1
     end do
     ! Interpolate expansion factor
     g%aexp = g%aexp_frw(i  )*(g%t-g%tau_frw(i-1))/(g%tau_frw(i  )-g%tau_frw(i-1))+ &
            & g%aexp_frw(i-1)*(g%t-g%tau_frw(i  ))/(g%tau_frw(i-1)-g%tau_frw(i  ))
     g%hexp = g%hexp_frw(i  )*(g%t-g%tau_frw(i-1))/(g%tau_frw(i  )-g%tau_frw(i-1))+ &
            & g%hexp_frw(i-1)*(g%t-g%tau_frw(i  ))/(g%tau_frw(i-1)-g%tau_frw(i  ))
     g%texp =    g%t_frw(i  )*(g%t-g%tau_frw(i-1))/(g%tau_frw(i  )-g%tau_frw(i-1))+ &
            &    g%t_frw(i-1)*(g%t-g%tau_frw(i  ))/(g%tau_frw(i-1)-g%tau_frw(i  ))
  else
     g%aexp = 1.0
     g%hexp = 0.0
     g%texp = g%t
  end if

777 format(' Main step=',i6,' mcons=',1pe9.2,' econs=',1pe9.2, &
         & ' epot=',1pe9.2,' ekin=',1pe9.2)
778 format(' Main step=',i6,' mcons=',1pe9.2,' econs=',1pe9.2, &
         & ' epot=',1pe9.2,' ekin=',1pe9.2,' eint=',1pe9.2)
888 format(' Fine step=',i6,' t=',1pe12.5,' dt=',1pe10.3, &
         & ' a=',1pe10.3,' mem=',0pF4.1,'% ',0pF4.1,'%')
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')
 
end subroutine update_time_2
  
subroutine clean_stop
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::info

  call finalize_timer

#ifndef WITHOUTMPI
  call MPI_FINALIZE(info)
#endif
  stop
end subroutine clean_stop

subroutine clean_abort
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::info
#ifndef WITHOUTMPI
     call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
     stop
#endif
end subroutine clean_abort

subroutine writemem(usedmem)
  real::usedmem
  integer::getpagesize

#ifdef NOSYSTEM
!  call PXFSYSCONF(_SC_PAGESIZE,ipagesize,ierror)
  ipagesize=4096
#else
!  ipagesize = getpagesize()
  ipagesize=4096
#endif
  usedmem=dble(usedmem)*dble(ipagesize)

  if(usedmem>1024.**3.)then
     write(*,999)usedmem/1024.**3.
  else if (usedmem>1024.**2.) then
     write(*,998)usedmem/1024.**2
  else if (usedmem>1024.) then
     write(*,997)usedmem/1024.
  endif

997 format(' Used memory:',F6.1,' kb')
998 format(' Used memory:',F6.1,' Mb')
999 format(' Used memory:',F6.1,' Gb')

end subroutine writemem

subroutine getmem(outmem)
  real::outmem
  character(len=300) :: dir, dir2,  cmd, file
  integer::read_status
  file='/proc/self/stat'
  open(unit=1,file=file,form='formatted',err=101)
  read(1,'(A300)',IOSTAT=read_status)dir
101  close(1)
  if (read_status < 0)then
     outmem=dble(0.)
     write(*,*)'Problem in checking free memory'
  else
     ind=300
     j=0
     do while (j<23)
        ind=index(dir,' ')
        dir2=dir(ind+1:300)
        j=j+1
        dir=dir2
     end do
     ind=index(dir,' ')
     dir2=dir(1:ind)
     read(dir2,'(I12)')nmem
     outmem=dble(nmem)
  end if

end subroutine getmem

subroutine cmpmem(outmem)
  use amr_commons
  use hydro_commons
  implicit none

  real::outmem,outmem_int,outmem_dp,outmem_qdp
  outmem_int=0.0
  outmem_dp=0.0
  outmem_qdp=0.0

  outmem_dp =outmem_dp +ngridmax*ndim      ! xg
  outmem_int=outmem_int+ngridmax*twondim   ! nbor
  outmem_int=outmem_int+ngridmax           ! father
  outmem_int=outmem_int+ngridmax           ! next
  outmem_int=outmem_int+ngridmax           ! prev
  outmem_int=outmem_int+ngridmax*twotondim ! son 
  outmem_int=outmem_int+ngridmax*twotondim ! flag1
  outmem_int=outmem_int+ngridmax*twotondim ! flag2
  outmem_int=outmem_int+ngridmax*twotondim ! cpu_map1
  outmem_int=outmem_int+ngridmax*twotondim ! cpu_map2
  outmem_qdp=outmem_qdp+ngridmax*twotondim ! hilbert_key

  ! Add communicator variable here

  if(hydro)then
     
  outmem_dp =outmem_dp +ngridmax*twotondim*nvar ! uold
  outmem_dp =outmem_dp +ngridmax*twotondim*nvar ! unew

  if(pressure_fix)then

  outmem_dp =outmem_dp +ngridmax*twotondim ! uold
  outmem_dp =outmem_dp +ngridmax*twotondim ! uold

  endif

  endif

  write(*,*)'Estimated memory=',(outmem_dp*8.+outmem_int*4.+outmem_qdp*8.)/1024./1024.

end subroutine cmpmem






