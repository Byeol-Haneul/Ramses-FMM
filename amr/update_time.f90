!################################################################
!################################################################
!################################################################
!################################################################
subroutine m_update_time(pst,ilevel)
  use amr_parameters, only: dp,n_frw
  use ramses_commons, only: pst_t
  use mdl_module
  implicit none
  type(pst_t)::pst
  integer::ilevel

  ! Local variables
  real::ttend
  real,save::ttstart=0.0
  real(dp)::dt,econs,mcons
  integer::i,itest
  integer,dimension(1:4)::input_array,dummy
  
  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)

  ! Local constants
  dt=g%dtnew(ilevel)
  itest=0

  if(ttstart.eq.0.0)call cpu_time(ttstart)

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
        g%epot_tot_int=g%epot_tot_int+0.5D0*(g%epot_tot_old+g%epot_tot)*log(g%aexp/g%aexp_old)
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
           
        !-------------------------------
        ! Output AMR structure to screen
        !-------------------------------
        write(*,*)'Mesh structure'
        do i=r%levelmin,r%nlevelmax
           if(m%noct_tot(i)>0)write(*,999)i,m%noct_tot(i),m%noct_min(i),m%noct_max(i),m%noct_tot(i)/g%ncpu
        end do
999     format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

        !----------------------------------------------
        ! Output mass and energy conservation to screen
        !----------------------------------------------
        if(r%cooling.or.r%pressure_fix)then
           write(*,778)g%nstep_coarse,mcons,econs,g%epot_tot,g%ekin_tot,g%eint_tot
        else
           write(*,777)g%nstep_coarse,mcons,econs,g%epot_tot,g%ekin_tot
        end if
777     format(' Main step=',i6,' mcons=',1pe9.2,' econs=',1pe9.2,' epot=',1pe9.2,' ekin=',1pe9.2)
778     format(' Main step=',i6,' mcons=',1pe9.2,' econs=',1pe9.2,' epot=',1pe9.2,' ekin=',1pe9.2,' eint=',1pe9.2)

        !----------------------------------------------
        ! Output fine step information and used memory
        !----------------------------------------------
        if(r%pic)then
           write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(r%ngridmax)),&
                & real(100.0D0*dble(p%npart_max)/dble(r%npartmax+1))
        else
           write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(r%ngridmax))
        endif
        itest=1
     end if
     g%output_done=.false.

     !---------------
     ! Exit program
     !---------------
     if(g%t>=r%tout(r%noutput).or.g%aexp>=r%aout(r%noutput).or.g%nstep_coarse>=r%nstepmax)then
        write(*,*)'Run completed'
        do i=r%levelmin,r%nlevelmax
           call write_screen(m,i)
        end do
        call cpu_time(ttend)
        write(*,*)'Total elapsed time:',ttend-ttstart
        call mdl_abort(mdl)
     end if

  end if
  g%nstep_coarse_old=g%nstep_coarse

  !----------------------------
  ! Output controls to screen
  !----------------------------
  if(mod(g%nstep,r%ncontrol)==0)then
     if(itest==0)then
        if(r%pic)then
           write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(r%ngridmax)),&
                & real(100.0D0*dble(p%npart_max)/dble(r%npartmax+1))
        else
           write(*,888)g%nstep,g%t,dt,g%aexp,real(100.0D0*dble(m%noct_used_max)/dble(r%ngridmax))
        endif
     end if
  end if
888 format(' Fine step=',i6,' t=',1pe12.5,' dt=',1pe10.3,' a=',1pe10.3,' mem=',0pF4.1,'% ',0pF4.1,'%')
 
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

  ! Broadcast aexp and hexp to all CPUs
  input_array(1:2)=transfer(g%aexp,input_array)
  input_array(3:4)=transfer(g%hexp,input_array)
  call r_broadcast_aexp(pst,input_array,4,dummy,0)

  end associate

end subroutine m_update_time
!##############################################################
!##############################################################
!##############################################################
!##############################################################
recursive subroutine r_broadcast_aexp(pst,input_array,input_size,output_array,output_size)
  use ramses_commons, only: pst_t
  use mdl_parameters
  use mdl_module
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  real(kind=8)::aexp

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_BROADCAST_AEXP,pst%iUpper+1,input_size,output_size,input_array)
     call r_broadcast_aexp(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  else
     pst%s%g%aexp=transfer(input_array(1:2),aexp)
     pst%s%g%hexp=transfer(input_array(3:4),aexp)
  endif

end subroutine r_broadcast_aexp
!##############################################################
!##############################################################
!##############################################################
!##############################################################
recursive subroutine r_clean_stop(pst,input_array,input_size,output_array,output_size)
  use ramses_commons, only: pst_t
  use mdl_parameters
  use mdl_module
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  if(pst%nLower>0)then
     call mdl_send_request(pst%s%mdl,MDL_CLEAN_STOP,pst%iUpper+1,input_size,output_size,input_array)
     call r_clean_stop(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,pst%iUpper+1,output_size)
  endif
  
end subroutine r_clean_stop
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine writemem(usedmem)
  real::usedmem
!  integer::getpagesize
  integer::ipagesize

#ifdef NOSYSTEM
!  call PXFSYSCONF(_SC_PAGESIZE,ipagesize,ierror)
  ipagesize=4096
#else
!  ipagesize = getpagesize()
  ipagesize=4096
#endif
  usedmem=usedmem*ipagesize

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
  character(len=300) :: dir, dir2, file
  integer::ind,j,nmem,read_status
  file='/proc/self/stat'
  open(unit=1,file=file,form='formatted',err=101)
  read(1,'(A300)',IOSTAT=read_status)dir
101  close(1)
  if (read_status < 0)then
     outmem=0.
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
     outmem=nmem
  end if

end subroutine getmem





