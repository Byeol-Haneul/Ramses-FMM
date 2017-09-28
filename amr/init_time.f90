!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_time(r,g)
  use amr_parameters, only: n_frw
  use amr_commons, only: run_t,global_t
  implicit none
  type(run_t)::r
  type(global_t)::g

  ! Local variables
  integer::i

  if(r%verbose)write(*,*)'Entering init_time'

  if(r%nrestart==0)then
     if(r%cosmo)then
        ! Get cosmological parameters from input files
        call init_cosmo(r,g)
     else
        ! Get parameters from input files
        if(r%initfile(r%levelmin).ne.' '.and.r%filetype.eq.'grafic')then
           call init_file(r,g)
        endif
        g%t=0.0
        g%aexp=1.0
     end if
  end if

  if(r%cosmo)then

     ! Compute Friedman model look up table
     if(g%myid==1)write(*,*)'Computing Friedman model'
     call friedman(dble(g%omega_m),dble(g%omega_l),dble(g%omega_k), &
          & 1.d-6,dble(g%aexp_ini), &
          & g%aexp_frw,g%hexp_frw,g%tau_frw,g%t_frw,n_frw)

     ! Compute initial conformal time                                    
     ! Find neighboring expansion factors                                  
     i=1
     do while(g%aexp_frw(i)>g%aexp.and.i<n_frw)
        i=i+1
     end do                                                                
     ! Interploate time                                                    
     if(r%nrestart==0)then                                                   
        g%t=g%tau_frw(i)*(g%aexp-g%aexp_frw(i-1))/(g%aexp_frw(i)-g%aexp_frw(i-1))+ &   
             & g%tau_frw(i-1)*(g%aexp-g%aexp_frw(i))/(g%aexp_frw(i-1)-g%aexp_frw(i)) 
        g%aexp=g%aexp_frw(i)*(g%t-g%tau_frw(i-1))/(g%tau_frw(i)-g%tau_frw(i-1))+ &     
             & g%aexp_frw(i-1)*(g%t-g%tau_frw(i))/(g%tau_frw(i-1)-g%tau_frw(i))      
        g%hexp=g%hexp_frw(i)*(g%t-g%tau_frw(i-1))/(g%tau_frw(i)-g%tau_frw(i-1))+ &     
             & g%hexp_frw(i-1)*(g%t-g%tau_frw(i))/(g%tau_frw(i-1)-g%tau_frw(i))      
     end if
     g%texp=g%t_frw(i)*(g%t-g%tau_frw(i-1))/(g%tau_frw(i)-g%tau_frw(i-1))+ &           
          & g%t_frw(i-1)*(g%t-g%tau_frw(i))/(g%tau_frw(i-1)-g%tau_frw(i))            
  else
     g%texp=g%t
  end if                                                                   

end subroutine init_time
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_file(r,g)
  use amr_parameters, only: sp
  use amr_commons, only: run_t,global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  !------------------------------------------------------
  ! Read geometrical parameters in the initial condition files.
  ! Initial conditions are supposed to be made by 
  ! Bertschinger's grafic version 2.0 code.
  !------------------------------------------------------
  integer:: ilevel
  real(sp)::dxini0,xoff10,xoff20,xoff30,astart0,omega_m0,omega_l0,h00
  character(LEN=80)::filename
  logical::ok

  if(r%verbose)write(*,*)'Entering init_file'

  ! Reading initial conditions parameters only
  g%nlevelmax_part=r%levelmin-1
  do ilevel=r%levelmin,r%nlevelmax
     if(r%initfile(ilevel).ne.' ')then
        filename=TRIM(r%initfile(ilevel))//'/ic_d'
        INQUIRE(file=filename,exist=ok)
        if(.not.ok)then
           if(g%myid==1)then
              write(*,*)TRIM(r%initfile(ilevel))
              write(*,*)'File '//TRIM(filename)//' does not exist'
           end if
           call clean_stop(g)
        end if
        open(10,file=filename,form='unformatted')
        if(g%myid==1)write(*,*)'Reading file '//TRIM(filename)
        rewind 10
        read(10)g%n1(ilevel),g%n2(ilevel),g%n3(ilevel),dxini0 &
             & ,xoff10,xoff20,xoff30 &
             & ,astart0,omega_m0,omega_l0,h00
        close(10)
        g%dxini(ilevel)=dxini0
        g%xoff1(ilevel)=xoff10
        g%xoff2(ilevel)=xoff20
        g%xoff3(ilevel)=xoff30
        g%nlevelmax_part=g%nlevelmax_part+1
     endif
  end do

  ! Check compatibility with run parameters
  if(         g%n1(r%levelmin).NE.2**r%levelmin &
       & .or. g%n2(r%levelmin).NE.2**r%levelmin &
       & .or. g%n3(r%levelmin).NE.2**r%levelmin) then 
     write(*,*)'coarser grid is not compatible with initial conditions file'
     write(*,*)'Found    n1=',g%n1(r%levelmin),&
          &            ' n2=',g%n2(r%levelmin),&
          &            ' n3=',g%n3(r%levelmin)
     write(*,*)'Expected n1=',2**r%levelmin &
          &           ,' n2=',2**r%levelmin &
          &           ,' n3=',2**r%levelmin
     call clean_stop(g)
  end if

  ! Write initial conditions parameters
  if(g%myid==1)then
     do ilevel=r%levelmin,g%nlevelmax_part
        write(*,'(" Initial conditions for level =",I4)')ilevel
        write(*,'(" n1=",I4," n2=",I4," n3=",I4)') &
             & g%n1(ilevel),&
             & g%n2(ilevel),&
             & g%n3(ilevel)
        write(*,'(" dx=",1pe10.3)')g%dxini(ilevel)
        write(*,'(" xoff=",1pe10.3," yoff=",1pe10.3," zoff=",1pe10.3)') &
             & g%xoff1(ilevel),&
             & g%xoff2(ilevel),&
             & g%xoff3(ilevel)
     end do
  end if

end subroutine init_file
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_cosmo(r,g)
  use amr_parameters, only: sp,dp,ndim
  use amr_commons, only: run_t,global_t
  use gadgetreadfilemod
  implicit none
  type(run_t)::r
  type(global_t)::g
  !------------------------------------------------------
  ! Read cosmological and geometrical parameters
  ! in the initial condition files.
  ! Initial conditions are supposed to be made by 
  ! Bertschinger's grafic version 2.0 code.
  !------------------------------------------------------
  integer:: ilevel,i
  real(sp)::dxini0,xoff10,xoff20,xoff30,astart0,omega_m0,omega_l0,h00
  real(dp)::d1a,fpeebl
  character(LEN=80)::filename
  character(LEN=5)::nchar
  logical::ok
  TYPE(gadgetheadertype) :: gadgetheader 

  if(r%verbose)write(*,*)'Entering init_cosmo'

  if(r%initfile(r%levelmin)==' ')then
     write(*,*)'You need to specifiy at least one level of initial condition'
     call clean_stop(g)
  end if

  SELECT CASE (r%filetype)
  case ('grafic', 'ascii')

     ! Reading initial conditions parameters only
     g%aexp=2.0
     g%nlevelmax_part=r%levelmin-1
     do ilevel=r%levelmin,r%nlevelmax
        if(r%initfile(ilevel).ne.' ')then
           if(r%multiple)then
              call title(g%myid,nchar)
              filename=TRIM(r%initfile(ilevel))//'/dir_deltab/ic_deltab.'//TRIM(nchar)
           else
              filename=TRIM(r%initfile(ilevel))//'/ic_deltab'
           endif
           INQUIRE(file=filename,exist=ok)
           if(.not.ok)then
              if(g%myid==1)then
                 write(*,*)'File '//TRIM(filename)//' does not exist'
              end if
              call clean_stop(g)
           end if
           open(10,file=filename,form='unformatted')
           if(g%myid==1)write(*,*)'Reading file '//TRIM(filename)
           rewind 10
           read(10)g%n1(ilevel),g%n2(ilevel),g%n3(ilevel),dxini0 &
                & ,xoff10,xoff20,xoff30 &
                & ,astart0,omega_m0,omega_l0,h00
           close(10)
           g%dxini(ilevel)=dxini0
           g%xoff1(ilevel)=xoff10
           g%xoff2(ilevel)=xoff20
           g%xoff3(ilevel)=xoff30
           g%astart(ilevel)=astart0
           g%omega_m=omega_m0
           g%omega_l=omega_l0
           if(r%hydro)g%omega_b=0.045 ! Be careful hard-coded !
           g%h0=h00
           g%aexp=MIN(g%aexp,g%astart(ilevel))
           g%nlevelmax_part=g%nlevelmax_part+1
           ! Compute SPH equivalent mass (initial gas mass resolution)
           r%mass_sph=g%omega_b/g%omega_m*0.5d0**(ndim*ilevel)           
        endif
     end do

     ! Compute initial expansion factor
     if(g%aexp_ini.lt.1.0)then
        g%aexp=g%aexp_ini
     else
        g%aexp_ini=g%aexp
     endif
     
     ! Check compatibility with run parameters
     if(.not. r%multiple) then
        if(         g%n1(r%levelmin).ne.2**r%levelmin &
             & .or. g%n2(r%levelmin).ne.2**r%levelmin &
             & .or. g%n3(r%levelmin).ne.2**r%levelmin) then 
           write(*,*)'coarser grid is not compatible with initial conditions file'
           write(*,*)'Found    n1=',g%n1(r%levelmin),&
                &            ' n2=',g%n2(r%levelmin),&
                &            ' n3=',g%n3(r%levelmin)
           write(*,*)'Expected n1=',2**r%levelmin &
                &           ,' n2=',2**r%levelmin &
                &           ,' n3=',2**r%levelmin
           call clean_stop(g)
        endif
     end if
     
     ! Compute box length in the initial conditions in units of h-1 Mpc
     g%boxlen_ini=2**r%levelmin*g%dxini(r%levelmin)*(g%h0/100.)
     
  CASE ('gadget')

     ! Reading gadget file header only
     if (r%verbose) write(*,*)'Reading in gadget format from '//TRIM(r%initfile(r%levelmin))
     call gadgetreadheader(TRIM(r%initfile(r%levelmin)), 0, gadgetheader, ok)
     if(.not.ok) call clean_stop(g)
     do i=1,6
        if (i .ne. 2) then
           if (gadgetheader%nparttotal(i) .ne. 0) then
              write(*,*) 'Non DM particles present in bin ', i
              call clean_stop(g)
           endif
        endif
     enddo
     if (gadgetheader%mass(2) == 0) then
        write(*,*) 'Particles have different masses, not supported'
        call clean_stop(g)
     endif
     g%omega_m = gadgetheader%omega0
     g%omega_l = gadgetheader%omegalambda
     if(r%hydro)g%omega_b=0.045 ! Be careful hard-coded !
     g%h0 = gadgetheader%hubbleparam * 100.d0
     g%boxlen_ini = gadgetheader%boxsize
     g%aexp = gadgetheader%time
     g%aexp_ini = g%aexp
     ! Compute SPH equivalent mass (initial gas mass resolution)
     r%mass_sph=g%omega_b/g%omega_m*0.5d0**(ndim*r%levelmin)
     g%nlevelmax_part = r%levelmin
     g%astart(r%levelmin) = g%aexp
     g%xoff1(r%levelmin)=0
     g%xoff2(r%levelmin)=0
     g%xoff3(r%levelmin)=0
     g%dxini(r%levelmin) = g%boxlen_ini/(2**r%levelmin*(g%h0/100.0))

  CASE DEFAULT
     write(*,*) 'Unsupported input format '//r%filetype
     call clean_stop(g)
  END SELECT

  ! Write cosmological parameters
  if(g%myid==1)then
     write(*,'(" Cosmological parameters:")')
     write(*,'(" aexp=",1pe10.3," H0=",1pe10.3," km s-1 Mpc-1")')g%aexp,g%h0
     write(*,'(" omega_m=",F7.3," omega_l=",F7.3)')g%omega_m,g%omega_l
     write(*,'(" box size=",1pe10.3," h-1 Mpc")')g%boxlen_ini
  end if
  g%omega_k=1.d0-g%omega_l-g%omega_m

  ! Save cosmological parameters into run parameters
  r%omega_b=g%omega_b
  r%omega_m=g%omega_m
  r%omega_l=g%omega_l
  r%omega_k=g%omega_k
  
  ! Compute linear scaling factor between aexp and astart(ilevel)
  do ilevel=r%levelmin,g%nlevelmax_part
     g%dfact(ilevel)=d1a(g%aexp,g%omega_m,g%omega_l)/d1a(g%astart(ilevel),g%omega_m,g%omega_l)
     g%vfact(ilevel)=g%astart(ilevel)*fpeebl(g%astart(ilevel),g%omega_m,g%omega_l) & ! Same scale factor as in grafic1
          & *sqrt(g%omega_m/g%astart(ilevel)+g%omega_l*g%astart(ilevel)*g%astart(ilevel)+g%omega_k) &
          & /g%astart(ilevel)*g%h0
  end do

  ! Write initial conditions parameters
  do ilevel=r%levelmin,g%nlevelmax_part
     if(g%myid==1)then
        write(*,'(" Initial conditions for level =",I4)')ilevel
        write(*,'(" dx=",1pe10.3," h-1 Mpc")')g%dxini(ilevel)*g%h0/100.
     endif
     if(.not.r%multiple)then
        if(g%myid==1)then
           write(*,'(" n1=",I4," n2=",I4," n3=",I4)') &
                & g%n1(ilevel),&
                & g%n2(ilevel),&
                & g%n3(ilevel)
           write(*,'(" xoff=",1pe10.3," yoff=",1pe10.3," zoff=",1pe10.3," h-1 Mpc")') &
                & g%xoff1(ilevel)*g%h0/100.,&
                & g%xoff2(ilevel)*g%h0/100.,&
                & g%xoff3(ilevel)*g%h0/100.
        endif
     else
        write(*,'(" myid=",I4," n1=",I4," n2=",I4," n3=",I4)') &
             & g%myid,g%n1(ilevel),g%n2(ilevel),g%n3(ilevel)
        write(*,'(" myid=",I4," xoff=",1pe10.3," yoff=",1pe10.3," zoff=",1pe10.3," h-1 Mpc")') &
             & g%myid,g%xoff1(ilevel)*g%h0/100.,&
             & g%xoff2(ilevel)*g%h0/100.,&
             & g%xoff3(ilevel)*g%h0/100.
     endif
  end do

  ! Scale displacement in Mpc to code velocity (v=dx/dtau)
  ! in coarse cell units per conformal time
  g%vfact(1)=g%aexp*fpeebl(g%aexp,g%omega_m,g%omega_l)*sqrt(g%omega_m/g%aexp+g%omega_l*g%aexp*g%aexp+g%omega_k)
  ! This scale factor is different from vfact in grafic by h0/aexp

end subroutine init_cosmo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
     & axp_out,hexp_out,tau_out,t_out,ntable)
  use amr_parameters
  implicit none
  integer::ntable
  real(kind=8)::O_mat_0, O_vac_0, O_k_0
  real(kind=8)::alpha,axp_min
  real(dp),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
  ! ######################################################!
  ! This subroutine assumes that axp = 1 at z = 0 (today) !
  ! and that t and tau = 0 at z = 0 (today).              !
  ! axp is the expansion factor, hexp the Hubble constant !
  ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
  ! time, and t the look-back time, both in unit of 1/H0. !
  ! alpha is the required accuracy and axp_min is the     !
  ! starting expansion factor of the look-up table.       !
  ! ntable is the required size of the look-up table.     !
  ! ######################################################!
  real(kind=8)::axp_tau, axp_t
  real(kind=8)::axp_tau_pre, axp_t_pre
  real(kind=8)::dadtau, dadt
  real(kind=8)::dtau,dt
  real(kind=8)::tau,t
  integer::nstep,nout,nskip

  if( (O_mat_0+O_vac_0+O_k_0) .ne. 1.0D0 )then
     write(*,*)'Error: non-physical cosmological constants'
     write(*,*)'O_mat_0,O_vac_0,O_k_0=',O_mat_0,O_vac_0,O_k_0
     write(*,*)'The sum must be equal to 1.0, but '
     write(*,*)'O_mat_0+O_vac_0+O_k_0=',O_mat_0+O_vac_0+O_k_0
     stop
  end if

  axp_tau = 1.0D0
  axp_t = 1.0D0
  tau = 0.0D0
  t = 0.0D0
  nstep = 0
  
  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau
     
     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
  end do

!  write(*,666)-t
!  666 format(' Age of the Universe (in unit of 1/H0)=',1pe10.3)

  nskip=nstep/ntable
  
  axp_t = 1.d0
  t = 0.d0
  axp_tau = 1.d0
  tau = 0.d0
  nstep = 0
  nout=0
  t_out(nout)=t
  tau_out(nout)=tau
  axp_out(nout)=axp_tau
  hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau

     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
     if(mod(nstep,nskip)==0)then
        nout=nout+1
        t_out(nout)=t
        tau_out(nout)=tau
        axp_out(nout)=axp_tau
        hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
     end if

  end do
  t_out(ntable)=t
  tau_out(ntable)=tau
  axp_out(ntable)=axp_tau
  hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

end subroutine friedman
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0) 
  use amr_parameters
  real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
  dadtau = axp_tau*axp_tau*axp_tau *  &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
       &     O_k_0   * axp_tau )
  dadtau = sqrt(dadtau)
  return
end function dadtau
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
  use amr_parameters
  real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
  dadt   = (1.0D0/axp_t)* &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_t*axp_t*axp_t + &
       &     O_k_0   * axp_t )
  dadt = sqrt(dadt)
  return
end function dadt
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function fy(a,omega_m,omega_l)
  use amr_parameters, only: dp
  implicit none
  !      Computes the integrand
  real(dp)::fy
  real(dp)::y,a,omega_m,omega_l
  
  y=omega_m*(1.d0/a-1.d0) + omega_l*(a*a-1.d0) + 1.d0
  fy=1.d0/y**1.5d0
  
  return
end function fy
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function d1a(a,omega_m,omega_l)
  use amr_parameters, only: dp
  implicit none
  real(dp)::d1a
  !     Computes the linear growing mode D1 in a Friedmann-Robertson-Walker
  !     universe. See Peebles LSSU sections 11 and 14.
  real(dp)::a,omega_m,omega_l,y12,y,eps,rombint
  
  eps=1.0d-6
  if(a .le. 0.0d0)then
     write(*,*)'a=',a
     stop
  end if
  y=omega_m*(1.d0/a-1.d0) + omega_l*(a*a-1.d0) + 1.d0
  if(y .lt. 0.0D0)then
     write(*,*)'y=',y
     stop
  end if
  y12=y**0.5d0
  d1a=y12/a*rombint(eps,a,eps,omega_m,omega_l)
  
  return
end function d1a
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function fpeebl(a,omega_m,omega_l)
  use amr_parameters, only: dp
  implicit none
  real(dp) :: fpeebl,a,omega_m,omega_l,rombint,fy
  !     Computes the growth factor f=d\log D1/d\log a.
  real(dp) :: fact,y,eps
  
  eps=1.0d-6
  y=omega_m*(1.d0/a-1.d0) + omega_l*(a*a-1.d0) + 1.d0
  fact=rombint(eps,a,eps,omega_m,omega_l)
  fpeebl=(omega_l*a*a-0.5d0*omega_m/a)/y - 1.d0 + a*fy(a,omega_m,omega_l)/fact
  return
end function fpeebl
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function rombint(a,b,tol,omega_m,omega_l)
  use amr_parameters, only: dp
  implicit none
  real(dp)::rombint
  !
  !     Rombint returns the integral from a to b of f(x)dx using Romberg 
  !     integration. The method converges provided that f(x) is continuous 
  !     in (a,b). The function f must be double precision and must be 
  !     declared external in the calling routine.  
  !     tol indicates the desired relative accuracy in the integral.
  !
  integer::maxiter=16,maxj=5
  real(dp),dimension(100):: g
  real(dp)::a,b,tol,fourj,omega_m,omega_l
  real(dp)::h,error,gmax,g0,g1,fy
  integer::nint,i,j,k,jmax
  
  h=0.5d0*(b-a)
  gmax=h*(fy(a,omega_m,omega_l)+fy(b,omega_m,omega_l))
  g(1)=gmax
  nint=1
  error=1.0d20
  i=0
10 i=i+1
  if(.not.  (i>maxiter.or.(i>5.and.abs(error)<tol)))then
     !	Calculate next trapezoidal rule approximation to integral.
     
     g0=0.0d0
     do k=1,nint
        g0=g0+fy(a+(k+k-1)*h,omega_m,omega_l)
     end do
     g0=0.5d0*g(1)+h*g0
     h=0.5d0*h
     nint=nint+nint
     jmax=min(i,maxj)
     fourj=1.0d0
     
     do j=1,jmax
        ! Use Richardson extrapolation.
        fourj=4.0d0*fourj
        g1=g0+(g0-g(j))/(fourj-1.0d0)
        g(j)=g0
        g0=g1
     enddo
     if (abs(g0).gt.tol) then
        error=1.0d0-gmax/g0
     else
        error=gmax
     end if
     gmax=g0
     g(jmax+1)=g0
     go to 10
  end if
  rombint=g0
  if (i>maxiter.and.abs(error)>tol) &
       &    write(*,*) 'Rombint failed to converge; integral, error=', &
       &    rombint,error
  return
end function rombint



