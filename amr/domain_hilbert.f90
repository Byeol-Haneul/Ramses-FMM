module domain_m
  use hilbert
  implicit none
  private
  public :: domain_t, get_domain

  type domain_t
     integer :: myid
     integer :: ncpu
     integer :: last_domain
!     integer :: n                                                 ! number of points in domain
!     integer :: overload                                          ! number of domains per cpu rank
     integer(kind=8), dimension(:,:), allocatable :: b            ! bound key
!     integer,         dimension(:),   allocatable :: r2d          ! rank to domain
!     integer,         dimension(:),   allocatable :: d2r          ! domain to rank
   contains
     procedure :: copy     => copy_domain
     procedure :: create   => create_domain
     procedure :: destroy  => destroy_domain
     procedure :: in_rank  => in_rank
     procedure :: get_rank => get_rank
     procedure, nopass :: pos2key  => hilbert_key
  end type domain_t
  
contains
 
  !================================================================
  !================================================================
  !================================================================
  !================================================================

  pure function in_rank(domain,key)
    use amr_parameters, only : nhilbert
    implicit none
    class(domain_t),                       intent(in) :: domain
    integer(kind=8), dimension(1:nhilbert), intent(in) :: key
    logical :: in_rank
    integer :: isub, myid, ncpu
    !
    in_rank=.false.
    myid = domain%myid
    ncpu = domain%ncpu
    in_rank = ge_keys(key,domain%b(1:nhilbert,myid-1)).and. &
         &    gt_keys(domain%b(1:nhilbert,myid),key)
!    do isub = 0,domain%overload-1
!       in_rank = ge_keys(key,domain%b(1:nhilbert,domain%r2d(myid+ncpu*isub)-1)).and. &
!                 gt_keys(domain%b(1:nhilbert,domain%r2d(myid+ncpu*isub)),key)
!    end do
  end function in_rank

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  function get_rank(domain,key) result(rank)
    use amr_parameters, only : nhilbert
    implicit none
    class(domain_t),                        intent(inout) :: domain
    integer(kind=8), dimension(1:nhilbert), intent(in) :: key
    integer :: rank, idom, ncpu
    !
    ncpu = domain%ncpu
!    idom = modulo(domain%d2r(idom)-1,ncpu)+1
!    rank = modulo(domain%d2r(idom)-1,ncpu)+1
!    rank = get_domain_lin(domain,key)
    rank = get_domain_bis(domain,key)
  end function get_rank

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  subroutine copy_domain(domain,source)
    implicit none
    class(domain_t)             :: domain
    class(domain_t), intent(in) :: source
    integer::myid
    !
    myid = domain%myid
    if (source%ncpu==0) then
       write(*,'(a,i7.7,a)') 'CPU=',myid,' ERROR copy_domain: Source domain not allocated.'
       stop
    endif
    !
    if (domain%ncpu .ne. source%ncpu) call create_domain(domain,source%myid,source%ncpu,source%ncpu)
    domain%b   = source%b
!    domain%r2d = source%r2d
!    domain%d2r = source%d2r
  end subroutine copy_domain

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  subroutine create_domain(domain,myid,ncpu,n)
    use amr_parameters, only : nhilbert
    implicit none
    class(domain_t) :: domain
    integer, intent(in) :: myid, ncpu, n
    integer :: i
    !
    if (domain%ncpu > 0) call destroy_domain(domain)
    domain%myid = myid
    domain%ncpu = ncpu
!    domain%n = n
!    domain%overload = n / ncpu
!    if (domain%ncpu .ne. domain%overload*ncpu) then
!       write(*,'(a,i7.7,a)') 'CPU=',myid,' ERROR create_domain: Number of domains not divisible by number of ranks.'
!       write(*,'(a,i7.7,a,i7,a,i7,a,f10.2)') 'CPU=',myid,' ERROR create_domain: Ndomain=',n,' Ncpu=', ncpu, ' overload=', real(n)/ncpu
!       stop
!    endif
    allocate(domain%b(1:nhilbert,0:domain%ncpu))
!    allocate(domain%d2r(1:domain%n))
!    allocate(domain%r2d(1:domain%n))
    ! Set default values for the domain:
    !   - All bound keys are zero
    !   - Ranks and domains map as the identity
    domain%b = 0
!    do i=1,domain%n
!       domain%d2r(i) = modulo(i-1,ncpu) + 1
!       domain%r2d(i) = i
!    enddo
  end subroutine create_domain

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  subroutine destroy_domain(domain)
    implicit none
    class(domain_t) :: domain
    !
!    domain%n = 0
!    domain%overload = 0
    domain%myid = 0
    domain%ncpu = 0
    if (allocated(domain%b))   deallocate(domain%b)
!    if (allocated(domain%r2d)) deallocate(domain%r2d)
!    if (allocated(domain%d2r)) deallocate(domain%d2r)
  end subroutine destroy_domain

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  function get_domain(domain,key) result(dom)
    use amr_parameters, only : nhilbert
    implicit none
    class(domain_t),                        intent(in) :: domain
    integer(kind=8), dimension(1:nhilbert), intent(in) :: key
    integer            :: dom
    !
    integer            :: idom, myid
    integer,      save :: last_domain=1, ncall=0
    logical,      save :: do_linear=.true.
    real(kind=8), save :: tlin,thunt
    real(kind=8)       :: t0=0.
    real(kind=8), external :: wallclock
    !
    myid = domain%myid
    ! The first 100 calls we use for timing
    if (ncall<=100) then
       if (ncall==0) then
          tlin=0; thunt=0
          do_linear = .true.
          ncall = ncall + 1
       end if
       t0 = wallclock()
    endif
    !
    if (do_linear) then                                                          ! do dumb linear search
       dom=1                                                                     ! default value
       do idom=1,domain%ncpu-1
          if (ge_keys(key,domain%b(1:nhilbert,idom))) then
             dom = dom + 1                                                       ! have to check lower bound
          else
             exit
          endif
       end do
    else                                                                         ! use hunt to predictively find the right domain
       dom = last_domain
       call hunt(domain%b,key,dom,domain%ncpu)
       last_domain = dom                                                         ! hopefully a close guess
       dom = dom + 1                                                             ! hunt returns lower bound, we need upper
    end if

    ! The first 100 calls we use for timing
    if (ncall<=100) then
       if (do_linear) then
          tlin  = tlin  + (wallclock() - t0)
       else
          thunt = thunt + (wallclock() - t0)
       end if
       if (ncall==100) then
          if (do_linear) then
             ncall=1; do_linear = .false.                                         ! first round done, now time hunting
          else
             do_linear = tlin < thunt                                             ! second round done, use whatever is fastest
             if (myid==1) then
                print *, ' Time to do linear search :', real(tlin)
                print *, ' Time to do hunt search   :', real(thunt)
                if (do_linear) then
                   print *, ' Using linear search'
                else
                   print *, ' Using hunt search'
                end if
             end if
          end if
       end if
    end if
    !
    ncall = modulo(ncall + 1,20000000)                                            ! Use 0.1% of the calls for timing
    !
  end function get_domain

  function get_domain_lin(domain,key) result(dom)
    use amr_parameters, only : nhilbert
    implicit none
    class(domain_t),                        intent(in) :: domain
    integer(kind=8), dimension(1:nhilbert), intent(in) :: key
    integer            :: dom
    !
    integer            :: idom, myid
    !
    myid = domain%myid
    !
    dom=1                                                                     ! default value
    do idom=1,domain%ncpu-1
       if (ge_keys(key,domain%b(1:nhilbert,idom))) then
          dom = dom + 1                                                       ! have to check lower bound
       else
          exit
       endif
    end do
    !
  end function get_domain_lin

  function get_domain_bis(domain,key) result(dom)
    use amr_parameters, only : nhilbert
    implicit none
    class(domain_t),                     intent(inout) :: domain
    integer(kind=8), dimension(1:nhilbert), intent(in) :: key
    integer            :: dom
    !
    integer            :: idom, myid
    !
    myid = domain%myid
    !
    dom = domain%last_domain
    call hunt(domain%b,key,dom,domain%ncpu)
    domain%last_domain = dom                                         ! hopefully a close guess
    dom = dom + 1                                                    ! hunt returns lower bound, we need upper
    !
  end function get_domain_bis

  !================================================================
  !================================================================
  ! Given an array xx(1:nhilbert,0:n), and given a key x(1:nhilbert), return a value jlo such that
  ! xx(:,jlo) <= x < xx(:,jlo+1). xx must be monotonic increasing. jlo on input is taken
  ! as the initial guess for jlo on output.
  !================================================================
  !================================================================

  subroutine hunt(xx,x,jlo,n)
    use amr_parameters, only : nhilbert
    implicit none
    integer :: jlo
    integer, intent(in) :: n
    integer(kind=8), dimension(1:nhilbert),     intent(in) :: x
    integer(kind=8), dimension(1:nhilbert,0:n), intent(in) :: xx
    integer :: inc,jhi,jm
  !...............................................................................
    if (jlo < 0 .or. jlo > n) then                                                ! Input guess not useful. Go immediately to bisection.
      jlo=-1
      jhi=n+1
    else
      inc=1                                                                       ! Set the hunting increment.
      if (ge_keys(x,xx(:,jlo))) then                                              ! Hunt up:
        do
          jhi=jlo+inc
          if (jhi > n) then                                                       ! Done hunting, since off end of table
            jhi= n+1
            exit
          else
            if (gt_keys(xx(:,jhi),x)) exit
            jlo=jhi                                                               ! Not done hunting,
            inc=inc+inc                                                           ! so double the increment
          end if
        end do                                                                    ! and try again.
      else                                                                        ! Hunt down:
        jhi=jlo
        do
          jlo=jhi-inc
          if (jlo < 0) then                                                       ! Done hunting, since off end of table.
            jlo=-1
            exit
          else
            if (ge_keys(x,xx(:,jlo))) exit
            jhi=jlo                                                               ! Not done hunting,
            inc=inc+inc                                                           ! so double the increment
          end if
        end do                                                                    ! and try again.
      end if
    end if                                                                        ! Done hunting, value bracketed.
    do                                                                            ! Hunt is done, so begin the final bisection phase:
      if (jhi-jlo <= 1) then
        if (jlo>n .and. eq_keys(x,xx(:,n))) jlo=n
        if (jlo<0 .and. eq_keys(x,xx(:,0))) jlo=0
        exit
      else
        jm=(jhi+jlo)/2
        if (ge_keys(x,xx(:,jm))) then
          jlo=jm
        else
          jhi=jm
        end if
      end if
    end do
  end subroutine hunt

end module domain_m
