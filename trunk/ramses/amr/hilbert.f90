module hilbert

  ! Some parameters used in many of the modules routines
  integer(kind=4),parameter::twotondim = 2**NDIM

#if NDIM == 3
  integer, parameter :: big_shift = -60, left_shift = 4, right_shift = -1
#endif

#if NDIM == 2
  integer, parameter :: big_shift = -60, left_shift = 4, right_shift = -2
#endif
  
#if NDIM == 1
  integer, parameter :: big_shift = -62, left_shift = 2, right_shift = -1
#endif

  ! How many levels of refinement can be used per 64bit integer Hilbert key
  integer(kind=4), dimension(1:3), parameter :: levels_per_key = (/63, 31, 21/)
  integer(kind=4), dimension(1:3), parameter :: bits_per_int = (/63, 62, 63/)
  ! State diagrams taken from:
  ! J K Lawder, Using State Diagrams for Hilbert Curve Mappings
  ! http://www.dcs.bbk.ac.uk/TriStarp/pubs/JL2_00.pdf

  ! Usage of state diagrams: 
  
  ! - cartesian index: position of a cell inside 
  !   its oct in cartesian order (0 to 7)
  
  ! - current state: integer encoding the "orientation" 
  !   of the hilbert curve inside the oct
  
  ! - three_digit_diagram(cartesion index + current state * 8)  
  !   hilbert index of the cell given by cartesion index
  
  ! - next_state_diagram(cartesian index + current state * 8)  
  !   orientation of curve inside the cell given by cartesian index

#if NDIM==3
  integer(kind=8),parameter,dimension(0:95)::next_digits_diagram=(/&
       &   0, 1, 3, 2, 7, 6, 4, 5,&
       &   0, 7, 1, 6, 3, 4, 2, 5,&
       &   0, 3, 7, 4, 1, 2, 6, 5,&
       &   2, 3, 1, 0, 5, 4, 6, 7,&
       &   4, 3, 5, 2, 7, 0, 6, 1,&
       &   6, 5, 1, 2, 7, 4, 0, 3,&
       &   4, 7, 3, 0, 5, 6, 2, 1,&
       &   6, 7, 5, 4, 1, 0, 2, 3,&
       &   2, 5, 3, 4, 1, 6, 0, 7,&
       &   2, 1, 5, 6, 3, 0, 4, 7,&
       &   4, 5, 7, 6, 3, 2, 0, 1,&
       &   6, 1, 7, 0, 5, 2, 4, 3 /)

  integer(kind=4),parameter,dimension(0:95)::next_state_diagram=(/&
       &   1, 2, 3, 2, 4, 5, 3, 5,&
       &   2, 6, 0, 7, 8, 8, 0, 7,&
       &   0, 9,10, 9, 1, 1,11,11,&
       &   6, 0, 6,11, 9, 0, 9, 8,&
       &  11,11, 0, 7, 5, 9, 0, 7,&
       &   4, 4, 8, 8, 0, 6,10, 6,&
       &   5, 7, 5, 3, 1, 1,11,11,&
       &   6, 1, 6,10, 9, 4, 9,10,&
       &  10, 3, 1, 1,10, 3, 5, 9,&
       &   4, 4, 8, 8, 2, 7, 2, 3,&
       &   7, 2,11, 2, 7, 5, 8, 5,&
       &  10, 3, 2, 6,10, 3, 4, 4 /)
    
  integer(kind=4),parameter,dimension(0:95, 1:3) :: one_digit_diagram=reshape((/&
       & 0,  0,  0,  0,  1,  1,  1,  1,&  
       & 0,  0,  1,  1,  1,  1,  0,  0,&  
       & 0,  1,  1,  0,  0,  1,  1,  0,&  
       & 0,  0,  0,  0,  1,  1,  1,  1,&  
       & 1,  1,  0,  0,  0,  0,  1,  1,&  
       & 1,  0,  0,  1,  1,  0,  0,  1,&  
       & 0,  1,  1,  0,  0,  1,  1,  0,&  
       & 1,  1,  1,  1,  0,  0,  0,  0,&  
       & 1,  1,  0,  0,  0,  0,  1,  1,&  
       & 1,  0,  0,  1,  1,  0,  0,  1,&  
       & 1,  1,  1,  1,  0,  0,  0,  0,&  
       & 0,  0,  1,  1,  1,  1,  0,  0,&  

       & 0,  0,  1,  1,  1,  1,  0,  0,&  
       & 0,  1,  1,  0,  0,  1,  1,  0,&  
       & 0,  0,  0,  0,  1,  1,  1,  1,&  
       & 1,  1,  0,  0,  0,  0,  1,  1,&  
       & 0,  1,  1,  0,  0,  1,  1,  0,&  
       & 1,  1,  1,  1,  0,  0,  0,  0,&  
       & 1,  1,  1,  1,  0,  0,  0,  0,&  
       & 0,  0,  1,  1,  1,  1,  0,  0,&  
       & 1,  0,  0,  1,  1,  0,  0,  1,&  
       & 0,  0,  0,  0,  1,  1,  1,  1,&  
       & 1,  1,  0,  0,  0,  0,  1,  1,&  
       & 1,  0,  0,  1,  1,  0,  0,  1,&  
  
       & 0,  1,  1,  0,  0,  1,  1,  0,&  
       & 0,  0,  0,  0,  1,  1,  1,  1,&  
       & 0,  0,  1,  1,  1,  1,  0,  0,&  
       & 1,  0,  0,  1,  1,  0,  0,  1,&  
       & 1,  1,  1,  1,  0,  0,  0,  0,&  
       & 0,  0,  1,  1,  1,  1,  0,  0,&  
       & 1,  1,  0,  0,  0,  0,  1,  1,&  
       & 1,  0,  0,  1,  1,  0,  0,  1,&  
       & 0,  0,  0,  0,  1,  1,  1,  1,&  
       & 1,  1,  0,  0,  0,  0,  1,  1,&  
       & 0,  1,  1,  0,  0,  1,  1,  0,&  
       & 1,  1,  1,  1,  0,  0,  0,  0/), (/96, 3/))  
  
  ! Next state diagram for reverse (hilbert to cartesian key) conversion
  integer(kind=4),parameter,dimension(0:95)::next_state_diagram_reverse=(/&
       & 1,   2,   2,   3,   3,   5,   5,   4,&  
       & 2,   0,   0,   8,   8,   7,   7,   6,&  
       & 0,   1,   1,   9,   9,  11,  11,  10,&  
       &11,   6,   6,   0,   0,   9,   9,   8,&  
       & 9,   7,   7,  11,  11,   0,   0,   5,&  
       &10,   8,   8,   6,   6,   4,   4,   0,&  
       & 3,  11,  11,   5,   5,   1,   1,   7,&  
       & 4,   9,   9,  10,  10,   6,   6,   1,&  
       & 5,  10,  10,   1,   1,   3,   3,   9,&  
       & 7,   4,   4,   2,   2,   8,   8,   3,&  
       & 8,   5,   5,   7,   7,   2,   2,  11,&  
       & 6,   3,   3,   4,   4,  10,  10,   2/)  

#endif
#if NDIM==2
  ! State diagrams for 2D case
  integer(kind=4),parameter,dimension(0:15)::next_state_diagram=(/&
       & 1, 0, 2, 0, &
       & 0, 3, 1, 1, &
       & 2, 2, 0, 3, &
       & 3, 1, 3, 2/)
  
  integer(kind=8),parameter,dimension(0:15)::next_digits_diagram=(/&
       & 0, 1, 3, 2, &
       & 0, 3, 1, 2, &
       & 2, 1, 3, 0, &
       & 2, 3, 1, 0/)

  integer(kind=4),parameter,dimension(0:15)::next_state_diagram_reverse=(/&
       & 1, 0, 0, 2, &
       & 0, 1, 1, 3, &
       & 3, 2, 2, 0, &
       & 2, 3, 3, 1/)

  integer(kind=4),parameter,dimension(0:15, 1:2)::one_digit_diagram=reshape((/&
       & 0, 0, 1, 1, &
       & 0, 1 ,1, 0, &
       & 1, 0, 0, 1, &
       & 1, 1, 0, 0, &
       & 0, 1, 1, 0, &
       & 0, 0 ,1, 1, &
       & 1, 1, 0, 0, &
       & 1, 0, 0, 1/), (/16,2/))
#endif

#if NDIM==1
  ! State diagrams for 1D case - a bit silly, but...
  integer(kind=4),parameter,dimension(0:1)::next_state_diagram=(/0, 0/)
  
  integer(kind=8),parameter,dimension(0:1)::next_digits_diagram=(/0, 1/)

  integer(kind=4),parameter,dimension(0:1)::next_state_diagram_reverse=(/0, 0/)

  integer(kind=4),parameter,dimension(0:1, 1:1)::one_digit_diagram=reshape((/0, 1/), (/2,1/))
#endif

  
contains
  
  !================================================================
  !================================================================
  !================================================================
  !================================================================
  subroutine hilbert1d(x,order,npoint)
    use amr_parameters, ONLY: qdp
    implicit none
    integer     ,INTENT(IN)                     ::npoint
    integer     ,INTENT(IN) ,dimension(1:npoint)::x
    real(qdp),INTENT(OUT),dimension(1:npoint)::order

    integer::ip

    do ip=1,npoint
#ifdef QUADHILBERT
       order(ip)=real(x(ip),kind=16)
#else
       order(ip)=real(x(ip),kind=8)
#endif
    end do

  end subroutine hilbert1d
  !================================================================
  !================================================================
  !================================================================
  !================================================================
  subroutine hilbert2d_orig(x,y,order,bit_length,npoint)

    use amr_parameters, ONLY: qdp
    implicit none

    integer     ,INTENT(IN)                     ::bit_length,npoint
    integer     ,INTENT(IN) ,dimension(1:npoint)::x,y
    real(qdp),INTENT(OUT),dimension(1:npoint)::order

    logical,dimension(0:2*bit_length-1)::i_bit_mask 
    logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask
    integer,dimension(0:3,0:1,0:3)::state_diagram
    integer::i,ip,cstate,nstate,b0,b1,sdigit,hdigit

    if(bit_length>bit_size(bit_length))then
       write(*,*)'Maximum bit length=',bit_size(bit_length)
       write(*,*)'stop in hilbert2d'
       call clean_stop
    endif

    state_diagram = RESHAPE( (/ 1, 0, 2, 0, &
         & 0, 1, 3, 2, &
         & 0, 3, 1, 1, &
         & 0, 3, 1, 2, &
         & 2, 2, 0, 3, &
         & 2, 1, 3, 0, &
         & 3, 1, 3, 2, &
         & 2, 3, 1, 0  /), &
         & (/ 4, 2, 4 /) )

    do ip=1,npoint

       ! convert to binary
       do i=0,bit_length-1
          x_bit_mask(i)=btest(x(ip),i)
          y_bit_mask(i)=btest(y(ip),i)
       enddo

       ! interleave bits
       do i=0,bit_length-1
          i_bit_mask(2*i+1)=x_bit_mask(i)
          i_bit_mask(2*i  )=y_bit_mask(i)
       end do

       ! build Hilbert ordering using state diagram
       cstate=0
       do i=bit_length-1,0,-1
          b1=0 ; if(i_bit_mask(2*i+1))b1=1
          b0=0 ; if(i_bit_mask(2*i)  )b0=1
          sdigit=b1*2+b0
          nstate=state_diagram(sdigit,0,cstate)
          hdigit=state_diagram(sdigit,1,cstate)
          i_bit_mask(2*i+1)=btest(hdigit,1)
          i_bit_mask(2*i  )=btest(hdigit,0)
          cstate=nstate
       enddo

       ! save Hilbert key as double precision real
       order(ip)=0.
       do i=0,2*bit_length-1
          b0=0 ; if(i_bit_mask(i))b0=1
#ifdef QUADHILBERT
          order(ip)=order(ip)+real(b0,kind=16)*real(2,kind=16)**i
#else
          order(ip)=order(ip)+real(b0,kind=8)*real(2,kind=8)**i
#endif
       end do

    end do

  end subroutine hilbert2d_orig
  !================================================================
  !================================================================
  !================================================================
  !================================================================
  subroutine hilbert3d_orig(x,y,z,order,bit_length,npoint)
    use amr_parameters, ONLY: qdp
    implicit none

    integer     ,INTENT(IN)                     ::bit_length,npoint
    integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
    real(qdp),INTENT(OUT),dimension(1:npoint)::order

    logical,dimension(0:3*bit_length-1)::i_bit_mask
    logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
    integer,dimension(0:7,0:1,0:11)::state_diagram
    integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

    if(bit_length>bit_size(bit_length))then
       write(*,*)'Maximum bit length=',bit_size(bit_length)
       write(*,*)'stop in hilbert3d'
       call clean_stop
    endif

    state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
         &   0, 1, 3, 2, 7, 6, 4, 5,&
         &   2, 6, 0, 7, 8, 8, 0, 7,&
         &   0, 7, 1, 6, 3, 4, 2, 5,&
         &   0, 9,10, 9, 1, 1,11,11,&
         &   0, 3, 7, 4, 1, 2, 6, 5,&
         &   6, 0, 6,11, 9, 0, 9, 8,&
         &   2, 3, 1, 0, 5, 4, 6, 7,&
         &  11,11, 0, 7, 5, 9, 0, 7,&
         &   4, 3, 5, 2, 7, 0, 6, 1,&
         &   4, 4, 8, 8, 0, 6,10, 6,&
         &   6, 5, 1, 2, 7, 4, 0, 3,&
         &   5, 7, 5, 3, 1, 1,11,11,&
         &   4, 7, 3, 0, 5, 6, 2, 1,&
         &   6, 1, 6,10, 9, 4, 9,10,&
         &   6, 7, 5, 4, 1, 0, 2, 3,&
         &  10, 3, 1, 1,10, 3, 5, 9,&
         &   2, 5, 3, 4, 1, 6, 0, 7,&
         &   4, 4, 8, 8, 2, 7, 2, 3,&
         &   2, 1, 5, 6, 3, 0, 4, 7,&
         &   7, 2,11, 2, 7, 5, 8, 5,&
         &   4, 5, 7, 6, 3, 2, 0, 1,&
         &  10, 3, 2, 6,10, 3, 4, 4,&
         &   6, 1, 7, 0, 5, 2, 4, 3 /), &
         & (/8 ,2, 12 /) )

    do ip=1,npoint

       ! convert to binary
       do i=0,bit_length-1
          x_bit_mask(i)=btest(x(ip),i)
          y_bit_mask(i)=btest(y(ip),i)
          z_bit_mask(i)=btest(z(ip),i)
       enddo

       ! interleave bits
       do i=0,bit_length-1
          i_bit_mask(3*i+2)=x_bit_mask(i)
          i_bit_mask(3*i+1)=y_bit_mask(i)
          i_bit_mask(3*i  )=z_bit_mask(i)
       end do

       ! build Hilbert ordering using state diagram
       cstate=0
       do i=bit_length-1,0,-1
          b2=0 ; if(i_bit_mask(3*i+2))b2=1
          b1=0 ; if(i_bit_mask(3*i+1))b1=1
          b0=0 ; if(i_bit_mask(3*i  ))b0=1
          sdigit=b2*4+b1*2+b0
          nstate=state_diagram(sdigit,0,cstate)
          hdigit=state_diagram(sdigit,1,cstate)
          i_bit_mask(3*i+2)=btest(hdigit,2)
          i_bit_mask(3*i+1)=btest(hdigit,1)
          i_bit_mask(3*i  )=btest(hdigit,0)
          cstate=nstate
       enddo

       ! save Hilbert key as double precision real
       order(ip)=0.
       do i=0,3*bit_length-1
          b0=0 ; if(i_bit_mask(i))b0=1
#ifdef QUADHILBERT
          order(ip)=order(ip)+real(b0,kind=16)*real(2,kind=16)**i
#else
          order(ip)=order(ip)+real(b0,kind=8)*real(2,kind=8)**i
#endif
       end do

    end do

  end subroutine hilbert3d_orig
  !================================================================
  !================================================================
  !================================================================
  !================================================================
!   subroutine hilbert3d(x,y,z,order,bit_length,npoint)
!     use amr_parameters, ONLY: qdp,nvector
!     implicit none
! #ifndef WITHOUTMPI
!     include 'mpif.h'
! #endif
!     integer     ,INTENT(IN)                     ::bit_length,npoint
!     integer     ,INTENT(IN) ,dimension(1:nvector)::x,y,z
!     real(qdp),INTENT(OUT),dimension(1:nvector)::order

!     integer::i,ip,info
!     integer(kind=8)::testint
!     integer(kind=8),dimension(1:nvector)::hkey
!     integer(kind=4),dimension(1:nvector)::cstate,nstate,hdigit,sdigit,ind


!     if(bit_length>bit_size(bit_length))then
!        write(*,*)'Maximum bit length=',bit_size(bit_length)
!        write(*,*)'stop in hilbert3d'
!        call clean_stop
!     endif


! #ifdef QUADHILBERT
!     call MPI_ABORT(MPI_COMM_WORLD,1,info)
! #endif

!     ! build Hilbert ordering using state diagram
!     cstate=0
!     hkey=0
!     do i=bit_length-1,0,-1
!        do ip=1,npoint
!           hkey(ip)=hkey(ip)*longeight
!        end do

!        sdigit=0
!        do ip=1,npoint
!           if(btest(x(ip),i))sdigit(ip)=sdigit(ip)+four
!           if(btest(y(ip),i))sdigit(ip)=sdigit(ip)+two
!           if(btest(z(ip),i))sdigit(ip)=sdigit(ip)+one
!        end do

!        do ip=1,npoint
!           ind(ip)=cstate(ip)*eight+sdigit(ip)
!        end do

!        do ip=1,npoint
!           nstate(ip)=next_state_diagram3d(ind(ip))
!        end do

!        do ip=1,npoint
!           hkey(ip)=hkey(ip)+three_digit_diagram(ind(ip))
!        end do

!        do ip=1,npoint
!           cstate(ip)=nstate(ip)
!        end do
!     enddo

!     do ip=1,npoint
!        order(ip)=real(hkey(ip),kind=8)             
!     end do

!   end subroutine hilbert3d
  !================================================================
  !================================================================
  !================================================================
  !================================================================
!   subroutine hilbert3d_nonvec(x,y,z,order,bit_length,npoint)
!     use amr_parameters, ONLY: qdp,nvector
!     implicit none
! #ifndef WITHOUTMPI
!     include 'mpif.h'
! #endif
!     integer  , INTENT(IN)                       :: bit_length,npoint
!     integer  , INTENT(IN) ,dimension(1:nvector) :: x,y,z
!     real(qdp), INTENT(OUT), dimension(1:nvector) :: order



!     integer(kind=4)::i,ip,xx,yy,zz,info
!     integer(kind=8)::hkey
!     integer(kind=4)::cstate,nstate,sdigit,ind
!     if(bit_length>bit_size(bit_length))then
!        write(*,*)'Maximum bit length=',bit_size(bit_length)
!        write(*,*)'stop in hilbert3d'
!        call clean_stop
!     endif

! #ifdef QUADHILBERT
! #ifndef WITHOUTMPI
!     call MPI_ABORT(MPI_COMM_WORLD,1,info)
! #else
!     stop
! #endif
! #endif


!     do ip=1,npoint  
!        ! build Hilbert ordering using state diagram
!        xx=x(ip); yy=y(ip); zz=z(ip)
!        cstate=0
!        hkey=0
!        do i=bit_length-1,0,-1
!           hkey=hkey*longeight
!           sdigit=zero
!           if(btest(xx,i))sdigit=sdigit+four
!           if(btest(yy,i))sdigit=sdigit+two
!           if(btest(zz,i))sdigit=sdigit+one
!           ind=cstate*eight+sdigit
!           nstate=next_state_diagram3d(ind)
!           hkey=hkey+three_digit_diagram(ind)
!           cstate=nstate        
!        enddo
!        order(ip)=real(hkey,kind=8)             
!     end do

!   end subroutine hilbert3d_nonvec
  !================================================================
  !================================================================
  !================================================================
  !================================================================
  ! subroutine hilbert2d(ix, iy, hkey1, hkey0, cstate, &
  !      initial_level, final_level, npoint)
  !   use amr_parameters, only: qdp, nvector
  !   implicit none
  !   integer        , intent(in)                          :: initial_level, final_level, npoint
  !   integer(kind=8), intent(in),    dimension(1:nvector) :: ix, iy
  !   integer(kind=4), intent(inout), dimension(1:nvector) :: cstate
  !   integer(kind=8), intent(inout), dimension(1:nvector) :: hkey1, hkey0
    
  !   ! Compute nvector 2-integer hilbert keys from the cartesian keys ix, iy
    
  !   ! Local vars
  !   integer :: ibit, ip
  !   integer(kind=4),dimension(1:nvector) :: nstate, sdigit, ind
    
  !   ! if no keys present yet
  !   if (initial_level==0) then
  !      cstate(1:npoint)=0
  !      hkey0(1:npoint)=0; hkey1(1:npoint)=0
  !   end if

  !   do ibit=final_level-initial_level-1,0,-1

  !      ! Use only 62 out of the 64 bits for 2d case 
  !      ! (no unsigned in in fortran...)
  !      if (final_level > 31)then
  !         do ip=1,npoint
  !            hkey1(ip) = ISHFT(hkey1(ip),4)
  !            hkey1(ip) = ISHFT(hkey1(ip),-2)
  !         end do
  !         do ip=1,npoint
  !            hkey1(ip) = hkey1(ip) + ISHFT(hkey0(ip),-60)
  !         end do
  !      end if

  !      do ip=1,npoint
  !         hkey0(ip) = ISHFT(hkey0(ip),4)
  !         hkey0(ip) = ISHFT(hkey0(ip),-2)
  !      end do

  !      sdigit=0
  !      do ip=1,npoint
  !         if(btest(ix(ip),ibit)) sdigit(ip) = sdigit(ip)+two
  !         if(btest(iy(ip),ibit)) sdigit(ip) = sdigit(ip)+one
  !      end do

  !      do ip=1,npoint
  !         ind(ip) = cstate(ip)*four + sdigit(ip)
  !      end do

  !      do ip=1,npoint
  !         nstate(ip) = next_state_diagram(ind(ip))
  !      end do

  !      do ip=1,npoint
  !         hkey0(ip) = hkey0(ip) + two_digit_diagram(ind(ip))
  !      end do

  !      do ip=1,npoint
  !         cstate(ip) = nstate(ip)
  !      end do
  !   enddo
  ! end subroutine hilbert2d
  ! !================================================================
  ! !================================================================
  ! !================================================================
  ! !================================================================
  ! subroutine hilbert2d_reverse(ix, iy, hkey1, hkey0, key_level, npoint)
  !   use amr_parameters, only: nvector
  !   implicit none

  !   ! Inpu/Output variables:
  !   integer        , intent(in)                                :: key_level, npoint
  !   integer(kind=8), intent(out), dimension(1:nvector)         :: ix, iy
  !   integer(kind=8), intent(in),  dimension(1:nvector), target :: hkey0, hkey1

  !   ! Descripton:
  !   ! Compute nvector cartesian keys from the corresponding 2-integer hilbert keys.

  !   ! pointer to one of the two hkey arrays
  !   integer(kind=8), pointer              :: use_key(:)   
  !   integer                               :: ip, ibit1, ikey, ilevel
  !   integer(kind=4), dimension(1:nvector) :: cstate, nstate, ind
  !   integer(kind=8), dimension(1:nvector) :: sdigit

  !   ! Build the cartesian key using the state diagrams
  !   cstate=0
  !   ix=0; iy=0

  !   do ilevel=1,key_level
  !      ibit1 = (key_level-ilevel)*2
  !      ikey = ibit1/62

  !      ! use a pointer here to define which of the three integer keys 
  !      ! must be accessed.
  !      if (ikey==0) use_key => hkey0
  !      if (ikey==1) use_key => hkey1

  !      ibit1 = mod(ibit1,62)

  !      ! leftshift the cartesian keys by one position
  !      do ip=1,npoint
  !         ix(ip)=ISHFT(ix(ip),1)
  !         iy(ip)=ISHFT(iy(ip),1)
  !      end do

  !      ! read the next two bits from the hilbert key
  !      do ip=1,npoint
  !         sdigit(ip) = ibits(use_key(ip),ibit1,2)
  !      end do

  !      ! Compute lookup index in (flat) state diagrams
  !      do ip=1,npoint
  !         ind(ip) = cstate(ip)*four + sdigit(ip)
  !      end do

  !      ! save next state
  !      do ip=1,npoint
  !         nstate(ip) = next_state_diagram_reverse(ind(ip))
  !      end do

  !      ! add one integer key digit each
  !      do ip=1,npoint
  !         ix(ip) = ix(ip)+x_digit_diagram2d(ind(ip))
  !      end do
  !      do ip=1,npoint
  !         iy(ip) = iy(ip)+y_digit_diagram2d(ind(ip))
  !      end do
  !      do ip=1,npoint
  !         cstate(ip) = nstate(ip)
  !      end do
  !   enddo
  !================================================================
  !================================================================
  !================================================================
  !================================================================
  subroutine hilbert(ix, hkey, cstate, initial_level, final_level, npoint)
    use amr_parameters, only: nvector, int_pre, nhilbert, ndim
    implicit none
    integer, intent(in) :: initial_level, final_level, npoint
    integer(int_pre), dimension(:,:), intent(in)    :: ix
    integer(kind=4),  dimension(:  ), intent(inout) :: cstate
    integer(kind=8),  dimension(:,:), intent(inout) :: hkey
    
    ! Compute nvector 3-integer hilbert keys from the cartesian keys ix
        
    ! Local vars
    integer :: ibit, ip, add_digit, idim, ikey, nkey_local
    integer(kind=4), dimension(1:nvector) :: nstate, sdigit, ind

    ! if no keys present yet
    if (initial_level == 0) then
       cstate(1:npoint) = 0
       hkey(1:npoint, 1:nhilbert) = 0
    end if

    nkey_local = ceiling(1.d0 * final_level / levels_per_key(ndim))
    
    do ibit = final_level - initial_level - 1, 0, -1
       
       do ikey = nkey_local, 1, -1
          do ip = 1, npoint
             hkey(ip, ikey) = ISHFT(hkey(ip, ikey), left_shift)
             hkey(ip, ikey) = ISHFT(hkey(ip, ikey), right_shift)
          end do
          if (ikey > 1)then
             do ip=1,npoint
                hkey(ip, ikey) = hkey(ip, ikey) + ISHFT(hkey(ip, ikey - 1), big_shift)
             end do
          end if
       end do


       sdigit=0
       do idim = 1, ndim
          add_digit = 2 ** (ndim - idim)
          do ip = 1, npoint
             if(btest(ix(ip, idim),ibit)) sdigit(ip) = sdigit(ip) + add_digit
          end do
       end do

       do ip=1,npoint
          ind(ip) = cstate(ip) * twotondim + sdigit(ip)
       end do

       do ip=1,npoint
          nstate(ip) = next_state_diagram(ind(ip))
       end do

       do ip=1,npoint
          hkey(ip, 1) = hkey(ip, 1) + next_digits_diagram(ind(ip))
       end do

       do ip=1,npoint
          cstate(ip) = nstate(ip)
       end do
    enddo
  end subroutine hilbert
  !================================================================
  !================================================================
  !================================================================
  !================================================================
  subroutine hilbert_reverse(ix, hkey, key_level, npoint)
    use amr_parameters, only: nvector, int_pre, ndim
    implicit none

    ! Inpu/Output variables:
    integer        , intent(in)                  :: key_level, npoint
    integer(int_pre), intent(out), dimension(:,:):: ix
    integer(kind=8), intent(in),  dimension(:,:) :: hkey

    ! Descripton:
    ! Compute nvector cartesian keys from the corresponding 3-integer hilbert keys.

    ! pointer to one of the three hkey arrays
    integer                               :: ip, ibit1, ikey, ilevel, idim
    integer(kind=4), dimension(1:nvector) :: cstate, nstate, ind
    integer(kind=8), dimension(1:nvector) :: sdigit

    ! Build the cartesian key using the state diagrams
    cstate = 0; ix = 0

    do ilevel = 1, key_level
       ibit1 = (key_level-ilevel) * ndim
       ikey = ibit1 / bits_per_int(ndim) + 1
       ibit1 = mod(ibit1, bits_per_int(ndim))

       ! leftshift the cartesian keys by one position
       do idim = 1, ndim
          do ip = 1, npoint
             ix(ip, idim) = ISHFT(ix(ip, idim),1)
          end do
       end do

       ! read the next ndim bits from the hilbert key
       do ip = 1, npoint
          sdigit(ip) = ibits(hkey(ip, ikey), ibit1, ndim)
       end do

       ! Compute lookup index in state diagrams
       do ip=1,npoint
          ind(ip) = cstate(ip) * twotondim + sdigit(ip)
       end do

       ! save next state
       do ip=1,npoint
          nstate(ip) = next_state_diagram_reverse(ind(ip))
       end do

       ! add one integer key digit each
       do idim = 1, ndim
          do ip=1,npoint
             ix(ip, idim) = ix(ip, idim) + one_digit_diagram(ind(ip), idim)
          end do
       end do

       do ip=1,npoint
          cstate(ip) = nstate(ip)
       end do
    enddo

  end subroutine hilbert_reverse

  !================================================================
  !================================================================
  !================================================================
  !================================================================
  subroutine hilbert_for_particle(offset, nparts, initial_level, final_level)
    use amr_parameters, only: nvector, boxlen, dp, int_pre, ndim, nhilbert
    use pm_commons,     only: part_hkey, current_state, xp
    implicit none

    integer, intent(in) :: initial_level, final_level
    integer, intent(in) :: offset, nparts

    ! Description:
    ! This subroutine computes 3D hilbert keys for particles 
    ! It assumes that the particles are stored as contiguous 
    ! arrays in memory and that positions, 3-integer hilbert keys
    ! and next_state are allocated as particle-based quantities.

    ! Iputs: 
    ! - Starting offset in particle arrays and number of particles
    !   to process (np)
    ! - Level of already computed hilbert key 
    ! - Desired level of hilbert key on exit

    ! Example: 
    ! call hilbert3d_for_particle(0, npart_levelmin, nlevelmax-1, nlevelmax) 
    ! will compute the last 3 bits of the hilbert key for the
    ! levelmin particles (resulting in a total of 3 * nlevelmax bits)

    ! Local variables
    integer :: ibit, ip, ind_part, idim, np, ioft
    integer(int_pre), dimension(1:nvector, 1:ndim) :: ix
    real(dp) :: ckey_factor
    
    ! Compute particle position to cartesian key factor
    ckey_factor = 2.0**final_level / dble(boxlen)

    do ioft = offset, offset + nparts - 1, nvector
       np = min(nvector, offset + nparts - ioft)

       ! compute cartesian keys
       do idim = 1, ndim
          do ip = 1, np
             ix(ip, idim) = floor(xp(ioft + ip, idim) * ckey_factor, kind=int_pre)
          end do
       end do

       ! Passing in array slices is ok (no copying) if the dummy argument has assumed shape and the interface is explicit!
       call hilbert_nd(ix, part_hkey(ioft + 1: ioft + np, 1:nhilbert), current_state(ioft + 1: ioft + np), initial_level, final_level, np)
    end do

  end subroutine hilbert_for_particle


  function ge_keys(key_a, key_b)
    implicit none
    logical::ge_keys

#if NHILBERT == 1
    integer(kind=8), intent(in), dimension(0:0) :: key_a, key_b

    ge_keys = (key_a(0) >= key_b(0))
#endif
  
#if NHILBERT == 2
    integer(kind=8), intent(in), dimension (0:1):: key_a, key_b
    if     (key_a(1) > key_b(1)) then
       ge_keys  =  .true.
    elseif (key_a(1) < key_b(1)) then
       ge_keys = .false.
    elseif (key_a(0) > key_b(0)) then
       ge_keys = .true.
    elseif (key_a(0) < key_b(0)) then
       ge_keys = .false.
    else
       ge_keys = .true.
    end if
#endif
    
#if NHILBERT == 3  
    integer(kind=8), intent(in), dimension (0:2):: key_a, key_b

    if     (key_a(2) > key_b(2)) then
       ge_keys = .true.
    elseif (key_a(2) < key_b(2)) then
       ge_keys = .false.
    elseif (key_a(1) > key_b(1)) then
       ge_keys = .true.
    elseif (key_a(1) < key_b(1)) then
       ge_keys = .false.
    elseif (key_a(0) > key_b(0)) then
       ge_keys = .true.
    elseif (key_a(0) < key_b(0)) then
       ge_keys = .false.
    else
       ge_keys = .true.
    end if
#endif

  end function ge_keys
 
  function gt_keys(key_a, key_b)
    implicit none
    logical::gt_keys

    integer(kind=8), intent(in), dimension (:):: key_a, key_b
#if NHILBERT == 1
    gt_keys = (key_a(1) > key_b(1))
#endif

#if NHILBERT == 2
    if     (key_a(2) > key_b(2)) then
       gt_keys  =  .true.
    elseif (key_a(2) < key_b(2)) then
       gt_keys = .false.
    elseif (key_a(1) > key_b(1)) then
       gt_keys = .true.
    elseif (key_a(1) < key_b(1)) then
       gt_keys = .false.
    else
       gt_keys = .false.
    end if
#endif

#if NHILBERT == 3
    if     (key_a(3) > key_b(3)) then
       gt_keys = .true.
    elseif (key_a(3) < key_b(3)) then
       gt_keys = .false.
    elseif (key_a(2) > key_b(2)) then
       gt_keys = .true.
    elseif (key_a(2) < key_b(2)) then
       gt_keys = .false.
    elseif (key_a(1) > key_b(1)) then
       gt_keys = .true.
    elseif (key_a(1) < key_b(1)) then
       gt_keys = .false.
    else
       gt_keys = .false.
    end if
#endif
  end function gt_keys


  
end module hilbert
