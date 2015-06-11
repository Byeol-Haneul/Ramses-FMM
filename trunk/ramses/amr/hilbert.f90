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
subroutine hilbert2d(x,y,order,bit_length,npoint)

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

end subroutine hilbert2d
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
subroutine hilbert3d_andreas(x,y,z,order,bit_length,npoint)
  use amr_parameters, ONLY: qdp,nvector
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer     ,INTENT(IN) ,dimension(1:nvector)::x,y,z
  real(qdp),INTENT(OUT),dimension(1:nvector)::order

  integer(kind=4),parameter,dimension(0:95)::state_diagram0=(/1, 2, 3, 2, 4, 5, 3, 5,&
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
  integer(kind=8),parameter,dimension(0:95)::state_diagram1=(/0, 1, 3, 2, 7, 6, 4, 5,&
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

  ! state diagrams:  state_diagram0(current cube/digit,  input orientation/state) gives new orientation/state
  !                  state_diagram1(current cube/digit,  input orientation/state) gives new cube/orientation



  integer::i,ip,info
  integer(kind=8)::testint
  integer(kind=8),dimension(1:nvector)::hkey
  integer(kind=4),dimension(1:nvector)::cstate,nstate,hdigit,sdigit,ind
  integer(kind=4),parameter::zero=0
  integer(kind=4),parameter::one=1
  integer(kind=4),parameter::two=2
  integer(kind=4),parameter::four=4
  integer(kind=4),parameter::eight=8
  integer(kind=8),parameter::longeight=8


  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert3d'
     call clean_stop
  endif


#ifdef QUADHILBERT
        call MPI_ABORT(MPI_COMM_WORLD,1,info)
#endif

  ! build Hilbert ordering using state diagram
  cstate=0
  hkey=0
  do i=bit_length-1,0,-1
     do ip=1,npoint
        hkey(ip)=hkey(ip)*longeight
     end do

     sdigit=0
     do ip=1,npoint
        if(btest(x(ip),i))sdigit(ip)=sdigit(ip)+four
        if(btest(y(ip),i))sdigit(ip)=sdigit(ip)+two
        if(btest(z(ip),i))sdigit(ip)=sdigit(ip)+one
     end do

     do ip=1,npoint
        ind(ip)=cstate(ip)*eight+sdigit(ip)
     end do

     do ip=1,npoint
        nstate(ip)=state_diagram0(ind(ip))
     end do

     do ip=1,npoint
        hkey(ip)=hkey(ip)+state_diagram1(ind(ip))
     end do

     do ip=1,npoint
        cstate(ip)=nstate(ip)
     end do
  enddo
  
  do ip=1,npoint
     order(ip)=real(hkey(ip),kind=8)             
  end do

end subroutine hilbert3d_andreas
!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert3d_nonvec(x,y,z,order,bit_length,npoint)
  use amr_parameters, ONLY: qdp,nvector
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer     ,INTENT(IN) ,dimension(1:nvector)::x,y,z
  real(qdp),INTENT(OUT),dimension(1:nvector)::order

  integer(kind=4),parameter,dimension(0:95)::state_diagram0=(/1, 2, 3, 2, 4, 5, 3, 5,&
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
  integer(kind=8),parameter,dimension(0:95)::state_diagram1=(/0, 1, 3, 2, 7, 6, 4, 5,&
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

  ! state diagrams:  state_diagram0(current cube/digit,  input orientation/state) gives new orientation/state
  !                  state_diagram1(current cube/digit,  input orientation/state) gives new cube/orientation


  integer(kind=4),parameter::zero=0
  integer(kind=4),parameter::one=1
  integer(kind=4),parameter::two=2
  integer(kind=4),parameter::four=4
  integer(kind=4),parameter::eight=8
  integer(kind=8),parameter::longeight=8
  integer(kind=4)::i,ip,xx,yy,zz,info
  integer(kind=8)::hkey
  integer(kind=4)::cstate,nstate,sdigit,hdigit,ind
  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert3d'
     call clean_stop
  endif

  
  
#ifdef QUADHILBERT
#ifndef WITHOUTMPI
        call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
        stop
#endif
#endif
  
  
  do ip=1,npoint  
     ! build Hilbert ordering using state diagram
     xx=x(ip); yy=y(ip); zz=z(ip)
     cstate=0
     hkey=0
     do i=bit_length-1,0,-1
        hkey=hkey*longeight
        sdigit=zero
        if(btest(xx,i))sdigit=sdigit+four
        if(btest(yy,i))sdigit=sdigit+two
        if(btest(zz,i))sdigit=sdigit+one
        ind=cstate*eight+sdigit
        nstate=state_diagram0(ind)
        hkey=hkey+state_diagram1(ind)
        cstate=nstate        
     enddo
     order(ip)=real(hkey,kind=8)             
  end do
  
end subroutine hilbert3d_nonvec
!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert3d_multiint(x,y,z,hkey2,hkey1,hkey0,bit_length,npoint)
  use amr_parameters, ONLY: qdp,nvector
  implicit none

  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer(kind=8),INTENT(IN) ,dimension(1:nvector)::x,y,z
!  real(qdp),INTENT(OUT),dimension(1:nvector)::order
  integer(kind=8),INTENT(OUT),dimension(1:nvector)::hkey2,hkey1,hkey0

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
  integer(kind=8),parameter,dimension(0:95)::three_digit_diagram=(/&
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

  ! state diagrams:  state_diagram0(current cube/digit,  input orientation/state) gives new orientation/state
  !                  state_diagram1(current cube/digit,  input orientation/state) gives new cube/orientation



  integer::i,ip
!#ifdef HILBERT>1
!#ifdef QUADHILBERT 
!  integer(kind=8),dimension(1:nvector)::hkey2
!#endif
! #ifdef HILBERT>2
!   integer(kind=8),dimension(1:nvector)::hkey2
! #endif
  integer(kind=4),dimension(1:nvector)::cstate,nstate,hdigit,sdigit,ind
  integer(kind=4),parameter::zero=0
  integer(kind=4),parameter::one=1
  integer(kind=4),parameter::two=2
  integer(kind=4),parameter::four=4
  integer(kind=4),parameter::eight=8
  integer(kind=8),parameter::longeight=8


  ! build Hilbert ordering using state diagram
  cstate=0
  hkey0=0; hkey1=0; hkey2=0
!#ifdef HILBERT>1
!#ifdef QUADHILBERT
!  hkey1=0
!#endif
! #ifdef HILBERT>2
!   hkey1=0
! #endif

  do i=bit_length-1,0,-1

     if (bit_length>42)then
        do ip=1,npoint
           hkey2(ip)=ISHFT(hkey2(ip),3)
        end do
        do ip=1,npoint
           hkey2(ip)=hkey2(ip)+ISHFT(hkey1(ip),-60)
        end do
     end if

     if (bit_length>21)then
        do ip=1,npoint
           hkey1(ip)=ISHFT(hkey1(ip),4)
           hkey1(ip)=ISHFT(hkey1(ip),-1)
        end do
        do ip=1,npoint
           hkey1(ip)=hkey1(ip)+ISHFT(hkey0(ip),-60)
        end do
     end if

     do ip=1,npoint
        hkey0(ip)=ISHFT(hkey0(ip),4)
        hkey0(ip)=ISHFT(hkey0(ip),-1)
     end do

     sdigit=0
     do ip=1,npoint
        if(btest(x(ip),i))sdigit(ip)=sdigit(ip)+four
        if(btest(y(ip),i))sdigit(ip)=sdigit(ip)+two
        if(btest(z(ip),i))sdigit(ip)=sdigit(ip)+one
     end do

     do ip=1,npoint
        ind(ip)=cstate(ip)*eight+sdigit(ip)
     end do

     do ip=1,npoint
        nstate(ip)=next_state_diagram(ind(ip))
     end do

     do ip=1,npoint
        hkey0(ip)=hkey0(ip)+three_digit_diagram(ind(ip))
     end do

     do ip=1,npoint
        cstate(ip)=nstate(ip)
     end do
  enddo
  
! #ifdef QUADHILBERT
!   do ip=1,npoint
!      order(ip)=real(hkey0(ip),kind=16)+real(hkey1(ip),kind=16)*real(2,kind=16)**63
!   end do
! #else
!   do ip=1,npoint
!      order(ip)=real(hkey0(ip),kind=8)             
!   end do  
! #endif

end subroutine hilbert3d_multiint
!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert3d_multiint_reverse(x,y,z,hkey2,hkey1,hkey0,bit_length,npoint)
  use amr_parameters, ONLY: nvector
  implicit none

  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer(kind=8),INTENT(OUT) ,dimension(1:nvector)::x,y,z
  integer(kind=8),dimension(1:nvector)::hkey2,hkey1,hkey0


  
  integer(kind=4),parameter,dimension(0:95)::x_digit_diagram=(/&
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
       & 0,  0,  1,  1,  1,  1,  0,  0/)  
  
  integer(kind=4),parameter,dimension(0:95)::y_digit_diagram=(/&
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
       & 1,  0,  0,  1,  1,  0,  0,  1/)  
  
  integer(kind=4),parameter,dimension(0:95)::z_digit_diagram=(/&
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
       & 1,  1,  1,  1,  0,  0,  0,  0/)  

  integer(kind=4),parameter,dimension(0:95)::next_state_diagram=(/&
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




  integer::i,ip,leading_zeroes
  integer(kind=4),dimension(1:nvector)::cstate,nstate,hdigit,ind
  integer(kind=8),dimension(1:nvector)::sdigit
  integer(kind=4),parameter::zero=0
  integer(kind=4),parameter::one=1
  integer(kind=4),parameter::two=2
  integer(kind=4),parameter::four=4
  integer(kind=4),parameter::eight=8
  integer(kind=8),parameter::longeight=8

  ! compute the leading zeroes in the integer key
  leading_zeroes=mod(3*(63-bit_length),63)
  
  ! shift by leading_zeroes bits to the left
  do ip=1,npoint
     hkey2(ip)=ISHFT(hkey2(ip),leading_zeroes)
  end do
  do ip=1,npoint
     hkey2(ip)=hkey2(ip)+ISHFT(hkey1(ip),leading_zeroes-63)
  end do
  do ip=1,npoint
     hkey1(ip)=ISHFT(hkey1(ip),leading_zeroes+1)
     hkey1(ip)=ISHFT(hkey1(ip),-1)
  end do
  do ip=1,npoint
     hkey1(ip)=hkey1(ip)+ISHFT(hkey0(ip),leading_zeroes-63)
  end do  
  do ip=1,npoint
     hkey0(ip)=ISHFT(hkey0(ip),leading_zeroes+1)
     hkey0(ip)=ISHFT(hkey0(ip),-1)
  end do

  


  ! build the cartesian key using the state diagrams
  cstate=0
  x=0; y=0; z=0
  do i=bit_length-1,0,-1     
     ! leftshift the cartesian keys by one position
     do ip=1,npoint
        x(ip)=ISHFT(x(ip),1)
        y(ip)=ISHFT(y(ip),1)
        z(ip)=ISHFT(z(ip),1)
     end do
     
     ! compute index for state diagrams from the 3 leading bits in the 
     ! relevant key integer
     if (bit_length>42)then
        do ip=1,npoint
           sdigit(ip)=ISHFT(hkey2(ip),-60)
        end do
     else if(bit_length>21)then
        do ip=1,npoint
           sdigit(ip)=ISHFT(hkey1(ip),-60)
        end do
     else
        do ip=1,npoint
           sdigit(ip)=ISHFT(hkey0(ip),-60)
        end do
     end if
     
     do ip=1,npoint
        ind(ip)=cstate(ip)*eight+sdigit(ip)
     end do

     ! save next state
     do ip=1,npoint
        nstate(ip)=next_state_diagram(ind(ip))
     end do
     
     ! add one integer key digit each
     do ip=1,npoint
        x(ip)=x(ip)+x_digit_diagram(ind(ip))
     end do
     do ip=1,npoint
        y(ip)=y(ip)+y_digit_diagram(ind(ip))
     end do
     do ip=1,npoint
        z(ip)=z(ip)+z_digit_diagram(ind(ip))
     end do

     do ip=1,npoint
        cstate(ip)=nstate(ip)
     end do


     ! shift the key by 3 bits to the left
     if (bit_length>42)then
        do ip=1,npoint
           hkey2(ip)=ISHFT(hkey2(ip),4)
           hkey2(ip)=ISHFT(hkey2(ip),-1)
        end do
        do ip=1,npoint
           hkey2(ip)=hkey2(ip)+ISHFT(hkey1(ip),-60)
        end do
     end if

     if (bit_length>21)then
        do ip=1,npoint
           hkey1(ip)=ISHFT(hkey1(ip),4)
           hkey1(ip)=ISHFT(hkey1(ip),-1)
        end do
        do ip=1,npoint
           hkey1(ip)=hkey1(ip)+ISHFT(hkey0(ip),-60)
        end do
     end if


     do ip=1,npoint
        hkey0(ip)=ISHFT(hkey0(ip),4)
        hkey0(ip)=ISHFT(hkey0(ip),-1)
     end do
     
  
  enddo
  
end subroutine hilbert3d_multiint_reverse







!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert3d_for_particle(offset, np, &
                                  initial_level, final_level)

  use amr_parameters, only: nvector, boxlen, dp
  use amr_commons,    only: myid
  use pm_commons,     only: part_hkey, current_state, xp_andreas
  implicit none

  integer, intent(in) :: initial_level, final_level
  integer, intent(in) :: offset, np

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
  
  
  ! State diagrams (lookup tables for hilbert key computation)

  integer(kind=8),parameter,dimension(0:95)::three_digit_diagram=(/&
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

  ! Usage of state diagrams: 

  ! - cartesian index: position of a cell inside 
  !   its oct in cartesian order (0 to 7)

  ! - current state: integer encoding the "orientation" 
  !   of the hilbert curve inside the oct

  ! - three_digit_diagram(cartesion index + current state * 8)  
  !   hilbert index of the cell given by cartesion index

  ! - next_state_diagram(cartesian index + current state * 8)  
  !   orientation of curve inside the cell given by cartesian index



  ! Local variables

  integer :: ibit, ip, ind_part
  integer :: sweep_size, sweep_offset, nsweep, isweep 
  integer(kind=4), dimension(1:nvector) :: nstate, sdigit, ind
  integer(kind=8), dimension(1:nvector) :: ix, iy, iz
  integer(kind=4), parameter :: one=1
  integer(kind=4), parameter :: two=2
  integer(kind=4), parameter :: four=4
  integer(kind=4), parameter :: eight=8
  integer(kind=8), parameter :: longeight=8
  real(dp) :: ckey_factor

  ! compute particle position to cartesian key factor
  ckey_factor = 2**final_level / dble(boxlen)

  ! loop particles in nvector sweeps
  nsweep = ceiling( 1.*np / nvector)

  do isweep = 1, nsweep
     sweep_offset = (isweep-1)*nvector
     sweep_size = min(nvector, np-sweep_offset)
     sweep_offset = sweep_offset + offset

     ! if starting from scratch, reset hilbert keys
     if (initial_level == 0) then
        part_hkey     (1+sweep_offset : sweep_size+sweep_offset, 0) = 0
        part_hkey     (1+sweep_offset : sweep_size+sweep_offset, 1) = 0 
        part_hkey     (1+sweep_offset : sweep_size+sweep_offset, 2) = 0         
        current_state (1+sweep_offset : sweep_size+sweep_offset   ) = 0
     end if

     ! compute cartesian keys
     do ip = 1, sweep_size 
        ix(ip) = int(xp_andreas(ip+sweep_offset,1)*ckey_factor, kind=8)
        iy(ip) = int(xp_andreas(ip+sweep_offset,2)*ckey_factor, kind=8)
        iz(ip) = int(xp_andreas(ip+sweep_offset,3)*ckey_factor, kind=8)
     end do
     
     ! loop over levels --> counter ibit determines which 
     ! bit of the cartesian key is processed into hilber key digits
     do ibit=final_level-initial_level-1,0,-1

        if (final_level>42)then
           do ind_part = 1+sweep_offset, sweep_size+sweep_offset
              part_hkey(ind_part,2) = ISHFT(part_hkey(ind_part,2), 3)
           end do
           do ind_part = 1+sweep_offset, sweep_size+sweep_offset
              part_hkey(ind_part,2) = part_hkey(ind_part,2) + ISHFT(part_hkey(ind_part,1), -60)
           end do
        end if

        if (final_level>21)then
           do ind_part = 1+sweep_offset, sweep_size+sweep_offset
              part_hkey(ind_part,1) = ISHFT(part_hkey(ind_part,1),  4)
              part_hkey(ind_part,1) = ISHFT(part_hkey(ind_part,1), -1)
           end do
           do ind_part = 1+sweep_offset, sweep_size+sweep_offset
              part_hkey(ind_part,1) = part_hkey(ind_part, 1) + ISHFT( part_hkey(ind_part, 0), -60)
           end do
        end if

        do ind_part = 1+sweep_offset, sweep_size+sweep_offset
           part_hkey(ind_part, 0) = ISHFT( part_hkey(ind_part, 0),  4)
           part_hkey(ind_part, 0) = ISHFT( part_hkey(ind_part, 0), -1)
        end do

        sdigit = 0
        do ip = 1, sweep_size
           if(btest(ix(ip),ibit)) sdigit(ip) = sdigit(ip) + four
           if(btest(iy(ip),ibit)) sdigit(ip) = sdigit(ip) + two
           if(btest(iz(ip),ibit)) sdigit(ip) = sdigit(ip) + one
        end do

        do ip=1,sweep_size
           ind(ip) = current_state(ip+sweep_offset)*eight + sdigit(ip)
        end do

        do ip=1,sweep_size
           nstate(ip) = next_state_diagram(ind(ip))
        end do

        do ip=1,sweep_size
           part_hkey(sweep_offset+ip,0) = part_hkey(sweep_offset+ip,0) + three_digit_diagram(ind(ip))
        end do

        do ip=1,sweep_size
           current_state(sweep_offset+ip) = nstate(ip)
        end do

     end do ! end loop over levels        
  end do ! end sweeps

end subroutine hilbert3d_for_particle
