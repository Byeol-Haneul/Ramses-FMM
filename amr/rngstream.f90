!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
!         Implementacion del algoritmo de L'Ecuyer en FORTRAN para el   !
!                        generador, MRG32k3a.                           !
!         Se trata de una copia casi directa del algoritmo en C del     !
!      propio Ecuyer, que se puede encontrar en su pagina web:          !
!            http://www.iro.umontreal.ca/~lecuyer                       !
!      http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c2010/      !
!        La documentacion esta adjunta en el documento x.pdf            !
!                                                                       !
! Pueden haber errores en este codigo fuente, se advierte que se tenga  !
! precaucion y se hagan los test correspondientes. No aseguro que este  !
!                        libre de fallos.                               !
!                                                                       !
!                                             Pablo Serna Martinez      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module rng
  implicit none
  
  private
  public :: norm,m1,m2,a12,a13n,a21,a23n,two17,two53,fact &
       ,InvA1, InvA2, A1p0, A2p0, A1p76, A2p76, A1p127, A2p127 &
       ,RngStream, state, rng_off

  !  Dimension del estado   
  integer, parameter :: ns = 6
  !  parametros   
  real (KIND=8), parameter :: norm = 2.328306549295727688d-10
  real (KIND=8), parameter :: m1 = 4294967087.0d0
  real (KIND=8), parameter :: m2 = 4294944443.0d0
  real (KIND=8), parameter :: a12 = 1403580.0d0
  real (KIND=8), parameter :: a13n = 810728.0d0
  real (KIND=8), parameter :: a21 = 527612.0d0
  real (KIND=8), parameter :: a23n = 1370589.0d0         
  real (KIND=8), parameter :: two17 = 131072.0d0
  real (KIND=8), parameter :: two53 = 9007199254740992.0d0
  real (KIND=8), parameter :: fact = 5.9604644775390625d-8    ! 2^{-24}
  integer(KIND=8), parameter, dimension(ns) :: default_seed &
       =(/ 12345,12345,12345,12345,12345,12345/)
  logical                  :: rng_off = .TRUE.
  
  !Transition matrices of two MRG components, elevadas a -1, 1, 2^76 y 2^127
  
  !Nota: column-major order
  real (KIND=8), dimension(3, 3) :: InvA1, InvA2, A1p0, A2p0, &
       A1p76, A2p76, A1p127, A2p127  
  
  integer (KIND=8), dimension(ns) :: state = default_seed
  
  type :: RngStream
     real (KIND=8), dimension(6)    :: Cg, Bg, Ig
     integer (KIND=4)               :: length
     character (LEN=256)            :: cname
     logical                        :: Anti=.FALSE.
     logical                        :: IncPrec=.FALSE.
  end type RngStream
end module rng

subroutine RngInit()
  use rng
  
  InvA1 = reshape( &
       (/184888585.0d0,1.0d0,0.0d0, &
       0.0d0,0.0d0,1.0d0, &
       1945170933.0d0,0.0d0,0.0d0/), shape(InvA1))
  InvA2 = reshape( &
       (/0.0d0,1.0d0,0.0d0, &
       360363334.0d0,0.0d0,1.0d0, &
       4225571728.0d0,0.0d0,0.0d0/), shape(InvA1))
  A1p0 = reshape( &
       (/0.0d0,0.0d0,-810728.0d0, &
       1.0d0,0.0d0,1403580.0d0, &
       0.0d0,1.0d0,0.0d0/), shape(InvA1))
  A2p0 = reshape( &
       (/0.0d0,0.0d0,-1370589.0d0, &
       1.0d0,0.0d0,0.0d0, &
       0.0d0,1.0d0,527612.0d0/), shape(InvA1))   
  A1p76 = reshape( &
       (/82758667.0d0, 3672831523.0d0,3672091415.0d0, &
       1871391091.0d0,69195019.0d0,3528743235.0d0, &
       4127413238.0d0,1871391091.0d0,69195019.0d0/), shape(InvA1))
  A2p76 = reshape( &
       (/ 1511326704.0d0, 4292754251.0d0, 3859662829.0d0, &
       3759209742.0d0, 1511326704.0d0, 4292754251.0d0, &
       1610795712.0d0, 3889917532.0d0, 3708466080.0d0/), shape(InvA1))
  A1p127 = reshape( &
       (/ 2427906178.0d0,  226153695.0d0, 1988835001.0d0, &
       3580155704.0d0, 1230515664.0d0,  986791581.0d0, &
       949770784.0d0, 3580155704.0d0, 1230515664.0d0/), shape(InvA1))
  A2p127 = reshape( &
       (/ 1464411153.0d0,   32183930.0d0, 2824425944.0d0, &
       277697599.0d0, 1464411153.0d0,   32183930.0d0, &
       1610723613.0d0, 1022607788.0d0, 2093834863.0d0/), shape(InvA1))   
  
  return
end subroutine RngInit

function MultModM(a, s, c, m) result(v)
  ! Calcula mod(a*s + c,m) . Con m <2^35. tambien va para s,c<0
  use rng
  implicit none
  real(KIND=8)    :: a, s, c, m
  real(KIND=8)    :: ap, sp, cp   
  real(KIND=8)    :: v
  integer(KIND=8) :: a1
  
  ap=a
  v = ap*s+c
  if(v.ge.two53.or.v.le.-two53) then
     a1 = int(ap/two17,KIND=8)
     ap = ap-dble(a1)*two17
     v = dble(a1) * s
     a1 = int( v/m , KIND=8)
     v = v - dble(a1)*m
     v = v*two17 + ap*s + c
  endif
  a1 = int(v/m,KIND=8);
  v = v - dble(a1)*m
  if(v.lt.0.0d0) v = v + m
  
  return
end function MultModM

subroutine MatVecModM (a, s, v, m) 
  ! Calcula v= mod(A*s ,m). Asume -m < s(i) < m. Va incluso si v=s
  implicit none
  real(KIND=8)  :: a(3,3), s(3), m, x(3)
  real(KIND=8)  :: v(3)
  real(KIND=8)  :: MultModM
  integer       :: i
  
  do i=1,3
     x(i) = MultModM(a(i,1),s(1),0.0d0,m)
     x(i) = MultModM(a(i,2),s(2),x(i),m)
     x(i) = MultModM(a(i,3),s(3),x(i),m)  
  enddo
  
  v(1:3)=x(1:3)
  
  return
end subroutine MatVecModM

subroutine MatMatModM (a, b, c, m) 
  ! Returns C = A*B % m. Work even if A = C or B = C or A = B = C.
  implicit none
  real(KIND=8)     :: a(3,3), b(3,3), c(3,3),m
  integer          :: i, j
  real(KIND=8)     :: v(3),v1(3),w(3,3)
  
  w=0.0d0
  
  do i=1,3
     do j=1,3
        v(j)=b(j,i)
     enddo
     v1=v
     
     call MatVecModM(a, v1, v, m)
     do j=1,3
        w(j,i)=v(j)
     enddo
  enddo
  c(:,:)=w(:,:)
  
  return
end subroutine MatMatModM

subroutine MatTwoPowModM (a, b, m, e)
  !   Compute matrix B = (A^(2^e) % m);  works even if A = B 
  implicit none
  real(KIND=8)     :: a(3,3), b(3,3), b1(3,3), m
  integer(KIND=8)  :: e
  integer(KIND=8)  :: i
  
  b(:,:)=a(:,:)
  b1=a
  do i=1,e
     call MatMatModM(b1, b1, b, m)
     b1=b
  enddo
  
  
  return
end subroutine MatTwoPowModM

subroutine MatPowModM (a, b, m, n)
  !     /* Compute matrix B = A^n % m ;  works even if A = B */
  implicit none
  real(KIND=8)    :: a(3,3),b(3,3),b1(3,3),m,w(3,3),w1(3,3)
  integer(KIND=8) :: n,n2
  integer         :: j
  
  n2 = n
  w(:,:) = a(:,:)
  b(:,:) = 0.0d0
  
  do j=1,3
     b(j,j) = 1.0d0
  enddo
  b1=b
  w1=w
  !       /* Compute B = A^n % m using the binary decomposition of n */
  do while(n2.gt.0)
     if(mod(n2,2).eq.1) then
        call MatMatModM(w, b1, b, m)
        b1=b
     endif
     call MatMatModM(w1, w1, w, m)
     w1=w
     n2=n2/2
  enddo
  
  return
end subroutine MatPowModM

double precision function  U01 (g)
  use rng 
  implicit none
  integer(KIND=8) :: k
  real(KIND=8)    :: p1,p2,u
  type(RngStream) :: g
  ! Componente 1  
  p1 = a12 * g%Cg(2) - a13n * g%Cg(1)
  k = p1/m1
  
  p1 = p1 - k*m1
  if(p1.lt.0.0d0) p1 = p1+m1
  
  g%Cg(1) = g%Cg(2)
  g%Cg(2) = g%Cg(3)
  g%Cg(3) = p1
  
  ! Componente 2
  p2 = a21 * g%Cg(6) - a23n * g%Cg(4)
  k = p2/m2
  p2 = p2-k*m2
  if(p2.lt.0.0d0) p2 = p2+m2
  
  g%Cg(4) = g%Cg(5)
  g%Cg(5) = g%Cg(6)
  g%Cg(6) = p2
  
  ! Combinacion
  
  if(p1.gt.p2) then
     u = (p1-p2)*norm
  else
     u = (p1-p2+m1)*norm
  endif
  
  if(g%Anti) then
     U01 = 1 - u
  else
     U01 = u
  endif
  return
end function U01

double precision function  U01d (g)
  use rng
  implicit none
  real (KIND=8)   :: u
  real (KIND=8)   :: U01
  type(RngStream) :: g
  external  U01
  
  u = U01(g)
  if(g%Anti) then
     ! Antithetic case
     u = u + (U01(g)-1.0d0)*fact  
     if(u.lt.0.0d0) u = u+1.0d0
  else
     u = u + U01(g)*fact
     if(u.gt.1.0d0) u = u-1.0d0
  endif
  
  U01d = u
  return
end function U01d

subroutine CheckSeed(seed,ierror)
  use rng
  integer (KIND=4)      :: i,ierror
  integer (KIND=8)      :: seed(6)
  
  
  ierror=0
  do i=1,3
     if(seed(i).gt.int(m1,KIND=8)) then
        write(*,'(A,/,A,I1,A,/,A,/)')  &
             '***************************************************' &
             ,'ERROR: Seed[',i,']>=m1, Seed is not set.' &
             ,'***************************************************' 
        ierror=1
     endif
  enddo
  do i=4,6
     if(seed(i).gt.int(m2,KIND=8)) then
        write(*,'(A,/,A,I1,A,/,A,/)')  &
             '***************************************************' &
             ,'ERROR: Seed[',i,']>=m2, Seed is not set.' &
             ,'***************************************************' 
        ierror=1     
     endif
  enddo
  
  if(seed(1).eq.0.and.seed(2).eq.0.and.seed(3).eq.0) then
     write(*,'(A,/,A,/,A,/)')  &
          '***************************************************' &
          ,'ERROR: First 3 seed are 0' &
          ,'***************************************************' 
  endif
  
  if(seed(4).eq.0.and.seed(5).eq.0.and.seed(6).eq.0) then
     write(*,'(A,/,A,/,A,/)')  &
          '***************************************************' &
          ,'ERROR: First 3 seed are 0' &
          ,'***************************************************' 
     ierror=1     
  endif
  return
end subroutine CheckSeed

!!!!!!!!!!!!!!! LA COSA PUBLICA

function RngStream_CreateStream (cname)  result(g)
  use rng
  integer (KIND=4)      :: i,length,ie
  real (KIND=8)         :: seed(3)
  character (LEN=*)     :: cname
  character (LEN=256)   :: cnameb
  type(RngStream) :: g
  
  if(rng_off) then
     call RngInit()
     rng_off=.FALSE.
  endif
  cnameb=cname
  g%length=LEN_TRIM(cnameb)
  if(g%length.gt.0) then
     g%cname = cname(1:g%length)
  else
     g%cname = '0'
  endif
  !  g%Anti = .FALSE.
  !  g%IncPrec = .FALSE.
  
  g%Bg(:)=state(:)
  g%Cg(:)=state(:)
  g%Ig(:)=state(:)    
  seed(1:3)=dble(state(1:3))
  call MatVecModM(A1p127,seed, seed, m1)
  state(1:3)=int(seed(1:3)+1d-14,KIND=8)
  seed(1:3)=dble(state(4:6))
  call MatVecModM(A2p127,seed, seed, m2)
  state(4:6)=int(seed(1:3)+1d-14,KIND=8)
  return
end function RngStream_CreateStream

!!! BORRAR STREAM --> deallocate(g)

subroutine RngStream_ResetStartStream (g)
  use rng
  type(RngStream) :: g
  
  g%Cg(:)=g%Ig(:)
  g%Bg(:)=g%Ig(:)   
  
  return
end subroutine RngStream_ResetStartStream

subroutine RngStream_ResetNextSubstream (g)
  use rng
  type(RngStream) :: g
  
  call MatVecModM(A1p76, g%Bg(1:3), g%Bg(1:3), m1)
  call MatVecModM(A2p76, g%Bg(4:6), g%Bg(4:6), m2)   
  g%Cg(:)=g%Bg(:)
  
  return
end subroutine RngStream_ResetNextSubstream

subroutine RngStream_ResetStartSubstream (g)
  use rng
  type(RngStream) :: g
  
  g%Cg(:)=g%Bg(:)
  
  return
end subroutine RngStream_ResetStartSubstream

subroutine RngStream_SetPackageSeed (seed)
  use rng
  integer(KIND=8)      :: seed(6)
  integer(KIND=4)      :: ierror
  call CheckSeed(seed,ierror)
  if(ierror.eq.1) stop 'ERROR'
  
  state(:)=seed(:)
  return
end subroutine RngStream_SetPackageSeed

subroutine RngStream_SetSeed (g,seed,ierror)
  use rng
  integer(KIND=8)      :: seed(6)
  integer              :: ierror
  type(RngStream)      :: g
  call CheckSeed(seed,ierror)
  if(ierror.eq.1) stop 'ERROR'
  
  g%Cg(:)=seed(:)   
  g%Bg(:)=seed(:)
  g%Ig(:)=seed(:)      
  
  return
end subroutine RngStream_SetSeed

subroutine RngStream_AdvanceState(g,e,c)
  use rng
  integer(KIND=8)              :: e, c
  real (KIND=8),dimension(3,3) :: b1, c1, b2, c2, ctemp
  real (KIND=8)                :: temp(3)
  type(RngStream)              :: g
    
  if(e.gt.0) then
     call MatTwoPowModM (A1p0, b1, m1, e)
     call MatTwoPowModM (A2p0, b2, m2, e)    
  else if(e.lt.0) then
     call MatTwoPowModM (InvA1, b1, m1, -e)
     call MatTwoPowModM (InvA2, b2, m2, -e)    
  endif

  if(c.gt.0) then
     call MatPowModM (A1p0, c1, m1, c)
     call MatPowModM (A2p0, c2, m2, c)    
  else 
     call MatPowModM (InvA1, c1, m1, -c)
     call MatPowModM (InvA2, c2, m2, -c)    
  endif

  
  if(e.ne.0) then
     ctemp=c1
     call MatMatModM (b1, ctemp, c1, m1)
     ctemp=c2    
     call MatMatModM (b2, ctemp, c2, m2)    
  endif
  
  temp=g%Cg(1:3)
  call MatVecModM (c1, g%Cg(1:3), temp,m1)
  g%Cg(1:3)=temp
  temp=g%Cg(4:6)
  call MatVecModM (c2, g%Cg(4:6), temp,m2)   
  g%Cg(4:6)=temp
  !stop
  return
end subroutine RngStream_AdvanceState


subroutine RngStream_GetState (g, seed)
  use rng
  type(RngStream) :: g
  integer(KIND=8) :: seed(6)
  
  seed(:) = g%Cg(:)
  return
end subroutine RngStream_GetState


subroutine RngStream_WriteState (g)
  use rng
  type(RngStream) :: g
  
  if(g%length.eq.0) return
  
  write(*,'(A,A,A)') 'The current state of the Rngstream ', &
       g%cname, ': '
  
  write(*,'(A,6(I20,X),A)') 'Cg = { ', g%Cg(1:6),'}'
  
end subroutine RngStream_WriteState


subroutine RngStream_WriteStateFull (g)
  use rng
  type(RngStream) :: g
  
  if(g%length.eq.0) return
  
  write(*,'(A,A,A)') 'The Rngstream ', &
       g%cname, ': '
  
  write(*,'(A,L1)') ' Anti = ',g%Anti
  write(*,'(A,L1)') ' IncPrec = ',g%IncPrec
  
  write(*,'(A,6(I20,X),A)') 'Ig = { ', g%Ig(1:6),'}'
  write(*,'(A,6(I20,X),A)') 'Bg = { ', g%Bg(1:6),'}'      
  write(*,'(A,6(I20,X),A)') 'Cg = { ', g%Cg(1:6),'}'
  
end subroutine RngStream_WriteStateFull

subroutine RngStream_IncreasedPrecis (g, incp)
  use rng
  type(RngStream) :: g
  logical         :: incp
  
  g%IncPrec = incp
  return
end subroutine RngStream_IncreasedPrecis

subroutine RngStream_SetAntithetic (g, a)
  use rng
  type(RngStream) :: g
  logical         :: a
  
  g%Anti = a
  return
end subroutine RngStream_SetAntithetic

function RngStream_RandUni (g) result(r)
  use rng
  type(RngStream) :: g
  real (KIND=8)   :: U01d,U01,r
  
  if(g%IncPrec) then
     r=U01d(g)
  else
     r=U01(g)
  endif
  
  return 
end function RngStream_RandUni

function RngStream_RandInt (g,i,j) result(l)
  use rng
  type(RngStream) :: g
  integer(KIND=4) :: i,j,l
  real (KIND=8)   :: RngStream_RandUni
  
  l= i+ int((dble(j)-dble(i)+1.0d0)*RngStream_RandUni(g))
  return 
end function RngStream_RandInt

SUBROUTINE gaussdev(p, normal_dev)
  
  ! ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
  
  ! Produces the normal deviate Z corresponding to a given lower
  ! tail area of P; Z is accurate to about 1 part in 10**16.
  
  ! The hash sums below are the sums of the mantissas of the
  ! coefficients.   They are included for use in checking
  ! transcription.

  ! This ELF90-compatible version by Alan Miller - 20 August 1996
  ! N.B. The original algorithm is as a function; this is a subroutine
  
  REAL (kind=8), INTENT(IN)   :: p
  REAL (kind=8), INTENT(OUT)  :: normal_dev
  
  ! Local variables
  
  REAL (kind=8) :: zero = 0.d0, one = 1.d0, half = 0.5d0, split1 = 0.425d0,  &
       split2 = 5.d0, const1 = 0.180625d0, const2 = 1.6d0, q, r
  
  ! Coefficients for P close to 0.5
  
  REAL (kind=8) :: a0 = 3.3871328727963666080D0, &
       a1 = 1.3314166789178437745D+2, &
       a2 = 1.9715909503065514427D+3, &
       a3 = 1.3731693765509461125D+4, &
       a4 = 4.5921953931549871457D+4, &
       a5 = 6.7265770927008700853D+4, &
       a6 = 3.3430575583588128105D+4, &
       a7 = 2.5090809287301226727D+3, &
       b1 = 4.2313330701600911252D+1, &
       b2 = 6.8718700749205790830D+2, &
       b3 = 5.3941960214247511077D+3, &
       b4 = 2.1213794301586595867D+4, &
       b5 = 3.9307895800092710610D+4, &
       b6 = 2.8729085735721942674D+4, &
       b7 = 5.2264952788528545610D+3
  ! HASH SUM AB    55.8831928806149014439
  
  ! Coefficients for P not close to 0, 0.5 or 1.
  
  REAL (kind=8) :: c0 = 1.42343711074968357734D0, &
       c1 = 4.63033784615654529590D0, &
       c2 = 5.76949722146069140550D0, &
       c3 = 3.64784832476320460504D0, &
       c4 = 1.27045825245236838258D0, &
       c5 = 2.41780725177450611770D-1, &
       c6 = 2.27238449892691845833D-2, &
       c7 = 7.74545014278341407640D-4, &
       d1 = 2.05319162663775882187D0, &
       d2 = 1.67638483018380384940D0, &
       d3 = 6.89767334985100004550D-1, &
       d4 = 1.48103976427480074590D-1, &
       d5 = 1.51986665636164571966D-2, &
       d6 = 5.47593808499534494600D-4, &
       d7 = 1.05075007164441684324D-9
  ! HASH SUM CD    49.33206503301610289036
  
  ! Coefficients for P near 0 or 1.
  
  REAL (kind=8) :: e0 = 6.65790464350110377720D0, &
       e1 = 5.46378491116411436990D0, &
       e2 = 1.78482653991729133580D0, &
       e3 = 2.96560571828504891230D-1, &
       e4 = 2.65321895265761230930D-2, &
       e5 = 1.24266094738807843860D-3, &
       e6 = 2.71155556874348757815D-5, &
       e7 = 2.01033439929228813265D-7, &
       f1 = 5.99832206555887937690D-1, &
       f2 = 1.36929880922735805310D-1, &
       f3 = 1.48753612908506148525D-2, &
       f4 = 7.86869131145613259100D-4, &
       f5 = 1.84631831751005468180D-5, &
       f6 = 1.42151175831644588870D-7, &
       f7 = 2.04426310338993978564D-15
  ! HASH SUM EF    47.52583317549289671629
  
  q = p - half
  IF (ABS(q) <= split1) THEN
     r = const1 - q * q
     normal_dev = q * (((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0) / &
          (((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + one)
     RETURN
  ELSE
     IF (q < zero) THEN
        r = p
     ELSE
        r = one - p
     END IF
     IF (r <= zero) THEN
        normal_dev = zero
        RETURN
     END IF
     r = SQRT(-LOG(r))
     IF (r <= split2) THEN
        r = r - const2
        normal_dev = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0) / &
             (((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + one)
     ELSE
        r = r - split2
        normal_dev = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0) / &
             (((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + one)
     END IF
     IF (q < zero) normal_dev = - normal_dev
     RETURN
  END IF
  RETURN
END SUBROUTINE gaussdev

subroutine poissdev(RandNum, AverNum, PoissNum)
  
  ! This is the original Knuth's algorithm for a Poisson law.  
  ! Choosing a large value for NPoissonLimit is more accurate but slower.
  ! For large N we fall back to the Gaussian deviate using the same
  ! input probability and the inverse CDF method.
  
  implicit none
  real(kind=8) :: RandNum
  real(kind=8) :: AverNum
  integer :: PoissNum

  integer,parameter :: NPoissonLimit = 10
  real(kind=8) :: Norm, Repar, Proba
  real(kind=8) :: GaussNum
  
  if(AverNum <= DBLE(NPoissonLimit)) then
     Norm=exp(-AverNum)
     Repar=1.0d0
     PoissNum=0
     Proba=1.0d0
     do while(Repar*Norm <= RandNum .and. PoissNum <= 10*NPoissonLimit )
        PoissNum=PoissNum+1
        Proba=Proba*AverNum/PoissNum
        Repar=Repar+Proba
     end do
  else
     call gaussdev(RandNum,GaussNum)
     GaussNum=GaussNum*sqrt(AverNum)-0.5d0+AverNum
     if(GaussNum<=0.0d0)GaussNum=0.0d0
     PoissNum=nint(GaussNum)
  endif
  
  return
end subroutine poissdev
