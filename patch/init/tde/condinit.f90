!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(r,g,x,q,dx,nn)
  use amr_parameters, only: ndim, nvector
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t, global_t
  use input_hydro_condinit_module, only: region_condinit
  implicit none
  type(run_t)::r
  type(global_t)::g
  integer ::nn                            ! Number of cells
  real(kind=8)::dx                            ! Cell size
#ifdef MHD
  real(kind=8),dimension(1:nvector,1:nvar+3-ndim)::q ! Primitive variables
#else
  real(kind=8),dimension(1:nvector,1:nvar)::q ! Primitive variables
#endif
  real(kind=8),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user (aka code) units:
  ! x(i,1:ndim) are in [0,boxlen]**ndim.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w and Q(i,5): P.
  ! If nvar >= 6, remaining variables are treated as passive
  ! scalars or non-thermal energies in the hydro solver.
  ! For 1D MHD, Q(i,nvar+1) is By and Q(i,nvar+2) is Bz.
  ! For 2D MHD, Q(i,nvar+1) is Bz.
  ! Q(:,:) are in user (aka code) units.
  !================================================================
  integer::i
  logical,save::read_flag=.false.
  integer,parameter::nrows=10000,ncols=2          ! CSV file parameters
  real(kind=8),dimension(1:nrows, 1:ncols),save::xx   ! Lane-Emden solutions (r, d, p)
  integer,save::nmax
  real(kind=8),save::rmax,dmin,mass
  real(kind=8)::r2,rx,ry,rz,rr,d,p,pi
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m

  ! Read Lane-Emden solutions into an array
  if (.not. read_flag) then
     xx=0d0
     open (unit=10,file="lane_emden.csv",action="read",status="old")
     nmax=0
     do
        nmax=nmax+1
        read (10,*,iostat=io) xx(nmax,:)
        if(io.ne.0)exit
     end do
     read_flag = .true.
     ! Maximum radius
     rmax=xx(nmax-1,1)
     ! Minimum density
     dmin=xx(nmax-1,2)
     ! Minimum pressure
     pmin=dmin**(4d0/3d0)
     ! Polytrope mass
     mass=0d0
     do i=1,nmax-1
        mass=mass*xx(i,2)*4d0*pi*xx(i,1)**2*(xx(i+1,1)-xx(i,1))
     end do
     write(*,*)'Lane Emden file read'
     write(*,*)'nmax=',nmax,' rmax=',rmax,' mass=',mass,' dmin=',dmin,' pmin=',pmin
  end if

  ! Scale factors
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3
  pi=ACOS(-1d0)

  ! vacuum as default
  do i=1,nn
     q(i,1)=dmin
     q(i,2)=0.0
     q(i,3)=0.0
     q(i,4)=0.0
     q(i,5)=pmin
  end do

  if(rmax>dx)then
     do i=1,nn
        rx=x(i,1)-r%boxlen/2.
        ry=x(i,2)-r%boxlen/2.
        rz=x(i,3)-r%boxlen/2.
        r2=rx**2+ry**2+rz**2
        rr=sqrt(r2)
        if(rr<rmax)then
           j=(rr/rmax)*(nmax-1)
           ! density
           d=xx(j,2)+(rr-xx(j,1))*((xx(j+1,2)-xx(j,2))/(xx(j+1,1)-xx(j,1)))
           ! pressure
           p=d**(4d0/3d0)
           ! store primitive variables
           q(i,1)=d
           q(i,5)=p
        endif
     end do
  else
     do i=1,nn
        rx=x(i,1)-r%boxlen/2.
        ry=x(i,2)-r%boxlen/2.
        rz=x(i,3)-r%boxlen/2.
        r2=rx**2+ry**2+rz**2
        rr=sqrt(r2)
        if(rr<dx)then
           ! density
           d=mass/dx**3
           ! pressure
           p=pmin
           ! store primitive variables
           q(i,1)=d
           q(i,5)=p
        end if
     end do
  endif

  ! Compute entropy if needed
  if(r%entropy)then
     q(1:nn,r%ientropy)=q(1:nn,5)/q(1:nn,1)**r%gamma
  endif

  ! Compute metallicity if needed
  if(r%metal)then
     q(1:nn,r%imetal)=r%z_ave*0.02
  endif

end subroutine condinit
