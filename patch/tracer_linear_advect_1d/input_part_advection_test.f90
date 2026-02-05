!====================================================================
! Tracer particle seeding for linear advection test (1D)
!====================================================================
subroutine input_trac_advection_test(r,g,p,npart_tot)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t, global_t
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  integer(kind=8)::npart_tot
  !----------------------------------------------------
  ! Seed tracer particles with variable counts per cell
  !----------------------------------------------------
  integer::icpu,ipart,ipercell,i1,i2,i3,idim
  integer::n1,n2,n3,n_low,n_high,n_per_cell
  integer,dimension(1:g%ncpu)::npart_loc
  integer(kind=8),dimension(1:g%ncpu+1)::start_ind
  integer(kind=8)::ipart_global
  real(kind=8)::dx,xx1,xx2,xx3,zz,cell_vol,rho
  real(kind=8),dimension(:),allocatable::tdx_low,tdy_low,tdz_low
  real(kind=8),dimension(:),allocatable::tdx_high,tdy_high,tdz_high

  !-------------------------------------------
  ! Mesh size at levelmin
  !-------------------------------------------
  n1 = 2**r%levelmin
  n2 = 1
  n3 = 1
  dx = 0.5d0**r%levelmin  ! Normalized cell size
  cell_vol = (r%boxlen*dble(dx))**ndim

  n_low = max(1,r%ntrac_per_cell)
  n_high = 2*n_low

  !--------------------------------------
  ! Compute starting index for each cpu
  !--------------------------------------
  p%npart_tot = npart_tot
  do icpu=1,g%ncpu+1
     start_ind(icpu)=((icpu-1)*npart_tot)/g%ncpu
     if(icpu>1)then
        npart_loc(icpu-1)=start_ind(icpu)-start_ind(icpu-1)
     endif
  end do
  p%npart = npart_loc(g%myid)

  ! Check that local number of tracer particles does not exceed maximum
  if(p%npart > r%ntracmax)then
     write(*,*)'ERROR: CPU ',g%myid,' has too many tracer particles: ',p%npart,' > ',r%ntracmax
     stop
  endif

  ! Precompute subcell offsets
  allocate(tdx_low(1:n_low),tdy_low(1:n_low),tdz_low(1:n_low))
  allocate(tdx_high(1:n_high),tdy_high(1:n_high),tdz_high(1:n_high))
  if(r%part_subcell_positions)then
     call part_subcell_positions(n_low,tdx_low,tdy_low,tdz_low)
     call part_subcell_positions(n_high,tdx_high,tdy_high,tdz_high)
  else
     tdx_low=0.d0 ; tdy_low=0.d0 ; tdz_low=0.d0
     tdx_high=0.d0 ; tdy_high=0.d0 ; tdz_high=0.d0
  endif

  ipart = 1
  ipart_global = 0
  do i3=0,n3-1
     do i2=0,n2-1
        do i1=0,n1-1
           zz = (dble(i1)+0.5d0)/dble(n1)
           if(zz >= 1.0d0/3.0d0 .and. zz < 2.0d0/3.0d0)then
              n_per_cell = n_high
              rho = 2.0d0
           else
              n_per_cell = n_low
              rho = 1.0d0
           endif
           do ipercell=1,n_per_cell
              if(ipart_global>=start_ind(g%myid).and.ipart_global<start_ind(g%myid+1))then
                 ! Normalized coordinates [0,1], then convert to code units
                 xx1=(dble(i1)+0.5d0)/dble(n1)
                 xx2=(dble(i2)+0.5d0)/dble(n2)
                 xx3=(dble(i3)+0.5d0)/dble(n3)
                 if(n_per_cell==n_high)then
                    p%xp(ipart,1)=(xx1+tdx_high(ipercell)*dx)*r%boxlen
                    if(ndim>1)p%xp(ipart,2)=(xx2+tdy_high(ipercell)*dx)*r%boxlen
                    if(ndim>2)p%xp(ipart,3)=(xx3+tdz_high(ipercell)*dx)*r%boxlen
                 else
                    p%xp(ipart,1)=(xx1+tdx_low(ipercell)*dx)*r%boxlen
                    if(ndim>1)p%xp(ipart,2)=(xx2+tdy_low(ipercell)*dx)*r%boxlen
                    if(ndim>2)p%xp(ipart,3)=(xx3+tdz_low(ipercell)*dx)*r%boxlen
                 endif
                 p%vp(ipart,1:ndim)=0.0d0
                 p%vp(ipart,ndim)=1.0d0
                 p%mp(ipart)=rho*cell_vol/dble(n_per_cell)
                 p%idp(ipart)=ipart_global+1
                 ipart=ipart+1
              endif
              ipart_global=ipart_global+1
           end do
        end do
     end do
  end do

  ! Deallocate temporary arrays
  if(allocated(tdx_low))deallocate(tdx_low,tdy_low,tdz_low)
  if(allocated(tdx_high))deallocate(tdx_high,tdy_high,tdz_high)

  ! Periodic box
  do ipart=1,p%npart
     do idim=1,ndim
        if(r%periodic(idim))then
           if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%boxlen
           if(p%xp(ipart,idim)>=r%boxlen)p%xp(ipart,idim)=p%xp(ipart,idim)-r%boxlen
        end if
     end do
  end do

  ! Compute particle initial level
  do ipart=1,p%npart
     p%levelp(ipart)=r%levelmin
  end do

  ! Put all particles in levelmin
  p%headp=p%npart+1
  p%tailp=p%npart
  p%headp(r%levelmin)=1
  p%tailp(r%levelmin)=p%npart
  if(ANY(.not.r%periodic(1:ndim)))then
     p%headp(r%levelmin-1)=1
     p%tailp(r%levelmin-1)=0
  endif

end subroutine input_trac_advection_test
