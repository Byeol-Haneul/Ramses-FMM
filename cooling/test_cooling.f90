program test_cool
  use constants
  use amr_parameters, only: dp, ndim, nvector, twotondim
  use hydro_parameters, only: nener, nion
  use rt_parameters, only: nrtgrp, smallnp
  use amr_commons, only: run_t, global_t, mesh_t
  use cooling_module, only: cooling_t, solve_cooling, T2_min_fix, set_table
  use coolrates_module, only: neq_cooling_t
  use neq_cooling_module, only: neq_solve_cooling, getmu
  use init_neq_chem_module, only: init_neq_chem

  type(run_t) :: r
  type(global_t) :: g
  type(cooling_t) :: c
  type(neq_cooling_t) :: tables
  integer, parameter :: nnH=6, nT=5, nXHi=5, nIt=1000
  real(dp),dimension(1:nvector):: Tmu, TK
  real(dp),dimension(1:nion, 1:nvector):: xion
  real(dp),dimension(1:nvector):: nH, Zsolar

  real(dp),dimension(nnH, nT, nXHi, nIt, 1+2*nion):: cells
  real(dp),dimension(nIt)::times

  real(dp) :: nH0, nH1, T0, T1, x0, x1, dlognH, dlogT, dx, mu
  real(dp) :: nHeCell, cooltime, dtcool, time0, time1, dlogTime
  integer :: i, j, k, l
  integer :: icount
  
  g%myid = 1

  r%isHe = .true.
  r%isH2 = .false.
  r%neq_isTconst = .true.
  
  ! Set species variable index
  iCount = 0
  ! HI fraction
  if(r%isH2) then
     iCount = iCount+1
     r%ixHI = iCount
  endif
  ! HII fraction
  iCount = iCount+1
  r%ixHII = iCount
  ! HeII and HeIII fractions
  if(r%isHe) then
     iCount = iCount+1
     r%ixHeII = iCount
     iCount = iCount+1
     r%ixHeIII = iCount
  endif

  ! Initiate non-eq tables
  call init_neq_chem(r, g, tables)

  ! Initiate time intervals-----------------------------------------------
  cooltime = Myr2sec*1d3
  time0 = 0.001 * Myr2sec ; time1=cooltime
  dlogTime = (log10(time1)-log10(time0))/(nIt-1)
  do i = 0, nIt-1
     times(i+1) = 10**(log10(time0) + dlogTime * i)
  end do

  ! Initiate nH intervals-------------------------------------------------
  nH0 = 1.d-4 ; nH1 = 1.d6
  dlognH = (log10(nH1)-log10(nH0))/(nnH-1)
  do i = 0, nnH-1
     nH(i+1) = 10**(log10(nH0) + dlognH * i)
  end do

  ! Initiate Tmu intervals------------------------------------------------
  T0 = 1.d1; T1 = 1.d7
  dlogT = (log10(T1)-log10(T0))/(nT-1)
  do i = 0, nT-1
     cells(:,i+1,:,1,1) = 10**(log10(T0) + dlogT * i)
  end do

  ! Initiate xHII, xHeIII intervals---------------------------------------
  X0 = 0;  X1 = 1
  dX = (X1-X0)/(nXHi-1)
  do i = 0, nXHi-1
     cells(:,:,i+1,1,r%ixHII+1) = X0 + dX * i
     cells(:,:,i+1,1,r%ixHII+1) = min(max(cells(:,:,i+1,1,r%ixHII+1),1e-6),0.99999)
     if(r%isH2)then
        cells(:,:,i+1,1,r%ixHI+1) = 1.-cells(:,:,i+1,1,r%ixHII+1)
     endif
     if(r%isHe) then
        cells(:,:,i+1,1,r%ixHeII+1)  = 0
        cells(:,:,i+1,1,r%ixHeII+1) = min(max(cells(:,:,i+1,1,r%ixHeII+1),1e-6),0.99999)
        cells(:,:,i+1,1,r%ixHeIII+1) = X0 + dX * i
        cells(:,:,i+1,1,r%ixHeIII+1) = min(max(cells(:,:,i+1,1,r%ixHeIII+1),1e-6),0.99999)
     endif
  end do

  ! Loop over times
  do i = 2, nIt
     dtcool = times(i)-times(i-1)
     print *,i,' time [Myr] =',times(i)/Myr2sec
     ! Loop over initial ionization states
     do j = 1, nXHi
        ! Loop over initial temperatures
        do k = 1, nT
           if(r%isH2)then
              xion(r%ixHI,1:nnH)=cells(1:nnH,k,j,i-1,r%ixHI+1) ! xHI
           endif
           xion(r%ixHII,1:nnH)=cells(1:nnH,k,j,i-1,r%ixHII+1) ! xHII
           if(r%isHe) then
              xion(r%ixHeII,1:nnH) =cells(1:nnH,k,j,i-1,r%ixHeII+1) ! xHeII
              xion(r%ixHeIII,1:nnH)=cells(1:nnH,k,j,i-1,r%ixHeIII+1) ! xHeIII
           endif
           do l = 1, nnh
              mu = getMu(r,xion(:,l), TK(l))
              TK(l) = cells(l,k,j,i-1,1)
              Tmu(l) = TK(l)/mu ! T/mu
              if(r%neq_isTconst) r%neq_Tconst = TK(l)
           end do

           ! Solve cooling over target time step for all densities
           call neq_solve_cooling(r, tables, Tmu, xion, nH, Zsolar, &
                & dtcool, nnH)

           ! Store results in data cube
           if(r%isH2)cells(1:nnH,k,j,i,r%ixHI+1) = xion(r%ixHI,1:nnH)
           cells(1:nnH,k,j,i,r%ixHII+1) = xion(r%ixHII,1:nnH)
           if(r%isHe) cells(1:nnH,k,j,i,r%ixHeII+1) = xion(r%ixHeII,1:nnH)
           if(r%isHe) cells(1:nnH,k,j,i,r%ixHeIII+1) = xion(r%ixHeIII,1:nnH)
           do l = 1,nnh
              mu = getMu(r,xion(:,l), TK(l))
              cells(l,k,j,i,1) = Tmu(l)*mu
           end do

        end do
        ! End loop over Tmu
     end do
     ! End loop over initial fraction
  end do
  ! End loop over time

  ! Convert time intervals from seconds to megayears
  times = times / Myr2sec

  open(10,file="onecell_cooling.bin",access="stream",action="write",form="unformatted")
  write(10)nnH, nT, nXHi, nIt, nion+1
  write(10)nH(1:nnH)
  write(10)times
  write(10)cells
  close(10)
  
end program test_cool
