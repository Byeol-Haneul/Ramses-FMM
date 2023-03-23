module star_formation_module

  type :: out_star_formation_t
     real(kind=8)::mass
  end type out_star_formation_t

contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_star_formation(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(out_star_formation_t)::output,next_output

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_STAR_FORMATION,pst%iUpper+1,input_size,output_size,ilevel)
     call r_star_formation(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
  else
     call star_formation(pst%s%r,pst%s%g,pst%s%m,pst%s%s,ilevel,output%mass)
  endif

end subroutine r_star_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine star_formation(r,g,m,s,ilevel,mstar_loc)
  use rng
  use constants
  use hydro_parameters, only:nvar
  use amr_parameters, only:dp,ndim,twotondim
  use amr_commons, only:run_t,global_t,mesh_t
  use pm_commons, only:part_t
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::s
  integer::ilevel
  real(kind=8)::mstar_loc
  !-------------------------------------------------------------------
  ! Spawn star particles according to various star formation models.
  ! We use a random Poisson process.
  !-------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:g%ncpu)::nsite_cpu_tot,nstar_cpu_tot
#endif
  integer(kind=8),dimension(0:g%ncpu)::nsite_cum,nstar_cum
  integer,dimension(1:g%ncpu)::nsite_cpu,nstar_cpu
  integer::i,ind,igrid,idim,icpu,ngrid,nleaf,nsite,nstar,nstar_loc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dx,vol,factG,n_star,nCOM,d,d0,mstar,dstar,tstar,mcell,mgas,PoissMean
#ifdef SOLVERmhd
  integer::neul=5
#else
  integer::neul=ndim+2
#endif
#if NENER>0
  integer::irad
#endif
  type(RngStream)::RngStream_CreateStream, gg
  real(kind=8)::RngStream_RandUni
  logical::ok

#ifdef HYDRO
  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Set some constants
  dx=r%boxlen/2**ilevel
  vol=dx**ndim
  factG=1d0
  if(r%cosmo)factG=3d0/4d0/twopi*g%omega_m*g%aexp

  ! Density threshold for star formation
  n_star = r%n_star

  ! This is used only for cosmological runs
  if(r%cosmo)then
     nCOM = 200.0*g%omega_b*rhoc*(g%h0/100)**2/g%aexp**3/mH
     n_star = MAX(nCOM,n_star)
  endif
  d0 = n_star/scale_nH
  mstar = r%m_star*r%mass_sph
  dstar = mstar/vol

  !---------------------------------------------------------
  ! Count potential star formation sites.
  ! This is important for the pseudo-random number generator
  !---------------------------------------------------------
  nsite=0
  ! Loop over octs with vector sweeps
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
        ! Select leaf cells
        ok = .not. m%grid(igrid)%refined(ind)
        ! Select dense enough cells
        d = m%grid(igrid)%uold(ind,1)
        ok = ok .and. d > d0
        ! Count number of random numbers
        if(ok)then
           nsite=nsite+1
        endif
     end do
  end do
  
  !---------------------------------------------------------
  ! Compute number of star formation sites across all CPUs.
  !---------------------------------------------------------
  nsite_cpu=0
  nsite_cpu(g%myid)=nsite
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nsite_cpu,nsite_cpu_tot,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nsite_cpu=nsite_cpu_tot
#endif
  nsite_cum=0
  do icpu=1,g%ncpu
     nsite_cum(icpu)=nsite_cum(icpu-1)+int(nsite_cpu(icpu),kind=8)
  end do

  !------------------------------------------------
  ! Initialize RNG stream based on the global seed
  !------------------------------------------------
  call RngStream_SetPackageSeed(r%seed)
  ! Get current rng stream state
  gg = RngStream_CreateStream('ramses_stream')
  ! Advance to the state corresponding to the current cpu
  if(g%myid>1)then
     if(nsite_cum(g%myid-1)>0)then
        call RngStream_AdvanceState(gg,0_8,nsite_cum(g%myid-1))
     endif
  endif

  nstar_loc=0
  mstar_loc=0.0d0
  ! Loop over octs
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
        ! Select leaf cells
        ok = .not. m%grid(igrid)%refined(ind)
        ! Select dense enough cells
        d = m%grid(igrid)%uold(ind,1)
        ok = ok .and. d > d0
        ! Draw Poisson process based on local SF rate
        if(ok)then
           mcell=d*vol
           tstar=0.5427d0/sqrt(factG*d)
           mgas=g%dtnew(ilevel)*(r%eps_star/tstar)*mcell
           ! Poisson mean
           PoissMean=mgas/mstar
           call poissdev(RngStream_RandUni(gg),PoissMean,nstar)
           ! Check that we don't empty the gas cell
           mgas=nstar*mstar
           if(mgas>0.9*mcell)then
              nstar=int(0.9*mcell/mstar)
           endif
           ! Create a new star particle
           if(nstar>0)then
              nstar_loc=nstar_loc+1
              mstar_loc=mstar_loc+nstar*mstar
              dstar=nstar*mstar/vol
              s%npart=s%npart+1
              if(s%npart>r%nstarmax)then
                 write(*,*)'Not enough memory for star particle'
                 write(*,*)'Increase nstarmax in the namelist'
                 stop
              endif
              ! Compute star particle coordinate from cell centers
              s%xp(s%npart,1)=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx-m%skip(1)
              s%xp(s%npart,2)=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx-m%skip(2)
              s%xp(s%npart,3)=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx-m%skip(3)
              ! Compute star particle velocity from gas velocity
              s%vp(s%npart,1)=m%grid(igrid)%uold(ind,2)/d
              s%vp(s%npart,2)=m%grid(igrid)%uold(ind,3)/d
              s%vp(s%npart,3)=m%grid(igrid)%uold(ind,4)/d
#ifdef GRAV
              ! Remove half a kick (will be added later)
              s%vp(s%npart,1)=s%vp(s%npart,1)-m%grid(igrid)%f(ind,1)*0.5d0*g%dtnew(ilevel)
              s%vp(s%npart,2)=s%vp(s%npart,2)-m%grid(igrid)%f(ind,2)*0.5d0*g%dtnew(ilevel)
              s%vp(s%npart,3)=s%vp(s%npart,3)-m%grid(igrid)%f(ind,3)*0.5d0*g%dtnew(ilevel)
#endif
              ! Compute star particle mass
              s%mp(s%npart)=nstar*mstar
              ! Compute star particle birth time
              if(r%cosmo)then
                 s%tp(s%npart)=g%texp
              else
                 s%tp(s%npart)=g%t
              endif
              ! Compute star particle metallicity
              if(r%metal)then
                 s%zp(s%npart)=m%grid(igrid)%uold(ind,r%imetal)/d
              else
                 s%zp(s%npart)=r%z_ave*0.02
              endif
              ! Compute level
              s%levelp(s%npart)=ilevel
              ! Remove corresponding mass from all conservative variables
              ! Careful: this will not work for MHD or non-thermal energies
              ! Best is to remove the corresponding energies before and add them back after
              m%grid(igrid)%uold(ind,1:nvar)=m%grid(igrid)%uold(ind,1:nvar)*(d-dstar)/d
           endif
        endif
     end do
  end do
  s%tailp(r%nlevelmax)=s%tailp(r%nlevelmax)+nstar_loc

  !----------------------------------------------
  ! Get the new global seed to the final state
  !----------------------------------------------
  call RngStream_SetPackageSeed(r%seed)
  gg = RngStream_CreateStream('ramses_stream')
  ! Advance to the state corresponding to the final global state
  if(nsite_cum(g%ncpu)>0)then
     call RngStream_AdvanceState(gg,0_8,nsite_cum(g%ncpu))
  endif
  ! Get the new seed and store it
  call RngStream_GetState(gg,r%seed)

  !---------------------------------------------------------
  ! Compute number of new stars across all CPUs.
  !---------------------------------------------------------
  nstar_cpu=0
  nstar_cpu(g%myid)=nstar_loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nstar_cpu,nstar_cpu_tot,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nstar_cpu=nstar_cpu_tot
#endif
  nstar_cum=0
  do icpu=1,g%ncpu
     nstar_cum(icpu)=nstar_cum(icpu-1)+int(nstar_cpu(icpu),kind=8)
  end do

  !--------------------------------------
  ! Compute new star particle index
  !--------------------------------------
  do i=s%npart-nstar_loc+1,s%npart
     s%idp(i)=s%npart_tot+nstar_cum(g%myid-1)+i
  end do
  s%npart_tot=s%npart_tot+nstar_cum(g%ncpu)

#endif

end subroutine star_formation
end module star_formation_module
