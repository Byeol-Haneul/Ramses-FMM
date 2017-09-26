!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine courant_fine_2(r,g,m,ilevel)
  use amr_parameters, only: dp,nvector,ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
  real(kind=8),dimension(3)::comm_buffin,comm_buffout
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::ivar,idim,ind,igrid
  real(dp)::dt_lev,dx,vol
  real(kind=8)::mass_loc,ekin_loc,eint_loc,dt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,dt_all
  real(dp),dimension(1:nvar)::uu
  real(dp),dimension(1:ndim)::gg

#ifdef HYDRO

  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
  dt_all=g%dtnew(ilevel); dt_loc=dt_all

  ! Mesh spacing at that level
  dx=0.5D0**ilevel*r%boxlen
  vol=dx**ndim

  ! Loop over active grids by vector sweeps
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim                

        ! Gather leaf cells
        if(.NOT. m%grid(igrid)%refined(ind))then

           ! Gather hydro variables
           do ivar=1,nvar
              uu(ivar)=m%grid(igrid)%uold(ind,ivar)
           end do

           ! Gather gravitational acceleration
           gg=0.0d0
#ifdef GRAV
           if(r%poisson)then
              do idim=1,ndim
                 gg(idim)=m%grid(igrid)%f(ind,idim)
              end do
           end if
#endif
           ! Compute total mass
           mass_loc=mass_loc+uu(1)*vol

           ! Compute total energy
           ekin_loc=ekin_loc+uu(ndim+2)*vol

           ! Compute total internal energy
           eint_loc=eint_loc+uu(ndim+2)*vol
           do ivar=1,ndim
              eint_loc=eint_loc-0.5d0*uu(1+ivar)**2/max(uu(1),r%smallr)*vol
           end do
#if NENER>0
           do ivar=1,nener
              eint_loc=eint_loc-uu(ndim+2+ivar)*vol
           end do
#endif
           ! Compute CFL time-step
           call cmpdt_2(r,uu,gg,dx,dt_lev)
           dt_loc=min(dt_loc,dt_lev)
        endif

     end do
     ! End loop over cells
  end do
  ! End loop over grids

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
  dt_all=dt_loc
#endif
  g%mass_tot=g%mass_tot+mass_all
  g%ekin_tot=g%ekin_tot+ekin_all
  g%eint_tot=g%eint_tot+eint_all
  g%dtnew(ilevel)=MIN(g%dtnew(ilevel),dt_all)

#endif

111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine_2
!###########################################################
!###########################################################
!###########################################################
!###########################################################

