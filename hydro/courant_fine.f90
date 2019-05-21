module courant_fine_module

type :: out_courant_fine_t
  real(kind=8)::mass,ekin,eint,dt
end type out_courant_fine_t

contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_courant_fine(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::input_size,output_size
  type(out_courant_fine_t) output, next_output

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_COURANT_FINE,pst%iUpper+1,input_size,output_size,ilevel)
     call r_courant_fine(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
     output%ekin=output%ekin+next_output%ekin
     output%eint=output%eint+next_output%eint
     output%dt=MIN(output%dt,next_output%dt)
  else
     call courant_fine(pst%s%r,pst%s%g,pst%s%m,ilevel,output%mass,output%ekin,output%eint,output%dt)
  endif

end subroutine r_courant_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine courant_fine(r,g,m,ilevel,mass,ekin,eint,dt)
  use amr_parameters, only: dp,nvector,ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  real(kind=8)::mass,ekin,eint,dt
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::ivar,idim,ind,igrid
  real(dp)::dt_lev,dx,vol
  real(dp),dimension(1:nvar)::uu
  real(dp),dimension(1:ndim)::gg

#ifdef HYDRO

  ! Mesh spacing at that level
  dx=0.5D0**ilevel*r%boxlen
  vol=dx**ndim

  mass=0.0D0
  ekin=0.0D0
  eint=0.0D0
  dt=dx/r%smallc

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
           mass=mass+uu(1)*vol

           ! Compute total energy
           ekin=ekin+uu(ndim+2)*vol

           ! Compute total internal energy
           eint=eint+uu(ndim+2)*vol
           do ivar=1,ndim
              eint=eint-0.5D0*uu(1+ivar)**2/max(uu(1),r%smallr)*vol
           end do
#if NENER>0
           do ivar=1,nener
              eint=eint-uu(ndim+2+ivar)*vol
           end do
#endif
           ! Compute CFL time-step
           call cmpdt(r,uu,gg,dx,dt_lev)
           dt=min(dt,dt_lev)
        endif

     end do
     ! End loop over cells
  end do
  ! End loop over grids

#endif

end subroutine courant_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module courant_fine_module