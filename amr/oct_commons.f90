module oct_commons
  use amr_parameters
  use hydro_parameters
  use rt_parameters, only: nrtvar, nrtgrp
  
  ! New type for oct structure
  type oct
#ifdef HYDRO
     real(dp),dimension(1:twotondim,1:nvar)::uold
     real(dp),dimension(1:twotondim,1:nvar)::unew
#endif
#if defined(HYDRO) && defined(GRADVPART)
     ! Per-cell PLM velocity slopes (dq) for (u,v,w) along (x,y,z)
     ! gradv(icell, ivel, idim) with ivel=1..NDIM, idim=1..NDIM
     real(dp),dimension(1:twotondim,1:ndim,1:ndim)::gradv
#endif
#ifdef MHD
     real(dp),dimension(1:twotondim,1:6)::bold
     real(dp),dimension(1:twotondim,1:6)::bnew
#endif
#ifdef RT
     real(dp),dimension(1:twotondim,1:nrtvar)::rtuold
     real(dp),dimension(1:twotondim,1:nrtvar)::rtunew
     real(dp),dimension(1:twotondim,1:nrtgrp)::emissivity
#endif
#ifdef TURB
     real(dp),dimension(1:twotondim,1:3)::fturb
#endif
#ifdef GRAV
     real(dp),dimension(1:twotondim,1:3)::f
     real(dp),dimension(1:twotondim)::rho
     real(dp),dimension(1:twotondim)::phi
     real(dp),dimension(1:twotondim)::phi_old
     real(dp),dimension(1:twotondim)::nref
#endif
     integer(kind=8),dimension(1:nhilbert)::hkey
     integer(kind=4),dimension(1:twotondim)::flag1
     integer(kind=4),dimension(1:twotondim)::flag2
     integer(kind=4),dimension(1:ndim)::ckey
     logical,dimension(1:twotondim)::refined
     integer(kind=4)::lev
     integer(kind=4)::superoct
  end type oct

  type nbor
     type(oct),pointer::p
  end type nbor
  
end module oct_commons

