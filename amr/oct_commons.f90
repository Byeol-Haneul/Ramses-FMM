module oct_commons
  use amr_parameters
  use hydro_parameters
  use rt_parameters, only: nrtvar
  
  ! New type for oct structure
  type oct
#ifdef HYDRO
     real(kind=dp),dimension(1:twotondim,1:nvar)::uold
     real(kind=dp),dimension(1:twotondim,1:nvar)::unew
#endif
#ifdef MHD
     real(kind=dp),dimension(1:twotondim,1:6)::bold
     real(kind=dp),dimension(1:twotondim,1:6)::bnew
#endif
#ifdef RT
     real(kind=dp),dimension(1:twotondim,1:nrtvar)::rtuold
     real(kind=dp),dimension(1:twotondim,1:nrtvar)::rtunew
     real(kind=dp),dimension(1:twotondim,1:nrtgrp)::emissivity
#endif
#ifdef GRAV
     real(kind=dp),dimension(1:twotondim,1:3)::f
     real(kind=dp),dimension(1:twotondim)::rho
     real(kind=dp),dimension(1:twotondim)::phi
     real(kind=dp),dimension(1:twotondim)::phi_old
     real(kind=dp),dimension(1:twotondim)::nref
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

