module oct_commons
  use amr_parameters
  use hydro_parameters
  
  ! New type for oct structure
  type oct
     integer(kind=4)::lev
     integer(kind=4),dimension(1:ndim)::ckey
     integer(kind=8),dimension(1:nhilbert)::hkey
     integer(kind=4),dimension(1:twotondim)::flag1
     integer(kind=4),dimension(1:twotondim)::flag2
     logical,dimension(1:twotondim)::refined
     integer(kind=4)::superoct
#ifdef GRAV
     real(kind=dp),dimension(1:twotondim)::rho
     real(kind=dp),dimension(1:twotondim)::phi
     real(kind=dp),dimension(1:twotondim)::phi_old
     real(kind=dp),dimension(1:twotondim,1:3)::f
#endif
#ifdef HYDRO
     real(kind=dp),dimension(1:twotondim,1:nvar)::uold
     real(kind=dp),dimension(1:twotondim,1:nvar)::unew
#endif
#ifdef DUALENER
     real(kind=dp),dimension(1:twotondim)::divu
     real(kind=dp),dimension(1:twotondim)::enew
#endif
  end type oct

  type nbor
     type(oct),pointer::p
  end type nbor
  
end module oct_commons

