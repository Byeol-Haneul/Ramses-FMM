module ramses_commons
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_module, only: mdl_t
  use cooling_module, only: cooling_t

  type ramses_t

     type(run_t)::r
     type(global_t)::g
     type(mesh_t)::m
     type(part_t)::p
     type(part_t)::s
     type(cooling_t)::cool
     type(mdl_t),pointer::mdl => null()

  end type ramses_t

  type pst_t
     
     type(ramses_t),pointer::s => null()
     type(pst_t),pointer::pLower => null()
     integer::iUpper = -1
     integer::nLower = 0
     integer::nUpper = 0
     
  end type pst_t

  
end module ramses_commons

