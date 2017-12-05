module call_back
  
  interface
     subroutine ramses_function(r,g,m,p,mdl,cpu_range,input_size,output_size,input,output)
       use amr_commons, only: run_t,global_t,mesh_t
       use pm_commons, only: part_t
       use mdl_commons, only: mdl_t
       type(run_t)::r
       type(global_t)::g
       type(mesh_t)::m
       type(part_t)::p
       type(mdl_t)::mdl
       integer::cpu_range,input_size,output_size
       integer,dimension(1:input_size),optional::input
       integer,dimension(1:output_size),optional::output
     end subroutine ramses_function
  end interface

  procedure(ramses_function)::r_clean_stop,r_broadcast_params,r_broadcast_global
  procedure(ramses_function)::r_init_amr,r_init_time,r_init_hydro,r_init_part
  procedure(ramses_function)::r_input_part_grafic,r_input_part_ascii,r_input_part_restart
  procedure(ramses_function)::r_init_flag,r_user_flag,r_ensure_ref_rules
  procedure(ramses_function)::r_collect_noct,r_noct_tot,r_noct_min,r_noct_max,r_gather_noct_max
  procedure(ramses_function)::r_init_refine_basegrid,r_init_refine_restart
  procedure(ramses_function)::r_broadcast_bound_key,r_collect_bound_key,r_load_balance
  procedure(ramses_function)::r_refine_fine,r_smooth_fine
  procedure(ramses_function)::r_input_hydro_condinit,r_input_hydro_grafic
  procedure(ramses_function)::r_upload_fine
  procedure(ramses_function)::r_multipole_leaf_cells,r_multipole_split_cells
  procedure(ramses_function)::r_reset_rho,r_cic_multipole,r_cic_part,r_split_part
  procedure(ramses_function)::r_collect_multipole,r_broadcast_multipole
  procedure(ramses_function)::r_output_amr,r_output_hydro,r_output_poisson,r_output_part
  procedure(ramses_function)::r_synchro_hydro_fine,r_force_analytic,r_gradient_phi
  procedure(ramses_function)::r_save_phi_old,r_compute_epot,r_compute_rhomax
  procedure(ramses_function)::r_broadcast_aexp,r_courant_fine,r_godunov_fine
  procedure(ramses_function)::r_set_unew,r_set_uold,r_gravity_hydro_fine
  procedure(ramses_function)::r_cooling_fine,r_newdt_part,r_broadcast_dt
  procedure(ramses_function)::r_make_initial_phi,r_init_mg,r_build_mg,r_cleanup_mg
  procedure(ramses_function)::r_make_mask,r_make_bc_rhs,r_restrict_mask
  procedure(ramses_function)::r_cmp_residual_mg,r_gauss_seidel_mg,r_reset_correction
  procedure(ramses_function)::r_restrict_residual,r_interpolate_and_correct
  procedure(ramses_function)::r_set_scan_flag,r_cmp_residual_norm2
  procedure(ramses_function)::r_output_frame
  
  type call_back_f
     procedure(ramses_function),pointer,nopass::proc
  end type call_back_f

end module call_back
