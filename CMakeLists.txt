cmake_minimum_required(VERSION 3.1)
project(ramses3d)
enable_language (Fortran)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "ifort.*")
  set (CMAKE_Fortran_FLAGS "-cpp")
endif()
add_executable(${PROJECT_NAME} "")

find_package(MPI REQUIRED)      # MPI support
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${MPI_Fortran_LIBRARIES})

if(EXISTS ${CMAKE_SOURCE_DIR}/mdl2)
  set(MDL2 true CACHE BOOL "Enable the new MDL version of RAMSES")
else()
  set(MDL2 false CACHE BOOL "Enable the new MDL version of RAMSES")
endif()
if(MDL2)
  add_subdirectory(mdl2)
  target_link_libraries(${PROJECT_NAME} mdl2 openpa)
  target_sources(${PROJECT_NAME} PRIVATE mdl2/mdl.f90)
  add_definitions(-DMDL2)
endif(MDL2)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set_target_properties( ${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-x f95-cpp-input -frecord-marker=4 -ffree-line-length-none -fcheck=all -g" )
endif()

#if(Git_FOUND)
#  message("Git found: ${GIT_EXECUTABLE}")
#endif()

set(NVECTOR 32)
set(NDIM 3)
set(NHILBERT 1)
set(NPRE 8)
set(NVAR 5)
set(NENER 0)
set(SOLVER hydro)
set(EXEC ramses)
#add_definitions(-DNVECTOR=${NVECTOR} -DNDIM=${NDIM} -DNPRE=${NPRE} -DNENER=${NENER} -DNVAR=${NVAR} -DNHILBERT=${NHILBERT} -DHYDRO -DGRAV -DMDL2)
add_definitions(-DNVECTOR=${NVECTOR} -DNDIM=${NDIM} -DNPRE=${NPRE} -DNENER=${NENER} -DNVAR=${NVAR} -DNHILBERT=${NHILBERT} -DHYDRO -DGRAV)

add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/write_gitinfo.f90
    COMMAND ${CMAKE_SOURCE_DIR}/utils/scripts/cr_write_gitinfo.sh ${CMAKE_BINARY_DIR}/write_gitinfo.f90
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    VERBATIM
)

add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/write_makefile.f90
    COMMAND ${CMAKE_SOURCE_DIR}/utils/scripts/cr_write_makefile.sh ${CMAKE_SOURCE_DIR}/CMakeLists.txt ${CMAKE_BINARY_DIR}/write_makefile.f90
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    VERBATIM
)

add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/write_patch.f90
    COMMAND ${CMAKE_SOURCE_DIR}/utils/scripts/cr_write_patch.sh "" ${CMAKE_BINARY_DIR}/write_patch.f90
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    VERBATIM
)

target_sources(${PROJECT_NAME}
PRIVATE	amr/ramses.f90
	amr/timer.f90
	amr/read_params.f90
	amr/init_amr.f90
	amr/init_time.f90
	amr/init_refine_basegrid.f90
	amr/init_refine_adaptive.f90
	amr/init_refine_restart.f90
	amr/adaptive_loop.f90
	amr/update_time.f90
	amr/title.f90
	patch/coeur/units.f90
	amr/refine_utils.f90
	amr/nbors_utils.f90
	amr/flag_utils.f90
	amr/smooth.f90
	amr/amr_step.f90
	amr/output_amr.f90
	amr/load_balance.f90
	amr/movie.f90
	amr/task_manager.f90

	amr/amr_parameters.f90
	hydro/hydro_parameters.f90
	amr/hash.f90
	pm/pm_parameters.f90
	pm/pm_commons.f90
	amr/hilbert.f90
	amr/domain_hilbert.f90
	hydro/hydro_commons.f90
	amr/amr_commons.f90
	poisson/poisson_parameters.f90
	poisson/poisson_commons.f90
	amr/mdl_commons.f90
	amr/ramses_commons.f90
	amr/call_back.f90
	amr/cache_commons.f90
	pm/gadgetreadfile.f90

        hydro/init_hydro.f90
	hydro/init_flow_fine.f90
	hydro/input_hydro_condinit.f90
	hydro/input_hydro_grafic.f90
	hydro/interpol_hydro.f90
        patch/coeur/condinit.f90
	hydro/hydro_flag.f90
	hydro/godunov_utils.f90
	hydro/write_screen.f90
	hydro/courant_fine.f90
        hydro/godunov_fine.f90
	hydro/umuscl.f90
	hydro/uplmde.f90
	hydro/output_hydro.f90
	hydro/synchro_hydro_fine.f90
	hydro/cooling_fine.f90

	pm/newdt_fine.f90
	pm/init_part.f90
	pm/input_part.f90
	pm/input_part_restart.f90
	pm/input_part_ascii.f90
	pm/input_part_grafic.f90
	pm/rho_fine.f90
	pm/move_fine.f90
	poisson/rho_ana.f90
	pm/output_part.f90

	poisson/output_poisson.f90
	poisson/force_fine.f90
	poisson/grav_ana.f90
	poisson/phi_fine_cg.f90
	hydro/upload.f90
	poisson/interpol_phi.f90
	poisson/poisson_flag.f90
	poisson/multigrid_fine_commons.f90
	poisson/multigrid_fine_coarse.f90

	${CMAKE_BINARY_DIR}/write_gitinfo.f90
	${CMAKE_BINARY_DIR}/write_makefile.f90
	${CMAKE_BINARY_DIR}/write_patch.f90

)

