cmake_minimum_required(VERSION 3.1)
project(ramses3d LANGUAGES Fortran CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Release build with debug info selected")
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

add_executable(${PROJECT_NAME})

find_package(MPI REQUIRED)      # MPI support
include_directories(${MPI_Fortran_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} ${MPI_Fortran_LIBRARIES})

if(EXISTS ${CMAKE_SOURCE_DIR}/mdl2)
  set(MDL2 true CACHE BOOL "Enable the new MDL version of RAMSES")
else()
  set(MDL2 false CACHE BOOL "Enable the new MDL version of RAMSES")
endif()

set(NVECTOR 32)
set(NDIM 3)
set(NHILBERT 1)
set(NPRE 8)
set(NVAR 5)
set(NENER 0)
set(SOLVER hydro)
set(UNITS coeur)
set(CONDINIT coeur)
set(EXEC ramses)

add_definitions(-DNVECTOR=${NVECTOR} -DNDIM=${NDIM} -DNPRE=${NPRE} -DNENER=${NENER} -DNVAR=${NVAR} -DNHILBERT=${NHILBERT} -DUNITS=${UNITS} -DCONDINIT=${CONDINIT} -DHYDRO -DGRAV)

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

set(fortran_files
    amr/ramses.f90
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
	amr/units.f90
	amr/refine_utils.f90
	amr/nbors_utils.f90
	amr/nbors_utils_p.f90
	amr/flag_utils.f90
	amr/smooth.f90
	amr/amr_step.f90
	amr/output_amr.f90
	amr/load_balance.f90
	amr/movie.f90
	amr/task_manager.f90
        amr/marshal.f90
	amr/amr_parameters.f90
	hydro/hydro_parameters.f90
	amr/hash.f90
	pm/pm_parameters.f90
	pm/pm_commons.f90
	amr/hilbert.f90
	amr/domain_hilbert.f90
	hydro/hydro_commons.f90
	amr/oct_commons.f90
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
        hydro/condinit.f90
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
	poisson/poisson_flag.f90
	poisson/multigrid_fine_commons.f90
	poisson/multigrid_fine_coarse.f90

	${CMAKE_BINARY_DIR}/write_gitinfo.f90
	${CMAKE_BINARY_DIR}/write_makefile.f90
	${CMAKE_BINARY_DIR}/write_patch.f90
)
if(MDL2)
  add_subdirectory(mdl2)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
  )
  list(APPEND fortran_files mdl2/mdl.f90 amr/mdlcache.f90 amr/nbors_utils_p.f90)
  target_link_libraries(${PROJECT_NAME} mdl2)
  target_sources(${PROJECT_NAME} PRIVATE amr/mdlcache.cxx)
  add_definitions(-DMDL2)
else()
  list(APPEND fortran_files mdl1/mdl.f90)
endif(MDL2)

set_source_files_properties(${fortran_files} ${MPI_Fortran_COMPILE_FLAGS})

# The Fortran source contains preprocessor directives which are enable differently for each compiler
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
 set_source_files_properties( ${fortran_files} PROPERTIES COMPILE_FLAGS "-cpp -recursive" )
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
 set_source_files_properties( ${fortran_files} PROPERTIES COMPILE_FLAGS "-x f95-cpp-input -frecursive -frecord-marker=4 -ffree-line-length-none -fallow-argument-mismatch" )
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
 set_source_files_properties( ${fortran_files} PROPERTIES COMPILE_FLAGS "-e ZR" )
endif()

target_sources(${PROJECT_NAME} PRIVATE ${fortran_files})

