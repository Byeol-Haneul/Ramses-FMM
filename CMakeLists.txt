cmake_minimum_required(VERSION 3.5)
cmake_policy(VERSION 3.5...3.14)

include(CMakePrintHelpers)

project(ramses LANGUAGES Fortran CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Release build with debug info selected")
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

set(NVECTOR 32)
set(NHILBERT 1)
set(NDIM 2 CACHE STRING "Set the number of dimensions")
set(EXEC ramses${NDIM}d)
set(NPRE 8)
set(NVAR 5)
set(NENER 0)
set(UNITS DEF CACHE STRING "Set the unit system")
set(INIT DEF CACHE STRING "Set the initial conditions")
set(GRAV false CACHE BOOL "Activate gravity solver")
set(HYDRO true CACHE BOOL "Activate hydro solver")
set(MHD false CACHE BOOL "Activate MHD solver")

add_executable(${EXEC} "")

find_package(MPI)      # MPI support
if(MPI_FOUND)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  target_link_libraries(${EXEC} ${MPI_Fortran_LIBRARIES})
else()
  message(STATUS "MPI library not found. Build without MPI")
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/mdl2)
  set(MDL2 true CACHE BOOL "Enable the new MDL version of RAMSES")
else()
  set(MDL2 false CACHE BOOL "Enable the new MDL version of RAMSES")
endif()

cmake_print_variables(NDIM NVAR UNITS INIT GRAV HYDRO MHD)

add_definitions(-DNVECTOR=${NVECTOR} -DNDIM=${NDIM} -DNPRE=${NPRE} -DNENER=${NENER} -DNVAR=${NVAR} -DNHILBERT=${NHILBERT} -DUNITS=${UNITS} -DINIT=${INIT})

if(NOT MPI_FOUND)
  add_definitions(-DWITHOUTMPI)
endif()

if(HYDRO)
  add_definitions(-DHYDRO)
endif()

if(GRAV)
  add_definitions(-DGRAV)
endif()

if(MHD)
  add_definitions(-DMHD)
endif()

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
  amr/flag_utils.f90
  amr/smooth.f90
  amr/amr_step.f90
  amr/output_amr.f90
  amr/load_balance.f90
  amr/movie.f90
  amr/task_manager.f90
  amr/cache.f90
  amr/sort.f90
  amr/boundaries.f90
  amr/marshal.f90
  amr/amr_parameters.f90
  hydro/hydro_parameters.f90
  pm/pm_parameters.f90
  pm/pm_commons.f90
  amr/hilbert.f90
  amr/domain_hilbert.f90
  hydro/hydro_commons.f90
  amr/oct_commons.f90
  amr/amr_commons.f90
  poisson/poisson_parameters.f90
  amr/mdl_commons.f90
  amr/ramses_commons.f90
  amr/call_back.f90
  amr/cache_commons.f90
  pm/gadgetreadfile.f90
  amr/rngstream.f90
  
  hydro/boundana.f90
  hydro/init_hydro.f90
  hydro/init_flow_fine.f90
  hydro/input_hydro_condinit.f90
  hydro/input_hydro_grafic.f90
  hydro/input_hydro_gadget.f90
  hydro/interpol_hydro.f90
  hydro/condinit.f90
  hydro/hydro_flag.f90
  hydro/godunov_utils.f90
  hydro/write_screen.f90
  hydro/courant_fine.f90
  hydro/godunov_fine.f90
  hydro/umuscl.f90
  hydro/output_hydro.f90
  hydro/synchro_hydro_fine.f90
  hydro/source_hydro_fine.f90
  
  pm/newdt_fine.f90
  pm/init_part.f90
  pm/input_part.f90
  pm/input_part_restart.f90
  pm/input_part_ascii.f90
  pm/input_part_grafic.f90
  pm/input_part_gadget.f90
  pm/input_part_zoom.f90
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

  cooling/constants.f90
  cooling/cooling_fine.f90
  cooling/cooling_module.f90
  cooling/cooling_module_ism.f90
  cooling/init_cooling.f90

  clump/clump_finder.f90
  clump/clump_merger.f90
  clump/clfind_commons.f90
  clump/sparse_mat.f90
  clump/output_clump.f90

  star/star_formation.f90
  feedback/feedback.f90

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
  list(APPEND fortran_files mdl2/mdl.f90 amr/mdlcache.f90)
  target_link_libraries(${PROJECT_NAME} mdl2)
  target_sources(${PROJECT_NAME} PRIVATE amr/mdlcache.cxx)
  add_definitions(-DMDL2)
else()
  list(APPEND fortran_files mdl1/mdl.f90 amr/hash.f90 amr/nbors_utils.f90)
endif(MDL2)

set_source_files_properties(${fortran_files} ${MPI_Fortran_COMPILE_FLAGS})

# The Fortran source contains preprocessor directives which are enable differently for each compiler
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set_source_files_properties( ${fortran_files} PROPERTIES COMPILE_FLAGS "-cpp -recursive" )
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10)
    set_source_files_properties( ${fortran_files} PROPERTIES COMPILE_FLAGS "-x f95-cpp-input -frecursive -frecord-marker=4 -ffree-line-length-none" )
  else()
   set_source_files_properties( ${fortran_files} PROPERTIES COMPILE_FLAGS "-x f95-cpp-input -frecursive -frecord-marker=4 -ffree-line-length-none -fallow-argument-mismatch" )
  endif()
  if(APPLE)
    target_link_options(${EXEC} PRIVATE -Wl,-no_compact_unwind)
  endif(APPLE)

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  set_source_files_properties( ${fortran_files} PROPERTIES COMPILE_FLAGS "-e ZR" )
endif()

target_sources(${EXEC} PRIVATE ${fortran_files})

