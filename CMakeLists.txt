cmake_minimum_required(VERSION 3.18)

project(FDSim3D
  VERSION 1.0.0
  DESCRIPTION "3D finite-difference simulation of seismic wave propagation in heterogeneous viscoelastic media"
  LANGUAGES Fortran)

# --------------------------------------------------------------------------- options
option(FDSIM_ENABLE_MPI  "Build FDMod3D and FDSim3D with MPI (-DUSE_MPI)" OFF)
option(FDSIM_MPI_MODULE  "Use 'USE MPI' instead of INCLUDE 'mpif.h' (-DMPI2)" ON)
set(FDMOD3D_MOD_FUNC "${PROJECT_SOURCE_DIR}/src/fdmod3d/mod_func.F90" CACHE FILEPATH
    "Model-definition module (mod_func.F90) compiled into FDMod3D")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# --------------------------------------------------------------------------- compiler flags
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # lanczos.F90 uses a cpp line splice ("\") that exceeds 132 columns after preprocessing
  set(FDSIM_FFLAGS -ffree-line-length-none)
  set(FDSIM_FFLAGS_RELEASE -O3)
  set(FDSIM_FFLAGS_DEBUG -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")   # Intel (ifort) and IntelLLVM (ifx)
  set(FDSIM_FFLAGS "")
  set(FDSIM_FFLAGS_RELEASE -O3)
  set(FDSIM_FFLAGS_DEBUG -O0 -g -check all -traceback -fpe0)
endif()

function(fdsim_set_flags target)
  target_compile_options(${target} PRIVATE
    ${FDSIM_FFLAGS}
    $<$<CONFIG:Release>:${FDSIM_FFLAGS_RELEASE}>
    $<$<CONFIG:Debug>:${FDSIM_FFLAGS_DEBUG}>)
  # each program has its own modules (some share names), keep them apart
  set_target_properties(${target} PROPERTIES
    Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod/${target})
endfunction()

if(FDSIM_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
endif()

function(fdsim_enable_mpi target)
  if(FDSIM_ENABLE_MPI)
    target_link_libraries(${target} PRIVATE MPI::MPI_Fortran)
    target_compile_definitions(${target} PRIVATE USE_MPI)
    if(FDSIM_MPI_MODULE)
      target_compile_definitions(${target} PRIVATE MPI2)
    endif()
  endif()
endfunction()

# --------------------------------------------------------------------------- FDMod3D
file(GLOB FDMOD3D_SRC CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/fdmod3d/*.F90)
list(FILTER FDMOD3D_SRC EXCLUDE REGEX "/mod_func\\.F90$")
list(APPEND FDMOD3D_SRC ${FDMOD3D_MOD_FUNC})
set_source_files_properties(${FDMOD3D_MOD_FUNC} PROPERTIES Fortran_PREPROCESS ON)
add_executable(FDMod3D ${FDMOD3D_SRC})
fdsim_set_flags(FDMod3D)
fdsim_enable_mpi(FDMod3D)

# --------------------------------------------------------------------------- FDSim3D
file(GLOB FDSIM3D_SRC CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/fdsim3d/*.F90)
add_executable(FDSim3D ${FDSIM3D_SRC})
fdsim_set_flags(FDSim3D)
fdsim_enable_mpi(FDSim3D)

# --------------------------------------------------------------------------- FDModPar3D (serial)
add_executable(FDModPar3D ${PROJECT_SOURCE_DIR}/src/fdmodpar3d/FDModPar3D.f90)
fdsim_set_flags(FDModPar3D)

# --------------------------------------------------------------------------- SourceTF (serial)
add_executable(SourceTF
  ${PROJECT_SOURCE_DIR}/src/sourcetf/fcoolr.f90
  ${PROJECT_SOURCE_DIR}/src/sourcetf/mod_prec.f90
  ${PROJECT_SOURCE_DIR}/src/sourcetf/sourcetf.f90)
fdsim_set_flags(SourceTF)

# --------------------------------------------------------------------------- install
include(GNUInstallDirs)
install(TARGETS FDMod3D FDSim3D FDModPar3D SourceTF
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

message(STATUS "FDSim3D: MPI=${FDSIM_ENABLE_MPI}  build type=${CMAKE_BUILD_TYPE}  mod_func=${FDMOD3D_MOD_FUNC}")
