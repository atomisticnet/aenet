#
# CMake build system for aenet
#
cmake_minimum_required(VERSION 3.15)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  message(FATAL_ERROR "Please run cmake from the project root, not from src/")
endif()

# ==============================================================================
# Build Options
# ==============================================================================

option(USE_MPI "Enable MPI parallelization" OFF)
option(USE_OPENBLAS "Use OpenBLAS instead of system BLAS/LAPACK" OFF)
option(USE_MKL "Use Intel MKL instead of system BLAS/LAPACK" OFF)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (Release/Debug)")

# ==============================================================================
# Construct Suffix for Executables
# ==============================================================================

set(AENET_SUFFIX "")
if(USE_MPI)
    set(AENET_SUFFIX "${AENET_SUFFIX}_mpi")
endif()

if(USE_MKL)
    set(AENET_SUFFIX "${AENET_SUFFIX}_mkl")
elseif(USE_OPENBLAS)
    set(AENET_SUFFIX "${AENET_SUFFIX}_openblas")
else()
    # Indicate default BLAS/LAPACK if MPI is also default
    if(NOT USE_MPI)
      set(AENET_SUFFIX "${AENET_SUFFIX}_serial")
    endif()
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(AENET_SUFFIX "${AENET_SUFFIX}_debug")
endif()

# ==============================================================================
# Directories in the build tree
# ==============================================================================

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# ==============================================================================
# Compiler and Fortran Settings
# ==============================================================================

set(CMAKE_Fortran_STANDARD 2008)

# Set compiler flags based on build type
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -fexternal-blas -pedantic")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -fexternal-blas")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -check all -traceback -fpe0")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")
endif()

# ==============================================================================
# Find Dependencies
# ==============================================================================

# Hint BLAS/LAPACK vendor based on user options
# (mkls: Intel10_64lp = threaded, Intel10_64lp_seq = sequential)
set(BLA_VENDOR "" CACHE STRING "Preferred BLAS/LAPACK vendor (e.g., OpenBLAS, Intel10_64lp_seq)")

if(USE_MKL)
  set(BLA_VENDOR "Intel10_64lp_seq" CACHE STRING "" FORCE)
elseif(USE_OPENBLAS)
  set(BLA_VENDOR "OpenBLAS" CACHE STRING "" FORCE)
endif()

# If vendor hint changes, clear cached results so CMake re-finds libs
unset(BLAS_LIBRARIES CACHE)
unset(LAPACK_LIBRARIES CACHE)

# Prefer static linkage
set(BLA_STATIC ON)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Find MPI if requested
if(USE_MPI)
  find_package(MPI REQUIRED)
endif()

# ==============================================================================
# L-BFGS-B Library Integration (extracted from tarball)
# ==============================================================================

set(LBFGSB_PARENT_DIR "${PROJECT_BINARY_DIR}/_third_party/lbfgsb")
file(MAKE_DIRECTORY "${LBFGSB_PARENT_DIR}")
set(LBFGSB_DIR "${LBFGSB_PARENT_DIR}/Lbfgsb.3.0")
set(LBFGSB_TARBALL "${PROJECT_SOURCE_DIR}/lib/Lbfgsb.3.0.tar.gz")
set(LBFGSB_SOURCES
    "${LBFGSB_DIR}/blas.f"
    "${LBFGSB_DIR}/lbfgsb.f"
    "${LBFGSB_DIR}/linpack.f"
    "${LBFGSB_DIR}/timer.f"
)

# Custom command to extract L-BFGS-B sources during build
add_custom_command(
    OUTPUT ${LBFGSB_SOURCES}
    COMMAND ${CMAKE_COMMAND} -E tar xzf "${LBFGSB_TARBALL}"
    WORKING_DIRECTORY ${LBFGSB_PARENT_DIR}
    DEPENDS ${LBFGSB_TARBALL}
    COMMENT "Extracting L-BFGS-B sources from ${LBFGSB_TARBALL}"
    VERBATIM
)

# Add the library target and its sources
add_library(lbfgsb_static STATIC)
target_sources(lbfgsb_static PRIVATE ${LBFGSB_SOURCES})
set_target_properties(lbfgsb_static PROPERTIES
    OUTPUT_NAME "lbfgsb"
    POSITION_INDEPENDENT_CODE ON
)

# ==============================================================================
# Source File Definitions
# ==============================================================================

set(SRC_FILES
    aenet.f90 aeio.f90 constants.f90
    geometry.f90 input.f90 optimize.f90 parallel.F90
    potential.f90 random.f90 sfsetup.f90
    trainset.f90 pytorchoutput.f90
)

set(EXT_FILES
    ext/arglib.f90 ext/chebyshev.f90 ext/feedforward.f90 ext/io.f90
    ext/lclist.f90 ext/sfbasis.f90 ext/sortlib.f90 ext/symmfunc.f90
    ext/timing.f90 ext/unittest.f90 ext/xsflib.f90
)

# ==============================================================================
# Library Targets (aenet static and shared)
# ==============================================================================

add_library(aenet_static STATIC ${SRC_FILES} ${EXT_FILES})
add_library(aenet_shared SHARED ${SRC_FILES} ${EXT_FILES})

set_target_properties(aenet_static PROPERTIES OUTPUT_NAME "aenet")
set_target_properties(aenet_shared PROPERTIES OUTPUT_NAME "aenet")

# Link libraries
target_link_libraries(aenet_static PRIVATE lbfgsb_static ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(aenet_shared PRIVATE lbfgsb_static ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
if(USE_MPI)
    target_link_libraries(aenet_static PRIVATE MPI::MPI_Fortran)
    target_link_libraries(aenet_shared PRIVATE MPI::MPI_Fortran)
endif()

# ==============================================================================
# Executable Targets
# ==============================================================================

# --- Main Executables ---
add_executable(generate.x generate.f90)
add_executable(train.x train.f90)
add_executable(predict.x predict.F90)

target_link_libraries(generate.x PRIVATE aenet_static)
target_link_libraries(train.x PRIVATE aenet_static)
target_link_libraries(predict.x PRIVATE aenet_static)

set_target_properties(generate.x PROPERTIES OUTPUT_NAME generate.x${AENET_SUFFIX})
set_target_properties(train.x PROPERTIES OUTPUT_NAME train.x${AENET_SUFFIX})
set_target_properties(predict.x PROPERTIES OUTPUT_NAME predict.x${AENET_SUFFIX})

# --- Tool Executables ---
add_executable(fingerprint.x tools/fingerprint.f90)
add_executable(neighbors.x tools/neighbors.f90)
add_executable(trnset_info.x tools/trnset_info.f90)
add_executable(trnset2ASCII.x tools/trnset2ASCII.f90)

target_link_libraries(fingerprint.x PRIVATE aenet_static)
target_link_libraries(neighbors.x PRIVATE aenet_static)
target_link_libraries(trnset_info.x PRIVATE aenet_static)
target_link_libraries(trnset2ASCII.x PRIVATE aenet_static)

set_target_properties(fingerprint.x neighbors.x trnset_info.x trnset2ASCII.x PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tools)

# --- Test Executables ---
enable_testing()
add_executable(test_feedforward tests/test_feedforward.f90)
add_executable(test_geometry tests/test_geometry.f90)
add_executable(test_lclist tests/test_lclist.f90)
add_executable(test_sfbasis tests/test_sfbasis.f90)
add_executable(test_sfsetup tests/test_sfsetup.f90)
add_executable(test_symmfunc tests/test_symmfunc.f90)
add_executable(test_xsflib tests/test_xsflib.f90)

target_link_libraries(test_feedforward PRIVATE aenet_static)
target_link_libraries(test_geometry PRIVATE aenet_static)
target_link_libraries(test_lclist PRIVATE aenet_static)
target_link_libraries(test_sfbasis PRIVATE aenet_static)
target_link_libraries(test_sfsetup PRIVATE aenet_static)
target_link_libraries(test_symmfunc PRIVATE aenet_static)
target_link_libraries(test_xsflib PRIVATE aenet_static)

set_target_properties(
  test_feedforward test_geometry test_lclist test_sfbasis
  test_sfsetup test_symmfunc test_xsflib
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/unittests
)

add_test(NAME feedforward:IO COMMAND test_feedforward test_IO)
add_test(NAME feedforward:eval COMMAND test_feedforward test_eval)
add_test(NAME feedforward:activation COMMAND test_feedforward test_activation)
add_test(NAME geometry:isolated COMMAND test_geometry test_isolated)
add_test(NAME geometry:periodic COMMAND test_geometry test_periodic)
add_test(NAME lclist:fcc_primitive COMMAND test_lclist test_fcc_primitive)
add_test(NAME lclist:fcc_conventional COMMAND test_lclist test_fcc_conventional)
add_test(NAME lclist:fcc_isolated COMMAND test_lclist test_fcc_isolated)
add_test(NAME sfbasis:setup COMMAND test_sfbasis test_setup)
add_test(NAME sfbasis:fcc_rotation COMMAND test_sfbasis test_fcc_rotation)
add_test(NAME sfbasis:derivatives COMMAND test_sfbasis test_derivatives)
add_test(NAME sfsetup:IO_chebyshev COMMAND test_sfsetup test_IO_chebyshev)
add_test(NAME sfsetup:IO_behler2011 COMMAND test_sfsetup test_IO_behler2011)
add_test(NAME symmfunc:setup COMMAND test_symmfunc test_setup)
add_test(NAME symmfunc:fcc_rotation COMMAND test_symmfunc test_fcc_rotation)
add_test(NAME symmfunc:derivatives COMMAND test_symmfunc test_derivatives)
add_test(NAME xsflib:isolated COMMAND test_xsflib test_isolated)
add_test(NAME xsflib:periodic COMMAND test_xsflib test_periodic)

# ==============================================================================
# Custom Targets (for convenience)
# ==============================================================================

add_custom_target(main DEPENDS generate.x train.x predict.x)
add_custom_target(lib DEPENDS aenet_static aenet_shared)
add_custom_target(tools DEPENDS fingerprint.x neighbors.x trnset_info.x trnset2ASCII.x)
add_custom_target(build_tests
   DEPENDS test_feedforward test_geometry test_lclist test_sfbasis
           test_sfsetup test_symmfunc test_xsflib)
add_custom_target(build_all DEPENDS main lib tools build_tests)

add_custom_target(list_targets ALL
   COMMAND ${CMAKE_COMMAND} -E echo ""
   COMMAND ${CMAKE_COMMAND} -E echo "Build with"
   COMMAND ${CMAKE_COMMAND} -E echo ""
   COMMAND ${CMAKE_COMMAND} -E echo "  cmake --build . --target <target>"
   COMMAND ${CMAKE_COMMAND} -E echo "  (OR make <target>)"
   COMMAND ${CMAKE_COMMAND} -E echo ""
   COMMAND ${CMAKE_COMMAND} -E echo "where <target> is one of the following:"
   COMMAND ${CMAKE_COMMAND} -E echo ""
   COMMAND ${CMAKE_COMMAND} -E echo "  main         (Builds generate.x, train.x, predict.x)"
   COMMAND ${CMAKE_COMMAND} -E echo "  lib          (Builds aenet static and shared libraries)"
   COMMAND ${CMAKE_COMMAND} -E echo "  tools        (Builds all tool executables)"
   COMMAND ${CMAKE_COMMAND} -E echo "  build_all    (Builds all executables, libraries, and tests)"
   COMMAND ${CMAKE_COMMAND} -E echo "  build_tests  (Builds all unit test executables)"
   COMMAND ${CMAKE_COMMAND} -E echo "  test         (Runs the unit tests)"
   COMMAND ${CMAKE_COMMAND} -E echo "  list_targets (Prints this message)"
   COMMAND ${CMAKE_COMMAND} -E echo ""
   VERBATIM)

# exclude all targets from the default ALL except for the list_target target
set_property(TARGET generate.x train.x predict.x
                    lbfgsb_static aenet_static aenet_shared
                    fingerprint.x neighbors.x trnset_info.x trnset2ASCII.x
                    test_feedforward test_geometry test_lclist test_sfbasis
                    test_sfsetup test_symmfunc test_xsflib
             PROPERTY EXCLUDE_FROM_ALL TRUE)

# add a "test" that builds all test binaries
add_test(test_build_test_binaries
  "${CMAKE_COMMAND}"
  --build "${CMAKE_BINARY_DIR}"
  --config "$<CONFIG>"
  --target build_tests)

# Make building the test binaries a required fixture for all tests.
# This ensures all test binaries are build when 'make test' is run.
set_tests_properties(test_build_test_binaries PROPERTIES FIXTURES_SETUP test_fixture)
set_tests_properties(feedforward:IO feedforward:eval feedforward:activation
                     geometry:isolated geometry:periodic
                     lclist:fcc_primitive lclist:fcc_conventional lclist:fcc_isolated
                     sfbasis:setup sfbasis:fcc_rotation sfbasis:derivatives
                     sfsetup:IO_chebyshev sfsetup:IO_behler2011
                     symmfunc:setup symmfunc:fcc_rotation symmfunc:derivatives
                     xsflib:isolated xsflib:periodic
  PROPERTIES FIXTURES_REQUIRED test_fixture)

# ==============================================================================
# Installation Rules
# ==============================================================================

# Install main executables with suffix
install(TARGETS generate.x train.x predict.x
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)

# Install libraries
install(TARGETS aenet_static aenet_shared
        ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib
        LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib)

# Install tool executables without suffix
install(TARGETS fingerprint.x neighbors.x trnset_info.x trnset2ASCII.x
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/tools)
