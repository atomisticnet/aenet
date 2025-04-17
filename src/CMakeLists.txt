cmake_minimum_required(VERSION 3.15)
project(AENET Fortran)

# Options
option(USE_MPI "Enable MPI parallelization" OFF)
option(USE_OPENBLAS "Use OpenBLAS instead of system BLAS/LAPACK" OFF)
option(USE_MKL "Use Intel MKL instead of system BLAS/LAPACK" OFF)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (Release/Debug)")

# Construct suffix based on options
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

# Fortran standard
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Set global output directories (relative to current build dir -> project root)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin) # Default runtime is bin

# Find BLAS/LAPACK
if(USE_MKL)
    find_package(MKL REQUIRED)
    set(BLAS_LIBRARIES ${MKL_LIBRARIES})
    set(LAPACK_LIBRARIES "")
elseif(USE_OPENBLAS)
    find_package(OpenBLAS REQUIRED)
    set(BLAS_LIBRARIES ${OpenBLAS_LIB})
    set(LAPACK_LIBRARIES "")
else()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif()

# --- L-BFGS-B Library Integration ---
set(LBFGSB_PARENT_DIR "${CMAKE_SOURCE_DIR}/../lib") # Directory containing the tarball
set(LBFGSB_DIR "${LBFGSB_PARENT_DIR}/Lbfgsb.3.0")   # Directory where sources will be extracted
set(LBFGSB_TARBALL "${LBFGSB_PARENT_DIR}/Lbfgsb.3.0.tar.gz")
set(LBFGSB_SOURCES
    "${LBFGSB_DIR}/blas.f"
    "${LBFGSB_DIR}/lbfgsb.f"
    "${LBFGSB_DIR}/linpack.f"
    "${LBFGSB_DIR}/timer.f"
)

# Custom command to extract L-BFGS-B sources during build
# This command declares the source files themselves as output
add_custom_command(
    OUTPUT ${LBFGSB_SOURCES}
    COMMAND ${CMAKE_COMMAND} -E tar xzf "${LBFGSB_TARBALL}"
    WORKING_DIRECTORY ${LBFGSB_PARENT_DIR}
    DEPENDS ${LBFGSB_TARBALL}
    COMMENT "Extracting L-BFGS-B sources from ${LBFGSB_TARBALL}"
    VERBATIM
)

# Add the library target initially without sources
add_library(lbfgsb_static STATIC)

# Add the generated sources to the library target
# CMake understands these sources come from the custom command above
target_sources(lbfgsb_static PRIVATE ${LBFGSB_SOURCES})

set_target_properties(lbfgsb_static PROPERTIES
    OUTPUT_NAME "lbfgsb"
    # ARCHIVE_OUTPUT_DIRECTORY inherits global setting (${CMAKE_BINARY_DIR}/../lib)
    POSITION_INDEPENDENT_CODE ON
)

# Find MPI if needed
if(USE_MPI)
    find_package(MPI REQUIRED)
endif()

# Source files
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

# Main executables
set(GENERATE_SRC generate.f90)
set(TRAIN_SRC train.f90)
set(PREDICT_SRC predict.F90)

# Tools
set(TOOLS_SRC
    tools/fingerprint.f90
    tools/neighbors.f90
    tools/trnset_info.f90
    tools/trnset2ASCII.f90
)

# Tests
set(TESTS_SRC
    tests/test_feedforward.f90
    tests/test_geometry.f90
    tests/test_lclist.f90
    tests/test_sfbasis.f90
    tests/test_sfsetup.f90
    tests/test_symmfunc.f90
    tests/test_xsflib.f90
)

# Static and shared libraries
add_library(aenet_static STATIC ${SRC_FILES} ${EXT_FILES})
add_library(aenet_shared SHARED ${SRC_FILES} ${EXT_FILES})
# ARCHIVE/LIBRARY_OUTPUT_DIRECTORY inherits global setting (${CMAKE_BINARY_DIR}/../lib)
set_target_properties(aenet_static PROPERTIES OUTPUT_NAME "aenet")
set_target_properties(aenet_shared PROPERTIES OUTPUT_NAME "aenet")

# Link libraries
target_link_libraries(aenet_static PRIVATE lbfgsb_static ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(aenet_shared PRIVATE lbfgsb_static ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
if(USE_MPI)
    target_link_libraries(aenet_static PRIVATE MPI::MPI_Fortran)
    target_link_libraries(aenet_shared PRIVATE MPI::MPI_Fortran)
endif()

# Main executables
add_executable(generate.x ${GENERATE_SRC})
add_executable(train.x ${TRAIN_SRC})
add_executable(predict.x ${PREDICT_SRC})
target_link_libraries(generate.x PRIVATE aenet_static)
target_link_libraries(train.x PRIVATE aenet_static)
target_link_libraries(predict.x PRIVATE aenet_static)

# Tool executables
add_executable(fingerprint.x tools/fingerprint.f90)
add_executable(neighbors.x tools/neighbors.f90)
add_executable(trnset_info.x tools/trnset_info.f90)
add_executable(trnset2ASCII.x tools/trnset2ASCII.f90)
target_link_libraries(fingerprint.x PRIVATE aenet_static)
target_link_libraries(neighbors.x PRIVATE aenet_static)
target_link_libraries(trnset_info.x PRIVATE aenet_static)
target_link_libraries(trnset2ASCII.x PRIVATE aenet_static)

# Test executables and CTest integration
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
add_test(NAME test_feedforward COMMAND test_feedforward)
add_test(NAME test_geometry COMMAND test_geometry)
add_test(NAME test_lclist COMMAND test_lclist)
add_test(NAME test_sfbasis COMMAND test_sfbasis)
add_test(NAME test_sfsetup COMMAND test_sfsetup)
add_test(NAME test_symmfunc COMMAND test_symmfunc)
add_test(NAME test_xsflib COMMAND test_xsflib)

# Custom targets for compatibility with legacy Makefile workflow
add_custom_target(lib DEPENDS aenet_static aenet_shared)
add_custom_target(tools DEPENDS fingerprint.x neighbors.x trnset_info.x trnset2ASCII.x)
add_custom_target(build_tests DEPENDS test_feedforward test_geometry test_lclist test_sfbasis test_sfsetup test_symmfunc test_xsflib) # Renamed from test_target
add_custom_target(main DEPENDS generate.x train.x predict.x) # Renamed from main_targets

# Install rules (optional)
# Install main executables with suffix
install(TARGETS generate.x train.x predict.x
        RUNTIME DESTINATION ../bin RENAME generate.x${AENET_SUFFIX})
install(TARGETS train.x
        RUNTIME DESTINATION ../bin RENAME train.x${AENET_SUFFIX})
install(TARGETS predict.x
        RUNTIME DESTINATION ../bin RENAME predict.x${AENET_SUFFIX})

install(TARGETS aenet_static aenet_shared
        ARCHIVE DESTINATION ../lib
        LIBRARY DESTINATION ../lib)
# Install tool executables without suffix to ../tools
install(TARGETS fingerprint.x neighbors.x trnset_info.x trnset2ASCII.x
        RUNTIME DESTINATION ../tools)

# Set output directories and names
# Set output names (RUNTIME_OUTPUT_DIRECTORY inherits global setting: ${CMAKE_BINARY_DIR}/../bin)
set_target_properties(generate.x PROPERTIES
    OUTPUT_NAME generate.x${AENET_SUFFIX}
)
set_target_properties(train.x PROPERTIES
    OUTPUT_NAME train.x${AENET_SUFFIX}
)
set_target_properties(predict.x PROPERTIES
    OUTPUT_NAME predict.x${AENET_SUFFIX}
)
# Set tool output directory explicitly to override global setting
set_target_properties(fingerprint.x neighbors.x trnset_info.x trnset2ASCII.x PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../tools
)
# Set test executable output directory explicitly to override global setting
set_target_properties(test_feedforward test_geometry test_lclist test_sfbasis test_sfsetup test_symmfunc test_xsflib PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../unittests
)
# Output directories for libs already set globally
