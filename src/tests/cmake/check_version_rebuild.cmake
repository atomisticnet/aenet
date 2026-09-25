# This file is part of the AENET package.  Copyright (C) 2012-2026
# Nongnuch Artrith and Alexander Urban  This Source Code Form is subject
# to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at
# http://mozilla.org/MPL/2.0/.  This program is distributed in the hope
# that it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the Mozilla Public License, v. 2.0, for more details.
# ---------------------------------------------------------------------
# If you make use of AENET for your publication, please cite: [1] N.
# Artrith and A. Urban, Comput. Mater. Sci. 114 (2016) 135-150. [2] J.
# Behler and M. Parrinello, Phys. Rev. Lett. 98 (2007) 146401.  If you
# used the Chebyshev descriptor, please cite: [3] N. Artrith, A. Urban,
# and G. Ceder, PRB 96 (2017) 014112.

# Separate integration check: configure/build a disposable source copy twice.
# Required: SOURCE, WORK, COMPILER, BLA_VENDOR, and EXECUTABLE_SUFFIX.
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}/source/lib")
file(COPY "${SOURCE}/CMakeLists.txt" DESTINATION "${WORK}/source")
file(COPY "${SOURCE}/src" DESTINATION "${WORK}/source"
  FILES_MATCHING PATTERN "*.f90" PATTERN "*.F90" PATTERN "*.cmake"
  PATTERN "CMakeLists.txt" PATTERN "VERSION")
file(COPY "${SOURCE}/lib/Lbfgsb.3.0.tar.gz"
  DESTINATION "${WORK}/source/lib")
execute_process(COMMAND "${CMAKE_COMMAND}" -S "${WORK}/source"
  -B "${WORK}/build" -DBUILD_AENET=ON "-DCMAKE_Fortran_COMPILER=${COMPILER}"
  "-DBLA_VENDOR=${BLA_VENDOR}" "-DUSE_OPENBLAS=${USE_OPENBLAS}"
  "-DCMAKE_PREFIX_PATH=${PREFIX_PATH}"
  "-DCMAKE_EXE_LINKER_FLAGS=${EXE_LINKER_FLAGS}"
  "-DCMAKE_SHARED_LINKER_FLAGS=${SHARED_LINKER_FLAGS}"
  -DCMAKE_BUILD_TYPE=Release
  RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "Initial configure failed: ${result}")
endif()
file(READ "${WORK}/source/src/VERSION" initial)
string(STRIP "${initial}" initial)
foreach(version "${initial}" "7.8.9")
  # No explicit reconfigure: the build must notice the changed source VERSION.
  file(WRITE "${WORK}/source/src/VERSION" "${version}\n")
  execute_process(COMMAND "${CMAKE_COMMAND}" --build "${WORK}/build"
    --target main lib --parallel 1 RESULT_VARIABLE result)
  if(NOT result STREQUAL "0")
    message(FATAL_ERROR "Build failed for ${version}: ${result}")
  endif()
  foreach(name generate train predict)
    set(PROGRAM "${WORK}/build/bin/${name}.x${EXECUTABLE_SUFFIX}")
    set(NAME "${name}.x")
    set(VERSION "${version}")
    set(MODE version)
    set(saved_work "${WORK}")
    set(WORK "${WORK}/run-${name}-${version}")
    include("${SOURCE}/src/tests/cmake/check_startup.cmake")
    set(WORK "${saved_work}")
  endforeach()
endforeach()
if(CMAKE_HOST_APPLE)
  execute_process(COMMAND otool -L "${WORK}/build/lib/libaenet.dylib"
    OUTPUT_VARIABLE metadata RESULT_VARIABLE result)
  if(NOT result STREQUAL "0" OR NOT metadata MATCHES "libaenet.7.dylib"
      OR NOT metadata MATCHES "current version 7.8.9")
    message(FATAL_ERROR "Incorrect Mach-O version metadata: ${metadata}")
  endif()
elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
  execute_process(COMMAND readelf -d "${WORK}/build/lib/libaenet.so"
    OUTPUT_VARIABLE metadata RESULT_VARIABLE result)
  if(NOT result STREQUAL "0" OR NOT metadata MATCHES "SONAME.*libaenet.so.7")
    message(FATAL_ERROR "Incorrect ELF SONAME: ${metadata}")
  endif()
  if(NOT EXISTS "${WORK}/build/lib/libaenet.so.7.8.9")
    message(FATAL_ERROR "Missing full-version ELF library")
  endif()
endif()
