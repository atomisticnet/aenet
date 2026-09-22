# This file is part of the AENET package.  Copyright (C) 2012-2019
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

# Exercise the installed interface from a clean prefix. The build-tree
# executable suffix is deliberately passed in so the test can reject it.
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
set(prefix "${WORK}/prefix")

execute_process(COMMAND "${CMAKE_COMMAND}" --build "${BUILD}"
  --config Release --target build_all --parallel 1
  RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "Artifact build failed: ${result}")
endif()

execute_process(COMMAND "${CMAKE_COMMAND}" --install "${BUILD}"
  --config Release --prefix "${prefix}" RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "Installation failed: ${result}")
endif()

set(programs generate.x train.x predict.x)
foreach(program IN LISTS programs)
  if(NOT EXISTS "${prefix}/bin/${program}")
    message(FATAL_ERROR "Missing installed program: bin/${program}")
  endif()
  if(EXECUTABLE_SUFFIX AND
      EXISTS "${prefix}/bin/${program}${EXECUTABLE_SUFFIX}")
    message(FATAL_ERROR
      "Build-tree suffix leaked into installed program: ${program}${EXECUTABLE_SUFFIX}")
  endif()
  execute_process(COMMAND "${prefix}/bin/${program}" --version
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(NOT result STREQUAL "0" OR
      NOT output STREQUAL "${program} ${VERSION}\n")
    message(FATAL_ERROR
      "Installed version query failed: ${program}; ${result}; ${output}; ${error}")
  endif()
endforeach()
file(GLOB installed_programs RELATIVE "${prefix}/bin" "${prefix}/bin/*")
list(SORT installed_programs)
set(expected_programs generate.x predict.x train.x)
if(NOT installed_programs STREQUAL expected_programs)
  message(FATAL_ERROR "Unexpected installed programs: ${installed_programs}")
endif()

foreach(tool fingerprint.x neighbors.x trnset_info.x trnset2ASCII.x)
  if(NOT EXISTS "${prefix}/tools/${tool}")
    message(FATAL_ERROR "Missing installed tool: tools/${tool}")
  endif()
endforeach()
file(GLOB installed_tools RELATIVE "${prefix}/tools" "${prefix}/tools/*")
list(SORT installed_tools)
set(expected_tools fingerprint.x neighbors.x trnset2ASCII.x trnset_info.x)
if(NOT installed_tools STREQUAL expected_tools)
  message(FATAL_ERROR "Unexpected installed tools: ${installed_tools}")
endif()

if(NOT EXISTS "${prefix}/include/aenet.h")
  message(FATAL_ERROR "Missing installed header: include/aenet.h")
endif()
if(NOT EXISTS "${prefix}/lib/${STATIC_PREFIX}aenet${STATIC_SUFFIX}")
  message(FATAL_ERROR "Missing installed static library")
endif()
if(NOT EXISTS "${prefix}/lib/${SHARED_PREFIX}aenet${SHARED_SUFFIX}")
  message(FATAL_ERROR "Missing installed shared-library link")
endif()

set(consumer "${WORK}/consumer")
file(MAKE_DIRECTORY "${consumer}")
file(WRITE "${consumer}/main.c" [=[
#include <aenet.h>

int main(void) {
  int status = -1;
  char name[] = "Cu";
  char *types[] = {name};
  aenet_init(1, types, &status);
  if (status != AENET_OK) return 1;
  aenet_final(&status);
  return status == AENET_OK ? 0 : 2;
}
]=])
file(WRITE "${consumer}/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.15)
project(aenet_installed_consumer LANGUAGES C)
find_path(AENET_INCLUDE_DIR aenet.h REQUIRED
  PATHS "${CMAKE_PREFIX_PATH}/include" NO_DEFAULT_PATH)
find_library(AENET_LIBRARY aenet REQUIRED
  PATHS "${CMAKE_PREFIX_PATH}/lib" NO_DEFAULT_PATH)
add_executable(aenet_consumer main.c)
target_include_directories(aenet_consumer PRIVATE "${AENET_INCLUDE_DIR}")
target_link_libraries(aenet_consumer PRIVATE "${AENET_LIBRARY}")
set_target_properties(aenet_consumer PROPERTIES
  BUILD_RPATH "${CMAKE_PREFIX_PATH}/lib")
]=])
execute_process(COMMAND "${CMAKE_COMMAND}" -S "${consumer}"
  -B "${consumer}/build" "-DCMAKE_PREFIX_PATH=${prefix}"
  RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "Installed C consumer configure failed: ${result}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" --build "${consumer}/build"
  RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "Installed C consumer build failed: ${result}")
endif()
execute_process(COMMAND "${consumer}/build/aenet_consumer"
  RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "Installed C consumer failed: ${result}")
endif()
