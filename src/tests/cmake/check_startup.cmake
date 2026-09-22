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

# Run in an isolated directory; MPI launchers may write diagnostics to stderr.
cmake_policy(SET CMP0054 NEW)

file(MAKE_DIRECTORY "${WORK}")
if(MODE STREQUAL "version")
  set(args --version)
elseif(MODE STREQUAL "missing")
  set(args nonexistent.in)
elseif(MODE STREQUAL "noargs")
  set(args)
else()
  message(FATAL_ERROR "Unknown startup test mode: ${MODE}")
endif()
execute_process(COMMAND ${LAUNCHER} "${PROGRAM}" ${POSTFLAGS} ${args}
  WORKING_DIRECTORY "${WORK}" TIMEOUT 20
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(MODE STREQUAL "version")
  if(NOT result STREQUAL "0" OR NOT output STREQUAL "${NAME} ${VERSION}\n")
    message(FATAL_ERROR
      "Version query: result=${result}; stdout=${output}; stderr=${error}")
  endif()
else()
  if(NOT result MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR
      "Input error must exit nonzero without hanging: ${result}; ${error}")
  endif()
  if(NOT error MATCHES "Error: (No input file|File not found)")
    message(FATAL_ERROR "Expected input diagnostic; got: ${error}")
  endif()
endif()
file(GLOB files "${WORK}/*" "${WORK}/.*")
if(files)
  message(FATAL_ERROR "Startup created unexpected files: ${files}")
endif()
