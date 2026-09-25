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

# Exercise release preparation only in a disposable source fixture.
file(MAKE_DIRECTORY "${WORK}")
file(COPY "${SOURCE}/prepare-release.sh" "${SOURCE}/license-header.txt"
  DESTINATION "${WORK}")
foreach(name sample.F90 tests/sample.f90 tools/sample.f90
    ext/sample.f90 Makefile Makefile.inc makefiles/Makefile.test)
  file(WRITE "${WORK}/${name}" "fixture\n")
endforeach()
file(WRITE "${WORK}/tests/sample.f90"
  "!+ Copyright (C) 2012-2019 Nongnuch Artrith and Alexander Urban\n"
  "!+ Copyright (C) 2018 John Kitchin\n"
  "!+ Copyright (C) 2025 Shusuke Kasamatsu\nprogram sample\n")
file(WRITE "${WORK}/sample.F90"
  "!+ Copyright (C) 2012-2019 Nongnuch Artrith and Alexander Urban\nprogram sample\n")
file(WRITE "${WORK}/VERSION" "2.0.4\n")
foreach(version v2.0.5 2.0 2.0.5-rc1 02.0.5 "2.0.5\n3.0.0")
  execute_process(COMMAND bash prepare-release.sh "${version}"
    WORKING_DIRECTORY "${WORK}" RESULT_VARIABLE result
    OUTPUT_VARIABLE output ERROR_VARIABLE error)
  file(READ "${WORK}/VERSION" actual)
  if(NOT result STREQUAL "1" OR NOT actual STREQUAL "2.0.4\n")
    message(FATAL_ERROR
      "Invalid version mutated fixture or succeeded: ${version}")
  endif()
endforeach()
execute_process(COMMAND bash prepare-release.sh 2.0.5
  WORKING_DIRECTORY "${WORK}" RESULT_VARIABLE result
  OUTPUT_VARIABLE output ERROR_VARIABLE error)
file(READ "${WORK}/VERSION" actual)
file(READ "${WORK}/tests/sample.f90" developer_header)
file(READ "${WORK}/sample.F90" founder_header)
if(NOT result STREQUAL "0" OR NOT actual STREQUAL "2.0.5\n"
    OR NOT output MATCHES "git tag -a v2.0.5"
    OR NOT developer_header MATCHES
      "Copyright \\(C\\) 2012-2026 Nongnuch Artrith and Alexander Urban"
    OR NOT developer_header MATCHES
      "Copyright \\(C\\) 2018 John Kitchin"
    OR NOT developer_header MATCHES
      "Copyright \\(C\\) 2025 Shusuke Kasamatsu"
    OR NOT founder_header MATCHES
      "Copyright \\(C\\) 2012-2026 Nongnuch Artrith and Alexander Urban")
  message(FATAL_ERROR
    "Release preparation failed: ${result}; ${output}; ${error}")
endif()
