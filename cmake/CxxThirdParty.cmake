# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 FoamAdapter authors

include(cmake/CPM.cmake)

if(FoamAdapter_BUILD_TESTS OR FoamAdapter_BUILD_BENCHMARKS)
  cpmaddpackage(NAME Catch2 GITHUB_REPOSITORY catchorg/Catch2 VERSION 3.4.0)
endif()

if(FOAMADAPTER_NEON_VIA_CPM)
  if(NOT DEFINED FOAMADAPTER_NEON_VERSION)
    set(FOAMADAPTER_NEON_VERSION
        "e2bd3c96fdfcde9156eae5a474d8c8b8c61bbaf1"
        CACHE INTERNAL "")
  endif()

  cpmaddpackage(
    NAME
    NeoN
    GITHUB_REPOSITORY
    exasim-project/NeoN
    GIT_TAG
    ${FOAMADAPTER_NEON_VERSION}
    SYSTEM
    YES)
endif()
