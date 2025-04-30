# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 FoamAdapter authors

include(cmake/CPM.cmake)

if(FoamAdapter_BUILD_TESTS OR FoamAdapter_BUILD_BENCHMARKS)
  cpmaddpackage(NAME Catch2 GITHUB_REPOSITORY catchorg/Catch2 VERSION 3.4.0)
endif()

if(NEOFOAM_NEON_VIA_CPM)
  # cpmaddpackage(NAME NeoN GITHUB_REPOSITORY exasim-project/NeoN SYSTEM YES VERSION 0.1)
  cpmaddpackage("gh:exasim-project/NeoN#main")
endif()
