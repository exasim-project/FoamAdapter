# SPDX-License-Identifier: Unlicense
#
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

add_library(FoamAdapter SHARED)

target_compile_definitions(FoamAdapter INTERFACE namespaceFoam=1)
target_link_libraries(FoamAdapter PUBLIC FoamAdapter_public_api OpenFOAM NeoN)

target_sources(
  FoamAdapter
  PRIVATE "conversion/convert.cpp" "setup.cpp" "meshAdapter.cpp" "readers/foamDictionary.cpp"
          "finiteVolume/cellCentred/pressureVelocityCoupling/pressureVelocityCoupling.cpp")
