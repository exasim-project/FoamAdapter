# SPDX-License-Identifier: Unlicense
#
# SPDX-FileCopyrightText: 2023 FoamAdapter authors

add_library(FoamAdapter SHARED)

target_compile_definitions(FoamAdapter INTERFACE namespaceFoam=1)
target_link_libraries(FoamAdapter PUBLIC FoamAdapter_public_api OpenFOAM)

target_sources(
  FoamAdapter
  PRIVATE "algorithms/pressureVelocityCoupling.cpp"
          "auxiliary/convert.cpp"
          "auxiliary/foamDictionary.cpp"
          "auxiliary/setup.cpp"
          "auxiliary/writers.cpp"
          "auxiliary/comparison.cpp"
          # "datastructures/foamMesh.cpp"
          "datastructures/meshAdapter.cpp"
          "compatibility/fvSolution.cpp")

install(TARGETS FoamAdapter)
