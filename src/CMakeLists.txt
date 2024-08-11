# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

add_library(FoamAdapter SHARED)

target_compile_definitions(
  FoamAdapter PUBLIC WM_LABEL_SIZE=$ENV{WM_LABEL_SIZE} NoRepository WM_$ENV{WM_PRECISION_OPTION}
                     OPENFOAM=2406 OMPI_SKIP_MPICXX)

target_include_directories(
  FoamAdapter PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
                     "$<INSTALL_INTERFACE:${CMAKE_INSTALL_FULL_INCLUDEDIR}>")

target_link_libraries(FoamAdapter PUBLIC OpenFOAM NeoFOAM Kokkos::kokkos)

target_sources(
  FoamAdapter PRIVATE "conversion/convert.cpp" "setup/setup.cpp" "readers/foamMesh.cpp"
                      "readers/foamDictionary.cpp" "fvcc/mesh/fvccNeoMesh.cpp")
