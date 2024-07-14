# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors
add_executable(test_treeAMR "test_treeAMR.cpp")

find_package(MPI REQUIRED)
find_package(ZLIB REQUIRED)

set_target_properties(test_treeAMR PROPERTIES CXX_STANDARD 20)

# target_compile_definitions(test_treeAMR PUBLIC WM_LABEL_SIZE=32 NoRepository WM_DP OPENFOAM=2306
# OMPI_SKIP_MPICXX)

target_include_directories(
  test_treeAMR SYSTEM
  PUBLIC ${CMAKE_SOURCE_DIR}/FoamAdapter/include ${CMAKE_SOURCE_DIR}/NeoFOAM/include
         ${MPI_INCLUDE_PATH} ${ZLIB_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include
         ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(
  test_treeAMR
  PUBLIC ${CMAKE_SOURCE_DIR}/FoamAdapter/lib/libt8.a
         ${CMAKE_SOURCE_DIR}/FoamAdapter/lib/libsc.a
         ${CMAKE_SOURCE_DIR}/FoamAdapter/lib/libp4est.a
         MPI::MPI_CXX
         NeoFOAM
         Kokkos::kokkos
         ${ZLIB_LIBRARIES})

install(TARGETS test_treeAMR DESTINATION bin)
