# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

include(${PROJECT_SOURCE_DIR}/cmake/Catch2.cmake)

add_executable(test_geometricFields "test_geometricFields.cpp")
set_target_properties(test_geometricFields PROPERTIES CXX_STANDARD 20)

target_compile_definitions(test_geometricFields PUBLIC WM_LABEL_SIZE=32 NoRepository WM_DP OPENFOAM=2306 OMPI_SKIP_MPICXX)

target_include_directories(
test_geometricFields
 SYSTEM
  PUBLIC $ENV{FOAM_SRC}/finiteVolume/lnInclude
         $ENV{FOAM_SRC}/meshTools/lnInclude
         $ENV{FOAM_SRC}/meshTools/lnInclude
         $ENV{FOAM_SRC}/OpenFOAM/lnInclude
         $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
         ${CMAKE_SOURCE_DIR}/include
         ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(
  test_geometricFields
  PUBLIC $ENV{FOAM_LIBBIN}/libOpenFOAM.so
         $ENV{FOAM_LIBBIN}/libfiniteVolume.so
         $ENV{FOAM_LIBBIN}/libmeshTools.so
         $ENV{FOAM_LIBBIN}/$ENV{FOAM_MPI}/libPstream.so
         NeoFOAM
         FoamAdapter
         Catch2::Catch2
         Kokkos::kokkos)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

include(CTest)
include(Catch)
enable_testing()
add_test(NAME geometricFields_test COMMAND test_geometricFields)

add_executable(test_operators "test_operators.cpp")
set_target_properties(test_operators PROPERTIES CXX_STANDARD 20)

target_compile_definitions(test_operators PUBLIC WM_LABEL_SIZE=32 NoRepository WM_DP OPENFOAM=2306 OMPI_SKIP_MPICXX)

target_include_directories(
test_operators
 SYSTEM
  PUBLIC $ENV{FOAM_SRC}/finiteVolume/lnInclude
         $ENV{FOAM_SRC}/meshTools/lnInclude
         $ENV{FOAM_SRC}/meshTools/lnInclude
         $ENV{FOAM_SRC}/OpenFOAM/lnInclude
         $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
         ${CMAKE_SOURCE_DIR}/include
         ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(
  test_operators
  PUBLIC $ENV{FOAM_LIBBIN}/libOpenFOAM.so
         $ENV{FOAM_LIBBIN}/libfiniteVolume.so
         $ENV{FOAM_LIBBIN}/libmeshTools.so
         $ENV{FOAM_LIBBIN}/$ENV{FOAM_MPI}/libPstream.so
         NeoFOAM
         FoamAdapter
         Catch2::Catch2
         Kokkos::kokkos)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

include(CTest)
include(Catch)
enable_testing()
add_test(NAME operator_test COMMAND test_operators)

add_executable(test_unstructuredMesh "test_unstructuredMesh.cpp")
set_target_properties(test_unstructuredMesh PROPERTIES CXX_STANDARD 20)

target_compile_definitions(test_unstructuredMesh PUBLIC WM_LABEL_SIZE=32 NoRepository WM_DP OPENFOAM=2306 OMPI_SKIP_MPICXX)

target_include_directories(
test_unstructuredMesh
 SYSTEM
  PUBLIC $ENV{FOAM_SRC}/finiteVolume/lnInclude
         $ENV{FOAM_SRC}/meshTools/lnInclude
         $ENV{FOAM_SRC}/meshTools/lnInclude
         $ENV{FOAM_SRC}/OpenFOAM/lnInclude
         $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
         ${CMAKE_SOURCE_DIR}/include
         ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(
  test_unstructuredMesh
  PUBLIC $ENV{FOAM_LIBBIN}/libOpenFOAM.so
         $ENV{FOAM_LIBBIN}/libfiniteVolume.so
         $ENV{FOAM_LIBBIN}/libmeshTools.so
         $ENV{FOAM_LIBBIN}/$ENV{FOAM_MPI}/libPstream.so
         NeoFOAM
         FoamAdapter
         Catch2::Catch2
         Kokkos::kokkos)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

include(CTest)
include(Catch)
enable_testing()
add_test(NAME unstructuredMesh_test COMMAND test_unstructuredMesh)
