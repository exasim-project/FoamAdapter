# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

add_library(FoamAdapter)

# find_package(Kokkos REQUIRED)
# find_package(MPI REQUIRED)
# find_package(ZLIB REQUIRED)
set_target_properties(FoamAdapter PROPERTIES CXX_STANDARD 17) # does not compile with c++20
target_compile_definitions(FoamAdapter PUBLIC WM_LABEL_SIZE=32 NoRepository WM_DP OPENFOAM=2306 OMPI_SKIP_MPICXX)

target_include_directories(
  FoamAdapter PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
                      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_FULL_INCLUDEDIR}>"
                      $ENV{FOAM_SRC}/OpenFOAM/lnInclude
                      $ENV{FOAM_SRC}/finiteVolume/lnInclude
                      $ENV{FOAM_SRC}/meshTools/lnInclude
                      $ENV{FOAM_SRC}/meshTools/lnInclude
                      $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
                      ${CMAKE_SOURCE_DIR}/NeoFOAM/include)
                    #   ${CMAKE_SOURCE_DIR}/include)
                    #   ${MPI_INCLUDE_PATH}
                    #   ${ZLIB_INCLUDE_DIRS})

target_link_libraries(
    FoamAdapter
    PUBLIC $ENV{FOAM_LIBBIN}/libOpenFOAM.so
           $ENV{FOAM_LIBBIN}/libfiniteVolume.so
           $ENV{FOAM_LIBBIN}/libmeshTools.so
           $ENV{FOAM_LIBBIN}/$ENV{FOAM_MPI}/libPstream.so)
            # ${CMAKE_SOURCE_DIR}/FoamAdapter/lib/libt8.a
            # ${CMAKE_SOURCE_DIR}/FoamAdapter/lib/libsc.a
            # ${CMAKE_SOURCE_DIR}/FoamAdapter/lib/libp4est.a
            # MPI::MPI_CXX
            # ${ZLIB_LIBRARIES})

target_sources(FoamAdapter PRIVATE "setup/setup.cpp")

add_subdirectory(benchmarks)
add_subdirectory(test)



