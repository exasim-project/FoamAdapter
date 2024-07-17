# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

add_executable(explicitOperators "explicitOperators.cpp")

target_compile_definitions(explicitOperators PUBLIC WM_LABEL_SIZE=32 NoRepository WM_DP
                                                    OPENFOAM=2306 OMPI_SKIP_MPICXX)

target_include_directories(
  explicitOperators SYSTEM
  PUBLIC $ENV{FOAM_SRC}/finiteVolume/lnInclude $ENV{FOAM_SRC}/meshTools/lnInclude
         $ENV{FOAM_SRC}/OpenFOAM/lnInclude $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
         ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(
  explicitOperators
  PUBLIC $ENV{FOAM_LIBBIN}/libOpenFOAM.so $ENV{FOAM_LIBBIN}/libfiniteVolume.so
         $ENV{FOAM_LIBBIN}/libmeshTools.so $ENV{FOAM_LIBBIN}/$ENV{FOAM_MPI}/libPstream.so
         $ENV{FOAM_LIBBIN}/libOSspecific.o)

install(TARGETS explicitOperators DESTINATION bin)
