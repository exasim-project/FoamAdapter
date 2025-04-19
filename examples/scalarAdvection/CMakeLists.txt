# SPDX-License-Identifier: Unlicense
#
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

add_executable(scalarAdvection scalarAdvection.cpp)

target_link_libraries(scalarAdvection FoamAdapter)
target_compile_definitions(scalarAdvection INTERFACE namespaceFoam=1)

if(NeoN_WITH_PETSC)
  target_link_libraries(scalarAdvection ${PETSc_LIBRARIES} fftw3 mpi)
  target_include_directories(scalarAdvection PUBLIC ${PETSc_INCLUDE_DIRS})
  target_link_directories(scalarAdvection PUBLIC ${PETSc_LIBRARY_DIRS})
  link_directories(${PETSc_LIBRARY_DIRS})
endif()
