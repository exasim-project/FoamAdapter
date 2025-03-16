# SPDX-License-Identifier: Unlicense
#
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

add_executable(icoNeoFoam icoNeoFoam.cpp)

target_link_libraries(icoNeoFoam FoamAdapter)

if(NEOFOAM_WITH_PETSC)
  message("Hallo")
  message("${PETSc_LIBRARIES}")
  target_link_libraries(icoNeoFoam ${PETSc_LIBRARIES} fftw3 mpi)
  target_include_directories(icoNeoFoam PUBLIC ${PETSc_INCLUDE_DIRS})
  target_link_directories(icoNeoFoam PUBLIC ${PETSc_LIBRARY_DIRS})
  link_directories(${PETSc_LIBRARY_DIRS})
endif()
