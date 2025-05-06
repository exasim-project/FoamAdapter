# SPDX-License-Identifier: Unlicense
#
# SPDX-FileCopyrightText: 2023 FoamAdapter authors

add_executable(heatTransfer heatTransfer.cpp)

target_link_libraries(heatTransfer FoamAdapter)

if(NeoN_WITH_PETSC)
  target_link_libraries(heatTransfer ${PETSc_LIBRARIES} fftw3 mpi)
  target_include_directories(heatTransfer PUBLIC ${PETSc_INCLUDE_DIRS})
  target_link_directories(heatTransfer PUBLIC ${PETSc_LIBRARY_DIRS})
  link_directories(${PETSc_LIBRARY_DIRS})
endif()
