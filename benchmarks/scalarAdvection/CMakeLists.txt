# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors
add_executable(scalarAdvection "scalarAdvection.cpp")

set_target_properties(scalarAdvection PROPERTIES CXX_STANDARD 20)

target_compile_definitions(scalarAdvection PUBLIC WM_LABEL_SIZE=32 NoRepository WM_DP OPENFOAM=2306
                                                  OMPI_SKIP_MPICXX)

target_link_libraries(scalarAdvection PUBLIC OpenFOAM NeoFOAM FoamAdapter Kokkos::kokkos)

install(TARGETS scalarAdvection DESTINATION bin)
