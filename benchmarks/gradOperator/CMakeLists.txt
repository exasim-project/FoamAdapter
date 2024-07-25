# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors
add_executable(gradOperator "gradOperator.cpp")

target_compile_definitions(gradOperator PUBLIC WM_LABEL_SIZE=32 NoRepository WM_DP OPENFOAM=2306
                                               OMPI_SKIP_MPICXX)

target_link_libraries(gradOperator PUBLIC OpenFOAM NeoFOAM FoamAdapter Kokkos::kokkos)

install(TARGETS gradOperator DESTINATION bin)
