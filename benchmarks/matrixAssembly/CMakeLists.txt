# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

add_executable(matrixAssembly "matrixAssembly.cpp")

target_compile_definitions(matrixAssembly PUBLIC WM_LABEL_SIZE=32 NoRepository WM_DP OPENFOAM=2306
                                                 OMPI_SKIP_MPICXX)

target_link_libraries(matrixAssembly PUBLIC OpenFOAM)
