# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

add_executable(explicitOperators "explicitOperators.cpp")

target_compile_definitions(explicitOperators PUBLIC WM_LABEL_SIZE=32 NoRepository WM_DP
                                                    OPENFOAM=2306 OMPI_SKIP_MPICXX)

target_link_libraries(explicitOperators PUBLIC OpenFOAM)

install(TARGETS explicitOperators DESTINATION bin)
