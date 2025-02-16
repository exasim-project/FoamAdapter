# SPDX-License-Identifier: Unlicense
#
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

add_executable(scalarAdvection scalarAdvection.cpp)

target_link_libraries(scalarAdvection FoamAdapter)
target_compile_definitions(scalarAdvection INTERFACE namespaceFoam=1)
