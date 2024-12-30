# SPDX-License-Identifier: Unlicense
#
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

add_executable(scalarAdvection scalarAdvection.cpp)

target_link_libraries(scalarAdvection FoamAdapter)

add_executable(scalarFoamAdvection scalarFoamAdvection.cpp)

target_link_libraries(scalarFoamAdvection FoamAdapter)
