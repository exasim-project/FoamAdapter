# SPDX-License-Identifier: Unlicense
#
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

add_executable(icoNeoFoam icoNeoFoam.cpp)

target_link_libraries(icoNeoFoam FoamAdapter)
