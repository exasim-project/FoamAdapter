# SPDX-License-Identifier: Unlicense
#
# SPDX-FileCopyrightText: 2023 FoamAdapter authors

add_executable(neoIcoFoam neoIcoFoam.cpp)

target_link_libraries(neoIcoFoam FoamAdapter)
