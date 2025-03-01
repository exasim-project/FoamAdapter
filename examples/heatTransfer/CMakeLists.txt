# SPDX-License-Identifier: Unlicense
#
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

add_executable(heatTransfer heatTransfer.cpp)

target_link_libraries(heatTransfer FoamAdapter)
