# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

add_library(FoamAdapter_public_api INTERFACE) # dummy target

set_property(
  TARGET FoamAdapter_public_api
  APPEND
  PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}")

set_property(
  TARGET FoamAdapter_public_api
  APPEND
  PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_BINARY_DIR}")
