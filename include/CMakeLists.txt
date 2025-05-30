# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 FoamAdapter authors

add_library(FoamAdapter_public_api INTERFACE) # dummy target

set_property(
  TARGET FoamAdapter_public_api
  APPEND
  PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}")

set_property(
  TARGET FoamAdapter_public_api
  APPEND
  PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_BINARY_DIR}")

target_link_libraries(FoamAdapter_public_api INTERFACE NeoN::NeoN NeoN::NeoN_public_api)

# Get list of some *.hpp files in folder include
file(GLOB_RECURSE include_files *.hpp)

# Convert the list of files into #includes
foreach(include_file ${include_files})
  set(include_statements "${include_statements}#include \"${include_file}\"\n")
endforeach()

configure_file(${CMAKE_CURRENT_LIST_DIR}/FoamAdapter/FoamAdapter.hpp.in
               ${CMAKE_CURRENT_BINARY_DIR}/FoamAdapter/FoamAdapter.hpp)
