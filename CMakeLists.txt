# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

cmake_minimum_required(VERSION 3.22.0)

project(
  FoamAdapter
  LANGUAGES C CXX
  VERSION 0.1.0
  DESCRIPTION "An implementation of FOAM")

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)

# Add the cmake folder so the find_package command finds custom packages
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Set the output directories for all binaries and libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Eventhough we compile NeoFOAM as SHARED library we need to explicitly enable PIC for all targets
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(cmake/FoamAdapaterConfig.cmake)

add_subdirectory(NeoFOAM)
add_subdirectory(src)

