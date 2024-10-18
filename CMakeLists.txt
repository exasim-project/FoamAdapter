# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2024 NeoFOAM authors

cmake_minimum_required(VERSION 3.22.0)

project(
  FoamAdapter
  LANGUAGES C CXX
  VERSION 0.1.0
  DESCRIPTION "An implementation of FOAM")

if(NOT DEFINED ENV{FOAM_SRC})
  message(FATAL_ERROR "You must source OpenFOAM before building FOAM_ADAPTER")
endif()
if(NOT ENV{FOAM_API} STREQUAL "2406")
  message(WARNING "OpenFOAM version != 2406")
endif()

option(FOAMADAPTER_BUILD_EXAMPLES "Build the NeoFOAM examples" OFF)
option(FOAMADAPTER_BUILD_TESTS "Build the unit tests" OFF)
option(FOAMADAPTER_BUILD_BENCHMARKS "Build benchmarks" OFF)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# Add the cmake folder so the find_package command finds custom packages
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Set the output directories for all binaries and libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

include(cmake/OpenFOAM.cmake)

add_subdirectory(include)
add_subdirectory(NeoFOAM)
add_subdirectory(src)

if(FOAMADAPTER_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(FOAMADAPTER_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()
