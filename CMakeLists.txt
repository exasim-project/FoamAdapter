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

option(NEOFOAM_DEVEL_TOOLS "Add development tools to the build system" FALSE)

# Add the cmake folder so the find_package command finds custom packages
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Set the output directories for all binaries and libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Eventhough we compile NeoFOAM as SHARED library we need to explicitly enable PIC for all targets
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(NeoFOAM)
add_subdirectory(src)

if(NEOFOAM_DEVEL_TOOLS)
  find_program(PRE_COMMIT pre-commit)
  if(NOT PRE_COMMIT)
    message(
      FATAL_ERROR
        "The pre-commit command was not found. "
        "It is necessary if you want to commit changes to NeoFOAM. "
        "If that is not the case, set NEOFOAM_DEVEL_TOOLS=OFF. "
        "Otherwise install pre-commit via pipx (or pip) using:\n"
        "    pipx install pre-commit")
  endif()

  execute_process(
    COMMAND "${PRE_COMMIT}" "install"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    RESULT_VARIABLE pre-commit-result
    OUTPUT_VARIABLE pre-commit-output
    ERROR_VARIABLE pre-commit-error)
  if(pre-commit-result)
    message(
      FATAL_ERROR "Failed to install the git hooks via pre-commit."
                  "Please check the error message:\n" "${pre-commit-output}\n${pre-commit-error}")
  endif()
  if(pre-commit-output MATCHES "^Running in migration mode with existing hooks")
    message(WARNING "An existing git hook was encountered during `pre-commit install`."
                    "The old git hook will also be executed."
                    "Consider removing it with `pre-commit install -f`")
  elseif(NOT pre-commit-output MATCHES "^pre-commit installed at")
    message(WARNING "`pre-commit install` did not exit normally."
                    "Please check the output message:\n" "${pre-commit-output}")
  endif()

  add_custom_target(
    pre-commit
    COMMENT "execute precommit"
    COMMAND bash -c "${PRE_COMMIT} run"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    VERBATIM)
endif()
