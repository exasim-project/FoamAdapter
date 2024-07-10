# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

cmake_minimum_required(VERSION 3.22.0)

project(
  FoamAdapter
  LANGUAGES C CXX
  VERSION 0.1.0
  DESCRIPTION "An implementation of FOAM")

add_subdirectory(NeoFOAM)

include(cmake/PreventInSourceBuilds.cmake)
include(cmake/StandardProjectSettings.cmake)

option(NEOFOAM_DP_SCALAR "double precision scalar" TRUE)
option(NEOFOAM_DP_LABELS "double precision label" TRUE)
option(NEOFOAM_BUILD_APPS "Build the applications folder" TRUE)
option(NEOFOAM_BUILD_TESTS "Build the unit tests" TRUE)

option(NEOFOAM_ENABLE_SANITIZE_ADDRESS "Enable address sanitizer" FALSE)
option(NEOFOAM_ENABLE_SANITIZE_LEAK "Enable leak sanitizer" FALSE)
option(NEOFOAM_ENABLE_SANITIZE_UB "Enable undefined behaviour sanitizer" FALSE)
option(NEOFOAM_ENABLE_SANITIZE_THREAD "Enable thread sanitizer" FALSE)
option(NEOFOAM_ENABLE_SANITIZE_MEMORY "Enable memory sanitizer" FALSE)
option(NEOFOAM_ENABLE_CPP_CHECK "Enable cpp check static analyzer" FALSE)
option(NEOFOAM_ENABLE_CLANG_TIDY "Enable clang tidy static analyzer" FALSE)
option(NEOFOAM_ENABLE_PCH "Enable precompiled header" FALSE)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 20)

# add_library(project_options INTERFACE) include(cmake/CompilerWarnings.cmake)
# add_library(project_warnings INTERFACE) myproject_set_project_warnings(project_warnings
# ${myproject_WARNINGS_AS_ERRORS} FALSE "" "" "")

# include(cmake/Sanitizer.cmake) enable_sanitizers( project_options NEOFOAM_ENABLE_SANITIZE_ADDRESS
# NEOFOAM_ENABLE_SANITIZE_LEAK NEOFOAM_ENABLE_SANITIZE_UB NEOFOAM_ENABLE_SANITIZE_THREAD
# NEOFOAM_ENABLE_SANITIZE_MEMORY)

if(NEOFOAM_ENABLE_PCH)
  target_precompile_headers(project_options INTERFACE <vector> <string> <utility>)
endif()

include(cmake/StaticAnalyzers.cmake)
if(NEOFOAM_ENABLE_CPP_CHECK)
  myproject_enable_cppcheck(FALSE "")
endif()
if(NEOFOAM_ENABLE_CLANG_TIDY)
  myproject_enable_cppcheck(FALSE "")
endif()

add_subdirectory(src)
add_subdirectory(applications)

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

if(NEOFOAM_ENABLE_IWYU)
  find_program(iwyu_path NAMES include-what-you-use iwyu REQUIRED)
  set(IWYU_PATH_AND_OPTIONS ${iwyu_path})
  set_property(TARGET NeoFOAM PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${IWYU_PATH_AND_OPTIONS})
endif()
