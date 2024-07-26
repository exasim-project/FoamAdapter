# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.4.0 # or a later release
)
FetchContent_MakeAvailable(Catch2)

include(CTest)
include(Catch)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

function(foam_adapter_unit_test TEST SETUP_DIRECTORY)
  if(NOT DEFINED "adapter_WORKING_DIRECTORY")
    set(adapter_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests)
  endif()

  add_executable(adapter_${TEST} "test_${TEST}.cpp")

  target_compile_definitions(adapter_${TEST} PUBLIC WM_LABEL_SIZE=32 NoRepository WM_DP
                                                    OPENFOAM=2306 OMPI_SKIP_MPICXX)

  target_include_directories(adapter_${TEST} PUBLIC ${CMAKE_SOURCE_DIR}/include)

  target_link_libraries(adapter_${TEST} NeoFOAM OpenFOAM FoamAdapter Catch2::Catch2WithMain)

  set_target_properties(adapter_${TEST} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                   ${adapter_WORKING_DIRECTORY})

  add_test(
    NAME adapter_${TEST}
    COMMAND ${adapter_WORKING_DIRECTORY}/adapter_${TEST}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test/${SETUP_DIRECTORY})

endfunction()

foam_adapter_unit_test(geometricFields setup_operator)
foam_adapter_unit_test(operators setup_operator)
foam_adapter_unit_test(unstructuredMesh setup_unstructuredMesh)
