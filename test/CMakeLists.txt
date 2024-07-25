# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.4.0 # or a later release
)
FetchContent_MakeAvailable(Catch2)

enable_testing()

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

function(foam_adapter_unit_test TEST WORKING_DIRECTORY)

  add_executable(${TEST} "test_${TEST}.cpp")

  target_compile_definitions(${TEST} PUBLIC WM_LABEL_SIZE=32 NoRepository WM_DP OPENFOAM=2306
                                            OMPI_SKIP_MPICXX)

  target_include_directories(${TEST} PUBLIC ${CMAKE_SOURCE_DIR}/include)

  target_link_libraries(${TEST} NeoFOAM OpenFOAM FoamAdapter Catch2::Catch2 Kokkos::kokkos)

  add_test(
    NAME adapter_${TEST}
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${TEST}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test/${WORKING_DIRECTORY})

endfunction()

foam_adapter_unit_test(geometricFields)
foam_adapter_unit_test(operators setup_operator)
foam_adapter_unit_test(unstructuredMesh setup_unstructuredMesh)
