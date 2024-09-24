# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

function(foam_adapter_benchmark BENCHMARK)
  if(NOT DEFINED "adapter_WORKING_DIRECTORY")
    set(adapter_benchmark_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/benchmark)
  endif()

  add_executable(adapter_benchmark_${BENCHMARK} "${BENCHMARK}.cpp")

  target_link_libraries(adapter_benchmark_${BENCHMARK} FoamAdapter)

  set_target_properties(adapter_benchmark_${BENCHMARK}
                        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${adapter_benchmark_WORKING_DIRECTORY})

endfunction()

add_subdirectory(explicitOperators)
add_subdirectory(matrixAssembly)
add_subdirectory(gradOperator)
add_subdirectory(scalarAdvection)
