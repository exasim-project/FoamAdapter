# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 FoamAdapter authors

add_subdirectory(catch2)

function(foam_adapter_unit_test TEST SETUP_DIRECTORY)
  if(NOT DEFINED "adapter_WORKING_DIRECTORY")
    set(adapter_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests)
  endif()

  add_executable(adapter_${TEST} "test_${TEST}.cpp")

  target_compile_definitions(adapter_${TEST} PUBLIC OMPI_SKIP_MPICXX)

  target_link_libraries(adapter_${TEST} foamadapter_catch_main OpenFOAM NeoN)

  set_target_properties(adapter_${TEST} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                   ${adapter_WORKING_DIRECTORY})

  if(NeoN_WITH_PETSC)
    message("Hallo")
    message("${PETSc_LIBRARIES}")
    target_link_libraries(adapter_${TEST} ${PETSc_LIBRARIES} fftw3 mpi)
    target_include_directories(adapter_${TEST} PUBLIC ${PETSc_INCLUDE_DIRS})
    target_link_directories(adapter_${TEST} PUBLIC ${PETSc_LIBRARY_DIRS})
    link_directories(${PETSc_LIBRARY_DIRS})
  endif()

  add_test(
    NAME adapter_${TEST}
    COMMAND ${adapter_WORKING_DIRECTORY}/adapter_${TEST}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test/${SETUP_DIRECTORY})

endfunction()

foam_adapter_unit_test(geometricFields setup_operator)
foam_adapter_unit_test(operators setup_operator)
foam_adapter_unit_test(stencils setup_stencil3D)
foam_adapter_unit_test(sparsityPattern setup_stencil3D)
foam_adapter_unit_test(implicitOperators setup_stencil3D)
foam_adapter_unit_test(readDict setup_operator)
foam_adapter_unit_test(pressureVelocityCoupling setup_pressureVelocityCoupling)
foam_adapter_unit_test(unstructuredMesh setup_unstructuredMesh)
foam_adapter_unit_test(advection setup_advection)
