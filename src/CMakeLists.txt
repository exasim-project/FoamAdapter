# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

add_subdirectory(benchmarks)
# add_library(NeoFOAM_GPL)

# # find_package(Kokkos REQUIRED)
# find_package(MPI REQUIRED)
# find_package(ZLIB REQUIRED)

# target_include_directories(
#   NeoFOAM_GPL PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
#                       "$<INSTALL_INTERFACE:${CMAKE_INSTALL_FULL_INCLUDEDIR}>"
#                       ${CMAKE_SOURCE_DIR}/NeoFOAM_GPL/include
#                       ${MPI_INCLUDE_PATH}
#                       ${ZLIB_INCLUDE_DIRS})

# # target_include_directories(
# # test_treeAMR
# #     SYSTEM
# #     PUBLIC ${CMAKE_SOURCE_DIR}/NeoFOAM_GPL/include
# #             ${MPI_INCLUDE_PATH}
# #             ${CMAKE_SOURCE_DIR}/include
# #             ${CMAKE_CURRENT_SOURCE_DIR})
# # include_directories(${ZLIB_INCLUDE_DIRS})
# # target_link_libraries(MyExecutable ${ZLIB_LIBRARIES})
# target_link_libraries(
#     NeoFOAM_GPL
#     PUBLIC ${CMAKE_SOURCE_DIR}/NeoFOAM_GPL/lib/libt8.a
#             ${CMAKE_SOURCE_DIR}/NeoFOAM_GPL/lib/libsc.a
#             ${CMAKE_SOURCE_DIR}/NeoFOAM_GPL/lib/libp4est.a
#             MPI::MPI_CXX
#             ${ZLIB_LIBRARIES})

# target_sources(NeoFOAM_GPL PRIVATE "treeAMR/treeAMRMesh.cpp")

