# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023-2024 FoamAdapter authors

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.4.0 # or a later release
)
FetchContent_MakeAvailable(Catch2)

include(CTest)
# include(Catch)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

add_library(foamadapter_catch_main test_main.cpp)
target_link_libraries(foamadapter_catch_main PUBLIC Catch2::Catch2 OpenFOAM FoamAdapter)
install(TARGETS foamadapter_catch_main Catch2)
