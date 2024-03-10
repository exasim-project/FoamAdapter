# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors
# include(FetchContent)
# FetchContent_Declare(
#   googletest
#   URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
# )
# FetchContent_MakeAvailable(googletest)

# include(FetchContent)
# FetchContent_Declare(
#   googlebenchmark
#   URL https://github.com/google/benchmark/archive/refs/tags/v1.8.3.zip
# )

# FetchContent_MakeAvailable(googlebenchmark)

add_subdirectory(explicitOperators)
add_subdirectory(matrixAssembly)
add_subdirectory(gradOperator)
add_subdirectory(scalarAdvection)
