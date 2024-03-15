#!/bin/bash

# Run CMake to generate the build files
# -DCMAKE_INSTALL_PREFIX=$PWD/NeoFOAM_GPL \
cmake  -S . -B build  -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=$PWD/NeoFOAM_GPL \
        -DNEOFOAM_BUILD_TESTS=ON \
        -DKokkos_ENABLE_SERIAL=ON \
        -DKokkos_ENABLE_OPENMP=ON \
        -DKokkos_ENABLE_CUDA=ON
        # -DKokkos_ARCH_NATIVE=ON \
        # -DKokkos_ENABLE_AGGRESSIVE_VECTORIZATION=ON

# Build the project using make
cmake --build build -j
# cmake --install build
