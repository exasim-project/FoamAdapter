#!/bin/bash

# Run CMake to generate the build files
cmake  -S . -B build  -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=$PWD/NeoFOAM_GPL \
        -DKokkos_ENABLE_SERIAL=ON \
        -DKokkos_ENABLE_OPENMP=ON \
        -DKokkos_ENABLE_CUDA=ON

# Build the project using make
cmake --build build --target install -j

PATH=$PWD/NeoFOAM_GPL/bin:$PATH