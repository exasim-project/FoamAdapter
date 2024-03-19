#!/bin/bash
PROJ_DIR=$PWD
git clone https://github.com/DLR-AMR/t8code
cd t8code
git checkout v1.6.1
git submodule init
git submodule update
./bootstrap
mkdir t8code_build
cd t8code_build
../configure CFLAGS="-O3" CXXFLAGS="-O3" --enable-mpi CC=mpicc CXX=mpicxx --prefix=$PROJ_DIR/FoamAdapter
make -j
make install
cd ..