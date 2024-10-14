**[Requirements](#requirements)** |
**[Compilation](#Compilation)** |
# FoamAdapter

This repository showcases how NeoFOAM can be used in combination with OpenFOAM.
It provides converters between OpenFOAM and NeoFOAM datastructures.

The repository is structured in the following way:
- src and include implement common functionality to copy data between OpenFOAM and NeoFOAM
- tests demonstrating that NeoFOAM and OpenFOAM deliver identical results are provided by this repository in the test folder.

## Requirements

FoamAdapter has the following requirements

*  _cmake 3.22+_
*  _gcc >= 10_ or  _clang >= 17+_
* OpenFOAM _2406_

## Compilation

We provide several Cmake presets to set commmonly required flags if you compile FoamAdapter

    cmake --list-presets # To list existing presets
    cmake --preset production # To compile for production use

## run

export PATH=$PWD/FoamAdapter/bin:$PATH
and go into benchmarks/
