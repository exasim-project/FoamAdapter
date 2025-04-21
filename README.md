**[Requirements](#requirements)** |
**[Compilation](#Compilation)** |
**[Documentation](https://exasim-project.com/FoamAdapter/latest)** |
# NeoFOAM

NeoFOAM is a re-implementation of common OpenFOAM algorithms and solvers using
[NeoN](https://github.com/exasim-project/NeoN) as a computational backend.
It provides converters between OpenFOAM and NeoN datastructures, examples, benchmarks and tests.

## Requirements

NeoFOAM has the following requirements

*  _cmake 3.22+_
*  _gcc >= 12_ or  _clang >= 18+_
* OpenFOAM _2406_+
* CUDA  _12.1_ (for GPU support)

## Compilation

We provide several Cmake presets to set commmonly required flags for building NeoFOAM

    cmake --list-presets # To list existing presets
    cmake --preset production # To build for production use

### NeoN as a submodule

This repository depends on NeoN which is automatically downloaded via the CPM package manager.
Alternatively, NeoN can be included as git submodule for development purposes.
To do so the following steps can be used to initialise the submodule.

    git submodule update --init --recursive

## Structure

The repository is structured in the following way:
- src and include implement common functionality to copy data between OpenFOAM and NeoN
- tests demonstrating that NeoFOAM and OpenFOAM deliver identical results are provided by this repository in the test folder.
- examples provides examples of how NeoFOAM and NeoFOAM can be used for writing applications
- tutorials provides tutorial cases which can be run like typical OpenFOAM cases
