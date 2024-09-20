**[Requirements](#requirements)** |
**[Compilation](#Compilation)** |
# FoamAdapter

## Requirements

FoamAdapter has the following requirements

*  _cmake 3.22+_
*  _gcc >= 10_ or  _clang >= 17+_

## Compilation

We provide several Cmake presets to set commmonly required flags if you compile FoamAdapter

    cmake --list-presets # To list existing presets
    cmake --preset production # To compile for production use

## run

export PATH=$PWD/FoamAdapter/bin:$PATH
and go into benchmarks/
