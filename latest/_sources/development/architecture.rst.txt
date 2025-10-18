Architecture
============

This document describes the overall architecture of FoamAdapter, including both the C++ core and Python interface components.

Overview
--------

FoamAdapter is designed as a multi-layered architecture that bridges OpenFOAM and NeoN computational backends while providing both C++ and Python interfaces.

.. mermaid::

   graph TB
       A[Python CLI/API] --> B[foamadapter Python Module]
       B --> C[pybFoam Bindings]
       C --> D[FoamAdapter C++ Core]
       D --> E[NeoN Backend]
       D --> F[OpenFOAM Libraries]
       E --> G[GPU/CPU Execution]
       F --> H[OpenFOAM Solvers]
