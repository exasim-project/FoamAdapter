Goals
=====

The Goal of FoamAdapter is to simply the usage of OpenFOAM for developers and users by providing a modern Python interface while leveraging the heterogeneous architecture.

This requires in our view a incremental rewrite of OpenFOAM solvers. The python interface improves the development experience and lowers the barrier of entry for engineers and scientists to implement their own CFD solvers. Additionally, the usage of python opens a new ecosysstem of tools for pre- and postprocessing.

Features
========


.. |TODO| unicode:: U+2610  .. ☐
.. |DONE| unicode:: U+2611  .. ☑

The current features of FoamAdapter include:

- |TODO| Multiphysics support
- |TODO| external coupling interface (e.g via precise)
- |TODO| GPU acceleration via NeoN
- |DONE| CLI to run tools and solvers
- |DONE| Support for OpenFOAM standard case structure
- |DONE| Support for common OpenFOAM solvers (icoFoam, pimpleFoam) via python interface
- |DONE| Simple solver extension
- |DONE| Unit tests written
