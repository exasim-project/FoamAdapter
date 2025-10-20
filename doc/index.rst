.. FoamAdapter documentation master file, created by
   sphinx-quickstart on Sat Dec 16 15:22:11 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FoamAdapter!
***********************

FoamAdapter an coupling interface, it provides platform-portable implementations of finite volume CFD algorithms and solvers using NeoN as a computational backend.
It leverages OpenFOAM for IO tasks such as reading OpenFOAM case files and configurations with the aim that OpenFOAM simulations can be run on accelerator devices.
Additionally to implementing a "classical" C++ interface for solvers like neoIcoFOAM, FoamAdapter provides a modern Python interface, which simplifies the development of new solvers and interfacing with other (python) libraries.

Features
--------

.. |TODO| unicode:: U+2610  .. ☐
.. |DONE| unicode:: U+2611  .. ☑

The current features of FoamAdapter include:

- |DONE| CLI to run tools and solvers
- |DONE| Support for OpenFOAM standard case structure
- |DONE| Support for common OpenFOAM solvers (icoFoam, pimpleFoam) via python interface
- |DONE| Simple solver extension
- |DONE| Unit tests written
- |TODO| GPU acceleration via NeoN
- |TODO| Multiphysics support
- |TODO| external coupling interface (e.g via precise)

Table of Contents
*****************
.. toctree::
   :maxdepth: 2

   self
   installation
   testcases

.. toctree::
   :maxdepth: 2
   :caption: Usage:

   usage/quickstart
   usage/cli

.. toctree::
   :maxdepth: 2
   :caption: Development:

   development/architecture
   development/contributing
   ci

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Indices and tables
==================
** :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
