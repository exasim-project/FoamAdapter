.. FoamAdapter documentation master file, created by
   sphinx-quickstart on Sat Dec 16 15:22:11 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FoamAdapter!
===================

The FoamAdapter is the coupling interface. It provides platform-portable implementations of common CFD algorithms and solvers using NeoN as a computational backend. It leverages standard OpenFOAM such that OpenFOAM simulations can be run on accelerator devices.

Goals
^^^^^

The Goal of FoamAdapter is to simply the usage of OpenFOAM for developers and users by providing a modern Python interface while leveraging the heterogeneous architecture.

This requires in our view a incremental rewrite of OpenFOAM solvers. The python interface improves the development experience and lowers the barrier of entry for engineers and scientists to implement their own CFD solvers. Additionally, the usage of python opens a new ecosysstem of tools for pre- and postprocessing.

Features
^^^^^^^^

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
^^^^^^^^^^^^^^^^^
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
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
