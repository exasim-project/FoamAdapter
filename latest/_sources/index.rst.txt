.. FoamAdapter documentation master file, created by
   sphinx-quickstart on Sat Dec 16 15:22:11 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FoamAdapter!
===================

The FoamAdapter is the coupling interface of the FoamAdapter CFD core to standard OpenFOAM. It provides application examples to run OpenFOAM simulations via the FoamAdapter backend.


Table of Contents
^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   self
   gettingStarted
   testcases

Running OpenFOAM with FoamAdapter backend
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We are aiming for a high level of compatibility with OpenFOAM. However, we don't expect binary or ABI compatibility.
This means that FoamAdapter won't produce a `libfiniteVolume.so` and `libOpenFOAM.so` which could serve as a plugin replacement for existing `libfiniteVolume.so` and `libOpenFOAM.so`.
Instead, we aim for source compatibility, i.e. the possibility to compile application OpenFOAM code like pimpleFoam and others against the FoamAdapter libraries.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
