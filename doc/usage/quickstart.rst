Quick Start Guide
=================

This guide will get you up and running with FoamAdapter Python module in minutes.
We'll walk through running your first CFD simulation using the Python interface.

Prerequisites
-------------

Before starting, ensure you have:

* FoamAdapter installed (see :doc:`installation`)
* OpenFOAM 2406+ properly sourced
* A basic understanding of OpenFOAM case structure

Your First Simulation
---------------------

Let's run the classic lid-driven cavity case using FoamAdapter's Python interface.

Step 1: Set Up a Case Directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a new directory for your case:

.. code-block:: bash

    mkdir my_cavity_case
    cd my_cavity_case

Step 2: Copy Tutorial Case
~~~~~~~~~~~~~~~~~~~~~~~~~~

Copy the cavity tutorial case (if available):

.. code-block:: bash

    # From FoamAdapter tutorials
    cd test/solver/pitzDaily


Step 3: Run the Solver
~~~~~~~~~~~~~~~~~~~~~~

Use FoamAdapter's CLI to run the icoFoam solver:

.. code-block:: bash

    foamadapter solver pimplefoam

That's it! Your first CFD simulation with FoamAdapter is running.

Understanding the CLI
---------------------

Basic Commands
~~~~~~~~~~~~~~

FoamAdapter provides several command-line tools:

.. code-block:: bash

    # Show all available commands
    foamadapter --help

    # Show solver options
    foamadapter solver --help

    # Get help for specific solver
    foamadapter solver icofoam --help

Available Solvers
~~~~~~~~~~~~~~~~~

Currently supported solvers:

* **icofoam**: Transient incompressible flow solver
* **pimplefoam**: Transient incompressible flow solver with PIMPLE algorithm

Running with Options
~~~~~~~~~~~~~~~~~~~~

You can pass additional options to solvers:

.. code-block:: bash

    # Check inputs before running
    foamadapter solver pimplefoam --check_inputs

    # Specify case directory
    foamadapter solver pimplefoam --case /path/to/case

Input Validation
----------------

FoamAdapter includes powerful input validation using Pydantic models.

Validating Your Case
~~~~~~~~~~~~~~~~~~~~

Before running expensive simulations, validate your inputs:

.. code-block:: bash

    foamadapter solver pimplefoam --check_inputs --case .

This will check:

* **controlDict** parameters (time settings, Courant number, etc.)
* **fvSchemes** configuration
* **fvSolution** settings
* **Boundary conditions** consistency
* **Transport properties**

Example Validation Output
~~~~~~~~~~~~~~~~~~~~~~~~~

If validation fails, you'll see detailed error messages:

.. code-block:: text

    Input validation failed:
    error in controlDict:
        error type     = validation_error
        affected key   = maxCo
        error message  = ensure this value is greater than 0
        provided value = -0.5

    error in fvSchemes:
        error type     = missing_key
        affected key   = divSchemes/div(phi,U)
        error message  = required field missing
        provided value = None

    Fix the above errors before running the solver.



Running the Complete Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Create mesh (assuming blockMeshDict exists)
    blockMesh

    # Run with input validation
    foamadapter solver icofoam --check_inputs

    # If validation passes, run the solver
    foamadapter solver icofoam
