Command Line Interface
======================

The FoamAdapter CLI provides a modern, user-friendly interface to OpenFOAM solvers with enhanced input validation, error reporting, and workflow automation.

Overview
--------

The CLI is built using Typer and provides:

* **Input validation** with detailed error messages
* **Structured command hierarchy** for different solver types
* **Integration with pybFoam** for seamless OpenFOAM interoperability
* **Pydantic-based configuration** for type safety and validation

Basic Usage
-----------

Main Command
~~~~~~~~~~~~

.. code-block:: bash

    foamadapter [OPTIONS] COMMAND [ARGS]...

Global Options
~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Option
     - Description
   * - ``--help``
     - Show help message and exit
   * - ``--version``
     - Show version information

Available Commands
~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 15 85
   :header-rows: 1

   * - Command
     - Description
   * - ``solver``
     - Run CFD solvers (icofoam, pimplefoam)

Solver Commands
---------------

The solver subcommand provides access to various CFD solvers.

Solver Command Structure
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    foamadapter solver [SOLVER_NAME] [OPTIONS]

Available Solvers
~~~~~~~~~~~~~~~~~

icoFoam Solver
^^^^^^^^^^^^^^

Transient incompressible Navier-Stokes solver for laminar flow.

.. code-block:: bash

    foamadapter solver icofoam [OPTIONS]

**Options:**

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Option
     - Description
   * - ``--help``
     - Show icoFoam specific help
   * - Standard OpenFOAM options
     - All standard icoFoam options are supported

Input Validation
----------------

Enhanced Error Checking
~~~~~~~~~~~~~~~~~~~~~~~

The ``--check_inputs`` flag provides comprehensive validation of your OpenFOAM case setup.


Validation Examples
~~~~~~~~~~~~~~~~~~~

**Successful Validation:**

.. code-block:: bash

    $ foamadapter solver pimplefoam --check_inputs
    All inputs are valid!

**Failed Validation:**

.. code-block:: bash

    $ foamadapter solver pimplefoam --check_inputs
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

    error in 0/U:
        error type     = boundary_mismatch
        affected key   = boundaryField/inlet/type
        error message  = boundary condition type not supported for this field
        provided value = "wrongType"

    Fix the above errors before running the solver.
