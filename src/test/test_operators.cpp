// SPDX-License-Identifier: MPL-2.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include "NeoFOAM/fields/Field.hpp"
#include "NeoFOAM/fields/FieldOperations.hpp"
#include "NeoFOAM/fields/FieldTypeDefs.hpp"

#include "NeoFOAM/fields/boundaryFields.hpp"
#include "NeoFOAM/fields/domainField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/fields/fvccVolField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/fvccBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/scalar/fvccScalarFixedValueBoundaryField.hpp"

#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/writers/writers.hpp"

#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"

Foam::Time* timePtr;      // A single time object
Foam::argList* argsPtr;   // Some forks want argList access at createMesh.H
Foam::fvMesh* meshPtr;    // A single mesh object

int main(int argc, char* argv[])
{

    // Initialize Catch2
    Kokkos::initialize(argc, argv);
    Catch::Session session;

    // Specify command line options
    int returnCode = session.applyCommandLine(argc, argv);
    if (returnCode != 0) // Indicates a command line error
        return returnCode;

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    argsPtr = &args;
    timePtr = &runTime;

    int result = session.run();

    // Run benchmarks if there are any
    Kokkos::finalize();

    return result;
}

TEST_CASE("GradOperator")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;
    #include "createMesh.H"
    NeoFOAM::executor exec = NeoFOAM::CPUExecutor();
    Foam::Info << "reading mesh" << Foam::endl;
    NeoFOAM::unstructuredMesh uMesh = readOpenFOAMMesh(exec, mesh);

    Foam::Info<< "Reading field T\n" << Foam::endl;

    Foam::volScalarField T
    (
        Foam::IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::AUTO_WRITE
        ),
        mesh
    );

    Foam::volVectorField ofGradT = Foam::fvc::grad(T);
    
    


    // NeoFOAM::fvccVolField<NeoFOAM::scalarField> gradT(exec, uMesh, T.internalField().size());

    // NeoFOAM::vectorField nofGradT = NeoFOAM::gaussGreenGrad(exec, uMesh).grad(Temperature);
    
}
