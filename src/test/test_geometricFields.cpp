// SPDX-License-Identifier: GPL-3.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include "NeoFOAM/fields/Field.hpp"
#include "NeoFOAM/fields/FieldOperations.hpp"
#include "NeoFOAM/fields/FieldTypeDefs.hpp"
#include "NeoFOAM/fields/comparisions/fieldComparision.hpp"

#include "NeoFOAM/fields/boundaryFields.hpp"
#include "NeoFOAM/fields/domainField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/fields/fvccVolField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/fields/fvccSurfaceField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/fvccBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/vol/scalar/fvccScalarFixedValueBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/vol/scalar/fvccScalarZeroGradientBoundaryField.hpp"

#include "NeoFOAM/cellCentredFiniteVolume/grad/gaussGreenGrad.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/surfaceInterpolation/linear.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/surfaceInterpolation/upwind.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/surfaceInterpolation/surfaceInterpolation.hpp"

#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/writers/writers.hpp"

#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"

Foam::Time *timePtr;    // A single time object
Foam::argList *argsPtr; // Some forks want argList access at createMesh.H
Foam::fvMesh *meshPtr;  // A single mesh object

int main(int argc, char *argv[])
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

TEST_CASE("fvccVolField")
{
    Foam::Time &runTime = *timePtr;
    Foam::argList &args = *argsPtr;
#include "createMesh.H"

    Foam::volScalarField T(
        Foam::IOobject(
            "T",
            runTime.timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::AUTO_WRITE),
        mesh);

    Foam::volVectorField U
    (
        Foam::IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::AUTO_WRITE
        ),
        mesh
    );

    NeoFOAM::executor exec = GENERATE(
        NeoFOAM::executor(NeoFOAM::CPUExecutor {}),
        NeoFOAM::executor(NeoFOAM::OMPExecutor {}),
        NeoFOAM::executor(NeoFOAM::GPUExecutor {})
    );
    std::string exec_name = std::visit([](auto e) { return e.print(); }, exec);


    SECTION("fvccVolField_[scalar]" + exec_name)
    {
        
        Foam::Info << "reading mesh with executor: " << exec_name << Foam::endl;
        NeoFOAM::unstructuredMesh uMesh = readOpenFOAMMesh(exec, mesh);

        NeoFOAM::fvccVolField<NeoFOAM::scalar> neoT = constructFrom(exec, uMesh, T);

        fill(neoT.internalField(), 1.0);
        REQUIRE(compare(neoT.internalField(), 1.0));
        neoT.correctBoundaryConditions();
        REQUIRE(compare(neoT.boundaryField().value(), 1.0));

        fill(neoT.internalField(), 2.0);
        neoT.correctBoundaryConditions();
        REQUIRE(compare(neoT.boundaryField().value(), 2.0));
    }

    SECTION("fvccVolField_[vector]" + exec_name)
    {
        
        Foam::Info << "reading mesh with executor: " << exec_name << Foam::endl;
        NeoFOAM::unstructuredMesh uMesh = readOpenFOAMMesh(exec, mesh);

        NeoFOAM::fvccVolField<NeoFOAM::Vector> neoU = constructFrom(exec, uMesh, U);

        fill(neoU.internalField(), NeoFOAM::Vector(1.0, 1.0, 1.0));
        REQUIRE(compare(neoU.internalField(), NeoFOAM::Vector(1.0, 1.0, 1.0)));
        neoU.correctBoundaryConditions();
        REQUIRE(compare(neoU.boundaryField().value(), NeoFOAM::Vector(1.0, 1.0, 1.0)));

        fill(neoU.internalField(), NeoFOAM::Vector(2.0, 2.0, 2.0));
        neoU.correctBoundaryConditions();
        REQUIRE(compare(neoU.boundaryField().value(), NeoFOAM::Vector(2.0, 2.0, 2.0)));
    }
}
