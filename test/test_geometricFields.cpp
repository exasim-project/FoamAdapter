// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "catch2/common.hpp"

#include "NeoFOAM/fields/field.hpp"

#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/writers/writers.hpp"

#define namespaceFoam // Suppress <using namespace Foam;>
#include "FoamAdapter/setup/setup.hpp"

namespace fvcc = NeoFOAM::finiteVolume::cellCentred;

extern Foam::Time* timePtr;   // A single time object
extern Foam::fvMesh* meshPtr; // A single mesh object


TEST_CASE("fvcc::VolumeField")
{
    Foam::Time& runTime = *timePtr;
    std::unique_ptr<Foam::fvMesh> meshPtr = Foam::createMesh(runTime);
    Foam::fvMesh& mesh = *meshPtr;

    Foam::volScalarField T(
        Foam::IOobject(
            "T", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
        ),
        mesh
    );

    Foam::volVectorField U(
        Foam::IOobject(
            "U", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
        ),
        mesh
    );

    NeoFOAM::Executor exec = GENERATE(
        NeoFOAM::Executor(NeoFOAM::CPUExecutor {}),
        NeoFOAM::Executor(NeoFOAM::SerialExecutor {}),
        NeoFOAM::Executor(NeoFOAM::GPUExecutor {})
    );
    std::string exec_name = std::visit([](auto e) { return e.print(); }, exec);


    SECTION("fvccVolField_[scalar]" + exec_name)
    {
        Foam::Info << "reading mesh with executor: " << exec_name << Foam::endl;
        NeoFOAM::UnstructuredMesh uMesh = readOpenFOAMMesh(exec, mesh);

        fvcc::VolumeField<NeoFOAM::scalar> neoT = constructFrom(exec, uMesh, T);

        REQUIRE(neoT.internalField().size() == T.internalField().size());
        fill(neoT.internalField(), 1.0);
        checkField(neoT.internalField(), 1.0);
        neoT.correctBoundaryConditions();
        checkField(neoT.boundaryField().value(), 1.0);

        fill(neoT.internalField(), 2.0);
        neoT.correctBoundaryConditions();
        checkField(neoT.boundaryField().value(), 2.0);
    }

    SECTION("fvccVolField_[vector]" + exec_name)
    {
        Foam::Info << "reading mesh with executor: " << exec_name << Foam::endl;
        NeoFOAM::UnstructuredMesh uMesh = readOpenFOAMMesh(exec, mesh);

        fvcc::VolumeField<NeoFOAM::Vector> neoU = constructFrom(exec, uMesh, U);

        REQUIRE(neoU.internalField().size() == U.internalField().size());
        fill(neoU.internalField(), NeoFOAM::Vector(1.0, 1.0, 1.0));
        checkField(neoU.internalField(), NeoFOAM::Vector(1.0, 1.0, 1.0));
        neoU.correctBoundaryConditions();
        checkField(neoU.boundaryField().value(), NeoFOAM::Vector(1.0, 1.0, 1.0));

        fill(neoU.internalField(), NeoFOAM::Vector(2.0, 2.0, 2.0));
        neoU.correctBoundaryConditions();
        checkField(neoU.boundaryField().value(), NeoFOAM::Vector(2.0, 2.0, 2.0));
    }
}
