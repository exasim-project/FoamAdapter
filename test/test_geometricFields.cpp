// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include "NeoFOAM/fields/field.hpp"

#include "NeoFOAM/fields/boundaryFields.hpp"
#include "NeoFOAM/fields/domainField.hpp"
#include "NeoFOAM/finiteVolume/cellCentred.hpp"

#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenGrad.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/interpolation/linear.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/interpolation/upwind.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/interpolation/surfaceInterpolation.hpp"

#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/writers/writers.hpp"

#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"
#include "FoamAdapter/setup/setup.hpp"

Foam::Time* timePtr;    // A single time object
Foam::argList* argsPtr; // Some forks want argList access at createMesh.H
Foam::fvMesh* meshPtr;  // A single mesh object

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
    std::unique_ptr<Foam::fvMesh> meshPtr = Foam::createMesh(runTime);
    Foam::fvMesh& mesh = *meshPtr;
    argsPtr = &args;
    timePtr = &runTime;

    int result = session.run();

    // Run benchmarks if there are any
    Kokkos::finalize();

    return result;
}

template<typename ValueType>
void checkField(const NeoFOAM::Field<ValueType>& field, ValueType value)
{
    auto field_host = field.copyToHost().span();
    for (int i = 0; i < field_host.size(); i++)
    {
        REQUIRE(field_host[i] == value);
    }
}

TEST_CASE("fvcc::VolumeField")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;
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

        NeoFOAM::fvcc::VolumeField<NeoFOAM::scalar> neoT = constructFrom(exec, uMesh, T);

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

        NeoFOAM::fvcc::VolumeField<NeoFOAM::Vector> neoU = constructFrom(exec, uMesh, U);

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
