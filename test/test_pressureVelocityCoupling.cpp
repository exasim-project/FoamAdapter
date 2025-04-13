// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "NeoFOAM/dsl/expression.hpp"
#include "NeoFOAM/dsl/solver.hpp"
#include "NeoFOAM/dsl/ddt.hpp"
#include "FoamAdapter/readers/foamDictionary.hpp"

#include "NeoFOAM/dsl/implicit.hpp"
#include "NeoFOAM/dsl/explicit.hpp"


#include "FoamAdapter/FoamAdapter.hpp"
#include "FoamAdapter/readers/foamDictionary.hpp"

#include "common.hpp"


using Foam::Info;
using Foam::endl;
using Foam::nl;
namespace fvc = Foam::fvc;
namespace fvm = Foam::fvm;

namespace dsl = NeoFOAM::dsl;
namespace fvcc = NeoFOAM::finiteVolume::cellCentred;

extern Foam::Time* timePtr; // A single time object


TEST_CASE("PressureVelocityCoupling")
{
    Foam::Time& runTime = *timePtr;

    NeoFOAM::Database db;
    fvcc::FieldCollection& fieldCollection = fvcc::FieldCollection::instance(db, "fieldCollection");

    NeoFOAM::Executor exec = GENERATE(NeoFOAM::Executor(NeoFOAM::SerialExecutor {})
                                      // NeoFOAM::Executor(NeoFOAM::CPUExecutor {}),
                                      // NeoFOAM::Executor(NeoFOAM::GPUExecutor {})
    );

    std::string execName = std::visit([](auto e) { return e.name(); }, exec);

    auto meshPtr = Foam::createMesh(exec, runTime);
    Foam::MeshAdapter& mesh = *meshPtr;
    auto nfMesh = mesh.nfMesh();

    auto ofU = randomVectorField(runTime, mesh, "ofU");
    auto& oldOfU = ofU.oldTime();
    oldOfU.primitiveFieldRef() = Foam::Vector(0.0, 0.0, 0.0);
    oldOfU.correctBoundaryConditions();

    Info << "creating NeoFOAM velocity fields" << endl;
    fvcc::VolumeField<NeoFOAM::Vector>& nfU =
        fieldCollection.registerField<fvcc::VolumeField<NeoFOAM::Vector>>(
            Foam::CreateFromFoamField<Foam::volVectorField> {
                .exec = exec,
                .nfMesh = nfMesh,
                .foamField = ofU,
                .name = "nfU"
            }
        );

    auto& nfOldU = fvcc::oldTime(nfU);
    NeoFOAM::fill(nfOldU.internalField(), NeoFOAM::Vector(0.0, 0.0, 0.0));
    nfOldU.correctBoundaryConditions();

    Foam::surfaceScalarField ofPhi(
        Foam::IOobject(
            "ofPhi",
            runTime.timeName(),
            mesh,
            Foam::IOobject::NO_READ,
            Foam::IOobject::NO_WRITE
        ),
        fvc::flux(ofU)
    );
    auto nfPhi = constructSurfaceField(exec, nfMesh, ofPhi);
    nfPhi.name = "nfPhi";

    Foam::surfaceScalarField ofNu(
        Foam::IOobject(
            "ofNu",
            runTime.timeName(),
            mesh,
            Foam::IOobject::NO_READ,
            Foam::IOobject::NO_WRITE
        ),
        mesh,
        Foam::dimensionedScalar("ofNu", Foam::dimensionSet(0, 2, -1, 0, 0), 0.01)
    );

    Info << "creating NeoFOAM nu fields" << endl;
    Info << "ofNu: " << ofNu << endl;
    auto nfNu = constructSurfaceField(exec, nfMesh, ofNu);
    nfNu.name = "nfNu";
    NeoFOAM::fill(nfNu.boundaryField().value(), 0.01);

    Foam::scalar t = runTime.time().value();
    Foam::scalar dt = runTime.deltaT().value();

    NeoFOAM::Dictionary fvSchemesDict = Foam::readFoamDictionary(mesh.schemesDict());
    NeoFOAM::Dictionary fvSolutionDict = Foam::readFoamDictionary(mesh.solutionDict());
    auto& solverDict = fvSolutionDict.get<NeoFOAM::Dictionary>("solvers");

    SECTION("discreteMomentumFields " + execName)
    {
        fvcc::Expression<NeoFOAM::Vector> nfUEqn(
            dsl::imp::ddt(nfU) + dsl::imp::div(nfPhi, nfU) - dsl::imp::laplacian(nfNu, nfU),
            // dsl::imp::ddt(nfU) + dsl::imp::laplacian(nfNu, nfU),
            // dsl::imp::ddt(nfU) - dsl::imp::div(nfPhi, nfU),
            nfU,
            fvSchemesDict,
            solverDict.get<NeoFOAM::Dictionary>("nfU")
        );

        nfUEqn.assemble(t, dt);


        Foam::fvVectorMatrix ofUEqn(
            fvm::ddt(ofU) + fvm::div(ofPhi, ofU) - fvm::laplacian(ofNu, ofU)
        );
        // Foam::fvVectorMatrix ofUEqn(fvm::ddt(ofU) + fvm::laplacian(ofNu, ofU));
        // Foam::fvVectorMatrix ofUEqn(fvm::ddt(ofU) - fvm::div(ofPhi, ofU));


        SECTION("rAU")
        {

            auto [nfrAU, nfHbyA] = fvcc::discreteMomentumFields(nfUEqn);
            Foam::volScalarField forAU("forAU", 1.0 / ofUEqn.A());
            forAU.write();

            auto hostnfRAU = nfrAU.internalField().copyToHost();
            write(nfrAU.internalField(), mesh, "nfrAU");

            for (size_t celli = 0; celli < hostnfRAU.size(); celli++)
            {
                REQUIRE(hostnfRAU[celli] == Catch::Approx(forAU[celli]).margin(1e-16));
            }
        }

        SECTION("rAU modified U")
        {
            ofU.primitiveFieldRef() *= 2.5;
            ofU.correctBoundaryConditions();
            nfU.internalField() *= 2.5;
            nfU.correctBoundaryConditions();
            auto [nfrAU, nfHbyA] = fvcc::discreteMomentumFields(nfUEqn);
            Foam::volScalarField forAU("forAU", 1.0 / ofUEqn.A());
            forAU.write();

            auto hostnfRAU = nfrAU.internalField().copyToHost();
            write(nfrAU.internalField(), mesh, "nfrAU");

            for (size_t celli = 0; celli < hostnfRAU.size(); celli++)
            {
                REQUIRE(hostnfRAU[celli] == Catch::Approx(forAU[celli]).margin(1e-16));
            }
        }

        SECTION("HbyA")
        {
            auto [nfrAU, nfHbyA] = fvcc::discreteMomentumFields(nfUEqn);
            Foam::volScalarField forAU("forAU", 1.0 / ofUEqn.A());
            Foam::volVectorField HbyA("HbyA", forAU * ofUEqn.H());
            HbyA.write();

            auto hostnfHbyA = nfHbyA.internalField().copyToHost();
            write(nfHbyA.internalField(), mesh, "nfHbyA");

            for (size_t celli = 0; celli < hostnfHbyA.size(); celli++)
            {
                REQUIRE(hostnfHbyA[celli][0] == Catch::Approx(HbyA[celli][0]).margin(1e-14));
                REQUIRE(hostnfHbyA[celli][1] == Catch::Approx(HbyA[celli][1]).margin(1e-14));
                REQUIRE(hostnfHbyA[celli][2] == Catch::Approx(HbyA[celli][2]).margin(1e-14));
            }
        }

        SECTION("HbyA modified U")
        {
            ofU.primitiveFieldRef() *= 2.5;
            ofU.correctBoundaryConditions();
            nfU.internalField() *= 2.5;
            nfU.correctBoundaryConditions();
            auto [nfrAU, nfHbyA] = fvcc::discreteMomentumFields(nfUEqn);
            Foam::volScalarField forAU("forAU", 1.0 / ofUEqn.A());
            Foam::volVectorField HbyA("HbyA", forAU * ofUEqn.H());
            HbyA.write();

            auto hostnfHbyA = nfHbyA.internalField().copyToHost();
            write(nfHbyA.internalField(), mesh, "nfHbyA");

            for (size_t celli = 0; celli < hostnfHbyA.size(); celli++)
            {
                REQUIRE(hostnfHbyA[celli][0] == Catch::Approx(HbyA[celli][0]).margin(1e-14));
                REQUIRE(hostnfHbyA[celli][1] == Catch::Approx(HbyA[celli][1]).margin(1e-14));
                REQUIRE(hostnfHbyA[celli][2] == Catch::Approx(HbyA[celli][2]).margin(1e-14));
            }
        }
    }
}
