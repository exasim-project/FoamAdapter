// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "NeoN/core/database/oldTimeCollection.hpp"
#include "NeoN/dsl/expression.hpp"
#include "NeoN/dsl/solver.hpp"
#include "NeoN/dsl/ddt.hpp"
#include "FoamAdapter/readers/foamDictionary.hpp"

#include "NeoN/dsl/implicit.hpp"
#include "NeoN/dsl/explicit.hpp"


#include "FoamAdapter/FoamAdapter.hpp"
#include "FoamAdapter/readers/foamDictionary.hpp"

#include "common.hpp"


using Foam::Info;
using Foam::endl;
using Foam::nl;
namespace fvc = Foam::fvc;
namespace fvm = Foam::fvm;

namespace dsl = NeoN::dsl;
namespace fvcc = NeoN::finiteVolume::cellCentred;

extern Foam::Time* timePtr; // A single time object


TEST_CASE("PressureVelocityCoupling")
{
    Foam::Time& runTime = *timePtr;

    NeoN::Database db;
    fvcc::VectorCollection& VectorCollection =
        fvcc::VectorCollection::instance(db, "VectorCollection");

    NeoN::Executor exec = GENERATE(NeoN::Executor(NeoN::SerialExecutor {})
                                   // NeoN::Executor(NeoN::CPUExecutor {}),
                                   // NeoN::Executor(NeoN::GPUExecutor {})
    );

    std::string execName = std::visit([](auto e) { return e.name(); }, exec);

    auto meshPtr = Foam::createMesh(exec, runTime);
    Foam::MeshAdapter& mesh = *meshPtr;
    auto nfMesh = mesh.nfMesh();

    auto ofU = randomVectorField(runTime, mesh, "ofU");
    auto& oldOfU = ofU.oldTime();
    oldOfU.primitiveFieldRef() = Foam::vector(0.0, 0.0, 0.0);
    oldOfU.correctBoundaryConditions();

    Info << "creating NeoFOAM velocity fields" << endl;
    fvcc::VolumeField<NeoN::Vec3>& nfU =
        VectorCollection.registerVector<fvcc::VolumeField<NeoN::Vec3>>(
            Foam::CreateFromFoamField<Foam::volVectorField> {
                .exec = exec,
                .nfMesh = nfMesh,
                .foamField = ofU,
                .name = "nfU"
            }
        );

    auto& nfOldU = fvcc::oldTime(nfU);
    NeoN::fill(nfOldU.internalVector(), NeoN::Vec3(0.0, 0.0, 0.0));
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
    NeoN::fill(nfNu.boundaryVector().value(), 0.01);

    Foam::scalar t = runTime.time().value();
    Foam::scalar dt = runTime.deltaT().value();

    NeoN::Dictionary fvSchemesDict = Foam::readFoamDictionary(mesh.schemesDict());
    NeoN::Dictionary fvSolutionDict = Foam::readFoamDictionary(mesh.solutionDict());
    auto& solverDict = fvSolutionDict.get<NeoN::Dictionary>("solvers");

    SECTION("discreteMomentumFields " + execName)
    {
        fvcc::Expression<NeoN::Vec3> nfUEqn(
            dsl::imp::ddt(nfU) + dsl::imp::div(nfPhi, nfU) - dsl::imp::laplacian(nfNu, nfU),
            // dsl::imp::ddt(nfU) + dsl::imp::laplacian(nfNu, nfU),
            // dsl::imp::ddt(nfU) - dsl::imp::div(nfPhi, nfU),
            nfU,
            fvSchemesDict,
            solverDict.get<NeoN::Dictionary>("nfU")
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

            auto hostnfRAU = nfrAU.internalVector().copyToHost();
            write(nfrAU.internalVector(), mesh, "nfrAU");

            for (size_t celli = 0; celli < hostnfRAU.size(); celli++)
            {
                REQUIRE(hostnfRAU.view()[celli] == Catch::Approx(forAU[celli]).margin(1e-16));
            }
        }

        SECTION("rAU modified U")
        {
            ofU.primitiveFieldRef() *= 2.5;
            ofU.correctBoundaryConditions();
            nfU.internalVector() *= 2.5;
            nfU.correctBoundaryConditions();
            auto [nfrAU, nfHbyA] = fvcc::discreteMomentumFields(nfUEqn);
            Foam::volScalarField forAU("forAU", 1.0 / ofUEqn.A());
            forAU.write();

            auto hostnfRAU = nfrAU.internalVector().copyToHost();
            write(nfrAU.internalVector(), mesh, "nfrAU");

            for (size_t celli = 0; celli < hostnfRAU.size(); celli++)
            {
                REQUIRE(hostnfRAU.view()[celli] == Catch::Approx(forAU[celli]).margin(1e-16));
            }
        }

        SECTION("HbyA")
        {
            auto [nfrAU, nfHbyA] = fvcc::discreteMomentumFields(nfUEqn);
            Foam::volScalarField forAU("forAU", 1.0 / ofUEqn.A());
            Foam::volVectorField HbyA("HbyA", forAU * ofUEqn.H());
            HbyA.write();

            auto hostnfHbyA = nfHbyA.internalVector().copyToHost();
            write(nfHbyA.internalVector(), mesh, "nfHbyA");

            for (size_t celli = 0; celli < hostnfHbyA.size(); celli++)
            {
                REQUIRE(hostnfHbyA.view()[celli][0] == Catch::Approx(HbyA[celli][0]).margin(1e-14));
                REQUIRE(hostnfHbyA.view()[celli][1] == Catch::Approx(HbyA[celli][1]).margin(1e-14));
                REQUIRE(hostnfHbyA.view()[celli][2] == Catch::Approx(HbyA[celli][2]).margin(1e-14));
            }
        }

        SECTION("HbyA modified U")
        {
            ofU.primitiveFieldRef() *= 2.5;
            ofU.correctBoundaryConditions();
            nfU.internalVector() *= 2.5;
            nfU.correctBoundaryConditions();
            auto [nfrAU, nfHbyA] = fvcc::discreteMomentumFields(nfUEqn);
            Foam::volScalarField forAU("forAU", 1.0 / ofUEqn.A());
            Foam::volVectorField HbyA("HbyA", forAU * ofUEqn.H());
            HbyA.write();

            auto hostnfHbyA = nfHbyA.internalVector().copyToHost();
            write(nfHbyA.internalVector(), mesh, "nfHbyA");

            for (size_t celli = 0; celli < hostnfHbyA.size(); celli++)
            {
                REQUIRE(hostnfHbyA.view()[celli][0] == Catch::Approx(HbyA[celli][0]).margin(1e-14));
                REQUIRE(hostnfHbyA.view()[celli][1] == Catch::Approx(HbyA[celli][1]).margin(1e-14));
                REQUIRE(hostnfHbyA.view()[celli][2] == Catch::Approx(HbyA[celli][2]).margin(1e-14));
            }
        }
    }
}
