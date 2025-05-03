// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 FoamAdapter authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "common.hpp"
#include "constrainHbyA.H"

using Foam::Info;
using Foam::endl;
using Foam::nl;
namespace fvc = Foam::fvc;
namespace fvm = Foam::fvm;

namespace dsl = NeoN::dsl;
namespace nnfvcc = NeoN::finiteVolume::cellCentred;
namespace nffvcc = FoamAdapter;

extern Foam::Time* timePtr; // A single time object


TEST_CASE("PressureVelocityCoupling")
{
    Foam::Time& runTime = *timePtr;

    NeoN::Database db;
    auto& VectorCollection = nnfvcc::VectorCollection::instance(db, "VectorCollection");

    auto [execName, exec] = GENERATE(allAvailableExecutor());

    auto meshPtr = FoamAdapter::createMesh(exec, runTime);
    FoamAdapter::MeshAdapter& mesh = *meshPtr;
    auto nfMesh = mesh.nfMesh();

    auto ofU = randomVectorField(runTime, mesh, "ofU");
    auto ofp = randomScalarField(runTime, mesh, "ofp");
    // a predictable field is simpler to debug
    // forAll(ofU, celli)
    // {
    //     ofU[celli] = Foam::vector(celli, celli, celli);
    // }
    ofU.correctBoundaryConditions();
    auto& oldOfU = ofU.oldTime();
    oldOfU.primitiveFieldRef() = Foam::vector(0.0, 0.0, 0.0);
    oldOfU.correctBoundaryConditions();

    Info << "creating FoamAdapter velocity fields" << endl;
    nnfvcc::VolumeField<NeoN::Vec3>& nfU =
        VectorCollection.registerVector<nnfvcc::VolumeField<NeoN::Vec3>>(
            FoamAdapter::CreateFromFoamField<Foam::volVectorField> {
                .exec = exec,
                .nfMesh = nfMesh,
                .foamField = ofU,
                .name = "nfU"
            }
        );
    Info << "creating FoamAdapter pressure fields" << endl;
    auto nfp = VectorCollection.registerVector<nnfvcc::VolumeField<NeoN::scalar>>(
        FoamAdapter::CreateFromFoamField<Foam::volScalarField> {
            .exec = exec,
            .nfMesh = nfMesh,
            .foamField = ofp,
            .name = "nfp"
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
    auto nfPhi = FoamAdapter::constructSurfaceField(exec, nfMesh, ofPhi);
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

    Info << "creating FoamAdapter nu fields" << endl;
    Info << "ofNu: " << ofNu << endl;
    auto nfNu = FoamAdapter::constructSurfaceField(exec, nfMesh, ofNu);
    nfNu.name = "nfNu";
    NeoN::fill(nfNu.boundaryData().value(), 0.01);

    Foam::scalar t = runTime.time().value();
    Foam::scalar dt = runTime.deltaT().value();

    NeoN::Dictionary fvSchemesDict = FoamAdapter::convert(mesh.schemesDict());
    NeoN::Dictionary fvSolutionDict = FoamAdapter::convert(mesh.solutionDict());
    auto& solverDict = fvSolutionDict.get<NeoN::Dictionary>("solvers");

    SECTION("discreteMomentumFields " + execName)
    {
        nffvcc::Expression<NeoN::Vec3> nfUEqn(
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

            auto [nfrAU, nfHbyA] = nffvcc::discreteMomentumFields(nfUEqn);
            Foam::volScalarField forAU("forAU", 1.0 / ofUEqn.A());
            forAU.write();

            auto hostnfRAU = nfrAU.internalVector().copyToHost();
            write(nfrAU.internalVector(), mesh, "nfrAU" + execName);

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
            auto [nfrAU, nfHbyA] = nffvcc::discreteMomentumFields(nfUEqn);
            Foam::volScalarField forAU("forAU", 1.0 / ofUEqn.A());
            forAU.write();

            auto hostnfRAU = nfrAU.internalVector().copyToHost();
            write(nfrAU.internalVector(), mesh, "nfrAU" + execName);

            for (size_t celli = 0; celli < hostnfRAU.size(); celli++)
            {
                REQUIRE(hostnfRAU.view()[celli] == Catch::Approx(forAU[celli]).margin(1e-16));
            }
        }

        SECTION("HbyA")
        {
            auto [nfrAU, nfHbyA] = nffvcc::discreteMomentumFields(nfUEqn);
            Foam::volScalarField forAU("forAU", 1.0 / ofUEqn.A());
            Foam::volVectorField HbyA("HbyA", forAU * ofUEqn.H());
            HbyA.write();

            auto hostnfHbyA = nfHbyA.internalVector().copyToHost();
            write(nfHbyA.internalVector(), mesh, "nfHbyA" + execName);

            for (size_t celli = 0; celli < hostnfHbyA.size(); celli++)
            {
                REQUIRE(hostnfHbyA.view()[celli][0] == Catch::Approx(HbyA[celli][0]).margin(1e-14));
                REQUIRE(hostnfHbyA.view()[celli][1] == Catch::Approx(HbyA[celli][1]).margin(1e-14));
                REQUIRE(hostnfHbyA.view()[celli][2] == Catch::Approx(HbyA[celli][2]).margin(1e-14));
            }

            SECTION("constrainHbyA")
            {
                Foam::volVectorField ofConstrainHbyA(
                    "ofConstrainHbyA",
                    Foam::constrainHbyA(forAU * ofUEqn.H(), ofU, ofp)
                );
                nffvcc::constrainHbyA(nfHbyA, nfU, nfp);
                auto hostBCnfHbyA = nfHbyA.boundaryData().value().copyToHost();

                forAll(ofConstrainHbyA.boundaryField(), patchi)
                {
                    REQUIRE(
                        ofConstrainHbyA.boundaryField()[patchi].size()
                        == nfHbyA.boundaryData().nBoundaryFaces(patchi)
                    );
                    const Foam::fvPatchVectorField& ofConstrainHbyAPatch =
                        ofConstrainHbyA.boundaryField()[patchi];
                    auto [start, end] = nfHbyA.boundaryData().range(patchi);

                    forAll(ofConstrainHbyAPatch, bfacei)
                    {
                        REQUIRE(
                            hostBCnfHbyA.view()[start + bfacei][0]
                            == Catch::Approx(ofConstrainHbyAPatch[bfacei][0]).margin(1e-14)
                        );
                        REQUIRE(
                            hostBCnfHbyA.view()[start + bfacei][1]
                            == Catch::Approx(ofConstrainHbyAPatch[bfacei][1]).margin(1e-14)
                        );
                        REQUIRE(
                            hostBCnfHbyA.view()[start + bfacei][2]
                            == Catch::Approx(ofConstrainHbyAPatch[bfacei][2]).margin(1e-14)
                        );
                    }
                }
            }
        }

        SECTION("HbyA modified U")
        {
            ofU.primitiveFieldRef() *= 2.5;
            ofU.correctBoundaryConditions();
            nfU.internalVector() *= 2.5;
            nfU.correctBoundaryConditions();
            auto [nfrAU, nfHbyA] = nffvcc::discreteMomentumFields(nfUEqn);
            Foam::volScalarField forAU("forAU", 1.0 / ofUEqn.A());
            Foam::volVectorField HbyA("HbyA", forAU * ofUEqn.H());
            HbyA.write();

            auto hostnfHbyA = nfHbyA.internalVector().copyToHost();
            write(nfHbyA.internalVector(), mesh, "nfHbyA" + execName);

            for (size_t celli = 0; celli < hostnfHbyA.size(); celli++)
            {
                REQUIRE(hostnfHbyA.view()[celli][0] == Catch::Approx(HbyA[celli][0]).margin(1e-14));
                REQUIRE(hostnfHbyA.view()[celli][1] == Catch::Approx(HbyA[celli][1]).margin(1e-14));
                REQUIRE(hostnfHbyA.view()[celli][2] == Catch::Approx(HbyA[celli][2]).margin(1e-14));
            }
        }
    }
}
