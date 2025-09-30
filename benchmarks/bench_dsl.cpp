// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2023 NeoN authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "FoamAdapter/FoamAdapter.hpp"

#include "NeoN/NeoN.hpp"
#include "benchmarks/catch_main.hpp"
#include "test/catch2/executorGenerator.hpp"
#include "common.hpp"

namespace nnfvcc = NeoN::finiteVolume::cellCentred;
namespace nffvcc = FoamAdapter;
namespace dsl = NeoN::dsl;

#include "fvc.H"
#include "fvm.H"
#include "fvMatrices.H"

TEST_CASE("advection–diffusion-equation_scalar")
{
    Foam::Time& runTime = *timePtr;

    SECTION("OpenFOAM")
    {
        std::unique_ptr<Foam::fvMesh> meshPtr = FoamAdapter::createMesh(runTime);
        Foam::fvMesh& mesh = *meshPtr;

        auto ofT = randomScalarField(runTime, mesh, "T");
        Foam::surfaceScalarField ofPhi(
            Foam::IOobject("phi", "0", mesh, Foam::IOobject::NO_READ, Foam::IOobject::AUTO_WRITE),
            mesh,
            Foam::dimensionedScalar("phi", Foam::dimensionSet(0, 3, -1, 0, 0), 0.0)
        );
        forAll(ofPhi, facei)
        {
            ofPhi[facei] = facei;
        }

        Foam::surfaceScalarField ofGamma(
            Foam::IOobject("Gamma", "0", mesh, Foam::IOobject::NO_READ, Foam::IOobject::AUTO_WRITE),
            mesh,
            Foam::dimensionedScalar("phi", Foam::dimensionSet(0, 2, -1, 0, 0), 0.0)
        );


        SECTION("explicit")
        {

            BENCHMARK(std::string("OpenFOAM"))
            {
                Foam::fvScalarMatrix advectDiffEqn(
                    Foam::fvm::ddt(ofT) + Foam::fvc::div(ofPhi, ofT)
                    - Foam::fvc::laplacian(ofGamma, ofT)
                );
                return;
            };
        }

        SECTION("implicit")
        {
            BENCHMARK(std::string("OpenFOAM"))
            {
                Foam::fvScalarMatrix advectDiffEqn(
                    Foam::fvm::ddt(ofT) + Foam::fvm::div(ofPhi, ofT)
                    - Foam::fvm::laplacian(ofGamma, ofT)
                );
                return;
            };
        }
    }


    SECTION("NeoN")
    {
        auto [execName, exec] = GENERATE(allAvailableExecutor());

        NeoN::Database db;

        auto& vectorCollection = nnfvcc::VectorCollection::instance(db, "VectorCollection");

        std::unique_ptr<FoamAdapter::MeshAdapter> meshPtr = FoamAdapter::createMesh(exec, runTime);
        FoamAdapter::MeshAdapter& mesh = *meshPtr;
        const auto& nfMesh = mesh.nfMesh();

        auto ofT = randomScalarField(runTime, mesh, "T");
        nnfvcc::VolumeField<NeoN::scalar>& nfT =
            vectorCollection.registerVector<nnfvcc::VolumeField<NeoN::scalar>>(
                FoamAdapter::CreateFromFoamField<Foam::volScalarField> {
                    .exec = exec,
                    .nfMesh = nfMesh,
                    .foamField = ofT,
                    .name = "nfT"
                }
            );
        nfT.correctBoundaryConditions();

        auto& nfOldT = fvcc::oldTime(nfT);
        nfOldT.internalVector() = nfT.internalVector();

        Foam::surfaceScalarField ofPhi(
            Foam::IOobject("phi", "0", mesh, Foam::IOobject::NO_READ, Foam::IOobject::AUTO_WRITE),
            mesh,
            Foam::dimensionedScalar("phi", Foam::dimensionSet(0, 3, -1, 0, 0), 0.0)
        );
        for (auto facei = 0; facei < ofPhi.size(); facei++)
        {
            ofPhi[facei] = facei;
        }

        auto nfPhi = FoamAdapter::constructSurfaceField(exec, nfMesh, ofPhi);

        Foam::surfaceScalarField ofGamma(
            Foam::IOobject("Gamma", "0", mesh, Foam::IOobject::NO_READ, Foam::IOobject::AUTO_WRITE),
            mesh,
            Foam::dimensionedScalar("Gamma", Foam::dimensionSet(0, 2, -1, 0, 0), 1.0)
        );

        auto nfGamma = FoamAdapter::constructSurfaceField(exec, nfMesh, ofGamma);

        NeoN::Dictionary fvSolutionExpDict;

        SECTION(std::string("explicit"))
        {
            NeoN::Dictionary fvSchemesExpDict;
            fvSchemesExpDict.insert(
                std::string("ddtSchemes"),
                NeoN::Dictionary(
                    {{std::string("type"), NeoN::TokenList({std::string("forwardEuler")})}}
                )
            );
            fvSchemesExpDict.insert(
                std::string("divSchemes"),
                NeoN::Dictionary(
                    {{std::string("div(phi,nfT)"),
                      NeoN::TokenList({std::string("Gauss"), std::string("upwind")})}}
                )
            );
            fvSchemesExpDict.insert(
                std::string("laplacianSchemes"),
                NeoN::Dictionary(
                    {{std::string("laplacian(Gamma,nfT)"),
                      NeoN::TokenList(
                          {std::string("Gauss"), std::string("linear"), std::string("uncorrected")}
                      )}}
                )
            );


            BENCHMARK(std::string(execName))
            {
                // Momentum predictor
                nffvcc::Expression<NeoN::scalar> advectDiffEqn(
                    dsl::imp::ddt(nfT) + dsl::exp::div(nfPhi, nfT)
                        - dsl::exp::laplacian(nfGamma, nfT),
                    nfT,
                    fvSchemesExpDict,
                    fvSolutionExpDict
                );
                advectDiffEqn.assemble(1.0, 1e-3);
                return;
            };
        }

        SECTION(std::string("implicit"))
        {
            NeoN::Dictionary fvSchemesImpDict;
            fvSchemesImpDict.insert(
                std::string("ddtSchemes"),
                NeoN::Dictionary(
                    {{std::string("type"), NeoN::TokenList({std::string("backwardEuler")})}}
                )
            );
            fvSchemesImpDict.insert(
                std::string("divSchemes"),
                NeoN::Dictionary(
                    {{std::string("div(phi,nfT)"),
                      NeoN::TokenList({std::string("Gauss"), std::string("upwind")})}}
                )
            );
            fvSchemesImpDict.insert(
                std::string("laplacianSchemes"),
                NeoN::Dictionary(
                    {{std::string("laplacian(Gamma,nfT)"),
                      NeoN::TokenList(
                          {std::string("Gauss"), std::string("linear"), std::string("uncorrected")}
                      )}}
                )
            );

            BENCHMARK(std::string(execName))
            {
                // Momentum predictor
                nffvcc::Expression<NeoN::scalar> advectDiffEqn(
                    dsl::imp::ddt(nfT) + dsl::imp::div(nfPhi, nfT)
                        - dsl::imp::laplacian(nfGamma, nfT),
                    nfT,
                    fvSchemesImpDict,
                    fvSolutionExpDict
                );
                advectDiffEqn.assemble(1.0, 1e-3);
                return;
            };
        }
    }
}
