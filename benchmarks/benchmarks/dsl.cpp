// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2023 NeoN authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "NeoN/NeoN.hpp"
#include "benchmarks/catch_main.hpp"
#include "NeoN/test/catch2/executorGenerator.hpp"

namespace fvcc = NeoN::finiteVolume::cellCentred;
namespace dsl = NeoN::dsl;

#include "fvc.H"
#include "fvm.H"
#include "fvMatrices.H"
#include "NeoN/dsl/explicit.hpp"
#include "NeoN/include/NeoN/linearAlgebra/linearSystem.hpp"

template<typename FieldType, typename RandomFunc>
FieldType createRandomField(
    const Foam::Time& runTime,
    const Foam::fvMesh& mesh,
    Foam::word name,
    RandomFunc rand
)
{
    FieldType t(
        Foam::IOobject(name, "0", mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE),
        mesh
    );

    forAll(t, celli)
    {
        t[celli] = rand();
    }

    t.correctBoundaryConditions();
    return t;
}


/* function to create a volScalarField filled with random values for test purposes */
auto randomScalarField(const Foam::Time& runTime, const Foam::fvMesh& mesh, Foam::word name)
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);
    return createRandomField<Foam::volScalarField>(runTime, mesh, name, [&]() { return dis(gen); });
}

TEST_CASE("advection–diffusion-equation_scalar")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    runTime.setDeltaT(1e-3);
    runTime++;

    SECTION("OpenFOAM")
    {
        std::unique_ptr<Foam::fvMesh> meshPtr = Foam::createMesh(runTime);
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

        fvcc::VectorCollection& vectorCollection =
            fvcc::VectorCollection::instance(db, "VectorCollection");

        std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, runTime);
        Foam::MeshAdapter& mesh = *meshPtr;
        const auto& nfMesh = mesh.nfMesh();
        // linear interpolation hardcoded for now


        auto ofT = randomScalarField(runTime, mesh, "T");
        fvcc::VolumeField<NeoN::scalar>& nfT =
            vectorCollection.registerVector<fvcc::VolumeField<NeoN::scalar>>(
                Foam::CreateFromFoamField<Foam::volScalarField> {
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
        forAll(ofPhi, facei)
        {
            ofPhi[facei] = facei;
        }

        auto nfPhi = constructSurfaceField(exec, nfMesh, ofPhi);

        Foam::surfaceScalarField ofGamma(
            Foam::IOobject("Gamma", "0", mesh, Foam::IOobject::NO_READ, Foam::IOobject::AUTO_WRITE),
            mesh,
            Foam::dimensionedScalar("Gamma", Foam::dimensionSet(0, 2, -1, 0, 0), 1.0)
        );

        auto nfGamma = constructSurfaceField(exec, nfMesh, ofGamma);


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
                fvcc::Expression<NeoN::scalar> advectDiffEqn(
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
                fvcc::Expression<NeoN::scalar> advectDiffEqn(
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
