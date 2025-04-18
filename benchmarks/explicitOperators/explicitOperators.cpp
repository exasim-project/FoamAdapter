// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2023 NeoN authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "NeoN/NeoN.hpp"
#include "benchmarks/catch_main.hpp"
#include "NeoN/test/catch2/executorGenerator.hpp"

namespace fvcc = NeoN::finiteVolume::cellCentred;
namespace dsl = NeoN::dsl;

#include "fv.H"
#include "fvc.H"
#include "gaussGrad.H"
#include "gaussConvectionScheme.H"
#include "gaussLaplacianScheme.H"
#include "NeoN/core/input.hpp"
#include "NeoN/dsl/explicit.hpp"

template<typename FieldType, typename RandomFunc>
FieldType createRandomField(
    const Foam::Time& runTime,
    const Foam::fvMesh& mesh,
    Foam::word name,
    RandomFunc rand
)
{
    FieldType t(
        Foam::IOobject(
            name,
            runTime.timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::AUTO_WRITE
        ),
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

TEST_CASE("DivOperator")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    SECTION("OpenFOAM")
    {
        std::unique_ptr<Foam::fvMesh> meshPtr = Foam::createMesh(runTime);
        Foam::fvMesh& mesh = *meshPtr;

        auto ofT = randomScalarField(runTime, mesh, "T");
        Foam::surfaceScalarField ofPhi(
            Foam::IOobject(
                "phi",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("phi", Foam::dimless, 0.0)
        );
        forAll(ofPhi, facei)
        {
            ofPhi[facei] = facei;
        }
        SECTION("with Allocation")
        {
            BENCHMARK(std::string("OpenFOAM"))
            {
                Foam::IStringStream is("linear");
                Foam::fv::gaussConvectionScheme<Foam::scalar> foamDivScalar(mesh, ofPhi, is);
                Foam::volScalarField ofDivT("ofDivT", foamDivScalar.fvcDiv(ofPhi, ofT));
                return;
            };
        }
    }


    SECTION("NeoN")
    {
        auto [execName, exec] = GENERATE(allAvailableExecutor());

        std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, runTime);
        Foam::MeshAdapter& mesh = *meshPtr;
        const auto& nfMesh = mesh.nfMesh();
        // linear interpolation hardcoded for now


        auto ofT = randomScalarField(runTime, mesh, "T");
        auto nfT = constructFrom(exec, nfMesh, ofT);
        nfT.correctBoundaryConditions();

        Foam::surfaceScalarField ofPhi(
            Foam::IOobject(
                "phi",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("phi", Foam::dimless, 0.0)
        );
        forAll(ofPhi, facei)
        {
            ofPhi[facei] = facei;
        }

        auto nfPhi = constructSurfaceField(exec, nfMesh, ofPhi);

        SECTION("with Allocation")
        {
            NeoN::TokenList scheme({std::string("linear")});

            BENCHMARK(std::string(execName))
            {
                fvcc::VolumeField<NeoN::scalar> divPhiT =
                    fvcc::GaussGreenDiv<NeoN::scalar>(exec, nfMesh, scheme)
                        .div(nfPhi, nfT, dsl::Coeff(1.0));
                return;
            };
        }

        SECTION("No allocation")
        {
            auto nfDivT = constructFrom(exec, nfMesh, ofT);
            NeoN::TokenList scheme({std::string("linear")});

            BENCHMARK(std::string(execName))
            {
                NeoN::fill(nfDivT.internalField(), 0.0);
                NeoN::fill(nfDivT.boundaryField().value(), 0.0);
                fvcc::GaussGreenDiv<NeoN::scalar>(exec, nfMesh, scheme)
                    .div(nfDivT, nfPhi, nfT, dsl::Coeff(1.0));
                Kokkos::fence();
                return;
            };
        }

        // TODO: dsl
    }
}


TEST_CASE("LaplacianOperator")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    SECTION("OpenFOAM")
    {
        std::unique_ptr<Foam::fvMesh> meshPtr = Foam::createMesh(runTime);
        Foam::fvMesh& mesh = *meshPtr;

        auto ofT = randomScalarField(runTime, mesh, "T");
        Foam::surfaceScalarField ofGamma(
            Foam::IOobject(
                "Gamma",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("Gamma", Foam::dimless, 1.0)
        );

        SECTION("with Allocation")
        {
            BENCHMARK(std::string("OpenFOAM"))
            {
                Foam::IStringStream is("linear uncorrected");
                Foam::fv::gaussLaplacianScheme<Foam::scalar, Foam::scalar> foamLapScalar(mesh, is);
                Foam::volScalarField ofLapT("ofLapT", foamLapScalar.fvcLaplacian(ofGamma, ofT));
                return;
            };
        }
    }


    SECTION("NeoN")
    {
        auto [execName, exec] = GENERATE(allAvailableExecutor());

        std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, runTime);
        Foam::MeshAdapter& mesh = *meshPtr;
        const auto& nfMesh = mesh.nfMesh();
        // linear interpolation hardcoded for now


        auto ofT = randomScalarField(runTime, mesh, "T");
        auto nfT = constructFrom(exec, nfMesh, ofT);
        nfT.correctBoundaryConditions();

        Foam::surfaceScalarField ofGamma(
            Foam::IOobject(
                "Gamma",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("Gamma", Foam::dimless, 1.0)
        );


        auto nfGamma = constructSurfaceField(exec, nfMesh, ofGamma);

        SECTION("with Allocation")
        {
            NeoN::TokenList scheme({std::string("linear")});

            BENCHMARK(std::string(execName))
            {
                fvcc::VolumeField<NeoN::scalar> lapT =
                    fvcc::GaussGreenLaplacian<NeoN::scalar>(exec, nfMesh, scheme)
                        .laplacian(nfGamma, nfT, dsl::Coeff(1.0));
                return;
            };
        }

        SECTION("No allocation")
        {
            auto nfLapT = constructFrom(exec, nfMesh, ofT);
            NeoN::TokenList scheme({std::string("linear"), std::string("uncorrected")});

            BENCHMARK(std::string(execName))
            {
                NeoN::fill(nfLapT.internalField(), 0.0);
                NeoN::fill(nfLapT.boundaryField().value(), 0.0);
                fvcc::GaussGreenLaplacian<NeoN::scalar>(exec, nfMesh, scheme)
                    .laplacian(nfLapT, nfGamma, nfT, dsl::Coeff(1.0));
                Kokkos::fence();
                return;
            };
        }

        // TODO: dsl
    }
}

TEST_CASE("GradOperator")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    SECTION("OpenFOAM")
    {
        std::unique_ptr<Foam::fvMesh> meshPtr = Foam::createMesh(runTime);
        Foam::fvMesh& mesh = *meshPtr;

        auto ofT = randomScalarField(runTime, mesh, "T");

        SECTION("with Allocation")
        {
            BENCHMARK(std::string("OpenFOAM"))
            {
                Foam::IStringStream is("linear");
                Foam::fv::gaussGrad<Foam::scalar> foamGrad(mesh, is);
                Foam::volVectorField ofGradT("ofGradT", foamGrad.calcGrad(ofT, "ofGradT"));
                return;
            };
        }
    }


    SECTION("NeoN")
    {
        auto [execName, exec] = GENERATE(allAvailableExecutor());

        std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, runTime);
        Foam::MeshAdapter& mesh = *meshPtr;
        const auto& nfMesh = mesh.nfMesh();
        // linear interpolation hardcoded for now


        auto ofT = randomScalarField(runTime, mesh, "T");
        auto nfT = constructFrom(exec, nfMesh, ofT);

        SECTION("with Allocation")
        {
            NeoN::TokenList scheme({std::string("linear")});

            BENCHMARK(std::string(execName))
            {
                fvcc::VolumeField<NeoN::Vec3> nfGradT =
                    fvcc::GaussGreenGrad(exec, nfMesh).grad(nfT);
                return;
            };
        }

        SECTION("No allocation")
        {
            fvcc::VolumeField<NeoN::Vec3> nfGradT = fvcc::GaussGreenGrad(exec, nfMesh).grad(nfT);

            BENCHMARK(std::string(execName))
            {
                NeoN::fill(nfGradT.internalField(), NeoN::Vec3(0, 0, 0));
                NeoN::fill(nfGradT.boundaryField().value(), NeoN::Vec3(0, 0, 0));
                fvcc::GaussGreenGrad(exec, nfMesh).grad(nfT, nfGradT);
                Kokkos::fence();
                return;
            };
        }

        // TODO: dsl
    }
}


TEST_CASE("FaceInterpolation")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    SECTION("OpenFOAM")
    {
        std::unique_ptr<Foam::fvMesh> meshPtr = Foam::createMesh(runTime);
        Foam::fvMesh& mesh = *meshPtr;

        auto ofT = randomScalarField(runTime, mesh, "T");
        Foam::surfaceScalarField ofPhi(
            Foam::IOobject(
                "phi",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("phi", Foam::dimless, 0.0)
        );
        forAll(ofPhi, facei)
        {
            ofPhi[facei] = facei;
        }

        SECTION("with Allocation")
        {
            BENCHMARK(std::string("OpenFOAM"))
            {
                Foam::IStringStream is("linear");
                auto Tf = Foam::fvc::interpolate(ofT, ofPhi, is);
                return;
            };
        }
    }


    SECTION("NeoN")
    {
        auto [execName, exec] = GENERATE(allAvailableExecutor());

        std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, runTime);
        Foam::MeshAdapter& mesh = *meshPtr;
        const auto& nfMesh = mesh.nfMesh();
        // linear interpolation hardcoded for now


        auto ofT = randomScalarField(runTime, mesh, "T");
        auto nfT = constructFrom(exec, nfMesh, ofT);
        nfT.correctBoundaryConditions();

        Foam::surfaceScalarField ofPhi(
            Foam::IOobject(
                "phi",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("phi", Foam::dimless, 0.0)
        );
        forAll(ofPhi, facei)
        {
            ofPhi[facei] = facei;
        }

        auto nfPhi = constructSurfaceField(exec, nfMesh, ofPhi);

        SECTION("with Allocation")
        {
            NeoN::TokenList scheme({std::string("linear")});

            BENCHMARK(std::string(execName))
            {
                // fvcc::SurfaceField<NeoN::scalar> Tf =
                fvcc::SurfaceInterpolation<NeoN::scalar>(exec, nfMesh, scheme)
                    .interpolate(nfPhi, nfT);
                return;
            };
        }
        SECTION("No allocation")
        {
            NeoN::TokenList scheme({std::string("linear")});

            fvcc::SurfaceField<NeoN::scalar> nfTf(nfPhi);
            nfTf.name = "Tf";

            BENCHMARK(std::string(execName))
            {
                NeoN::fill(nfTf.internalField(), 0.0);
                NeoN::fill(nfTf.boundaryField().value(), 0.0);
                fvcc::SurfaceInterpolation<NeoN::scalar>(exec, nfMesh, scheme)
                    .interpolate(nfPhi, nfT, nfTf);
                Kokkos::fence();
                return;
            };
        }

        // TODO: dsl
    }
}


TEST_CASE("FaceNormalGradient")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    SECTION("OpenFOAM")
    {
        std::unique_ptr<Foam::fvMesh> meshPtr = Foam::createMesh(runTime);
        Foam::fvMesh& mesh = *meshPtr;

        auto ofT = randomScalarField(runTime, mesh, "T");


        SECTION("with Allocation")
        {
            BENCHMARK(std::string("OpenFOAM"))
            {
                Foam::IStringStream is("uncorrected");
                auto snGradScheme = Foam::fv::snGradScheme<Foam::scalar>::New(mesh, is);
                snGradScheme->snGrad(ofT);
                return;
            };
        }
    }


    SECTION("NeoN")
    {
        auto [execName, exec] = GENERATE(allAvailableExecutor());

        std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, runTime);
        Foam::MeshAdapter& mesh = *meshPtr;
        const auto& nfMesh = mesh.nfMesh();
        // linear interpolation hardcoded for now


        auto ofT = randomScalarField(runTime, mesh, "T");
        auto nfT = constructFrom(exec, nfMesh, ofT);
        nfT.correctBoundaryConditions();


        SECTION("with Allocation")
        {
            NeoN::TokenList scheme({std::string("uncorrected")});

            BENCHMARK(std::string(execName))
            {
                // fvcc::SurfaceField<NeoN::scalar> Tf =
                fvcc::FaceNormalGradient<NeoN::scalar>(exec, nfMesh, scheme).faceNormalGrad(nfT);
                return;
            };
        }
        SECTION("No allocation")
        {
            NeoN::TokenList scheme({std::string("uncorrected")});

            std::string nameFaceGrad = "faceGrad_" + nfT.name;
            fvcc::SurfaceField<NeoN::scalar> faceGradT(
                exec,
                nameFaceGrad,
                nfMesh,
                fvcc::createCalculatedBCs<fvcc::SurfaceBoundary<NeoN::scalar>>(nfMesh)
            );

            BENCHMARK(std::string(execName))
            {
                NeoN::fill(faceGradT.internalField(), 0.0);
                NeoN::fill(faceGradT.boundaryField().value(), 0.0);
                fvcc::FaceNormalGradient<NeoN::scalar>(exec, nfMesh, scheme)
                    .faceNormalGrad(nfT, faceGradT);
                Kokkos::fence();
                return;
            };
        }

        // TODO: dsl
    }
}
