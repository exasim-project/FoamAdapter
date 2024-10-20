// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main
#include <random>
#include <span>

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/catch_approx.hpp>
#include "catch2/common.hpp"

#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenGrad.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenDiv.hpp"

#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/writers/writers.hpp"
#include "FoamAdapter/comparison/fieldComparison.hpp"
#include "FoamAdapter/setup/setup.hpp"

#define namespaceFoam // Suppress <using namespace Foam;>
#include "gaussConvectionScheme.H"

namespace fvcc = NeoFOAM::finiteVolume::cellCentred;

extern Foam::Time* timePtr;    // A single time object
extern Foam::argList* argsPtr; // Some forks want argList access at createMesh.H
extern Foam::fvMesh* meshPtr;  // A single mesh object


Foam::volScalarField createRandomField(const Foam::Time& runTime, const Foam::fvccNeoMesh& mesh)
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);

    Foam::volScalarField t(
        Foam::IOobject(
            "T", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
        ),
        mesh
    );

    forAll(t, celli)
    {
        t[celli] = dis(gen);
    }

    t.correctBoundaryConditions();
    return t;
}

template<typename NFFIELD, typename OFFIELD, typename Compare>
void compare(NFFIELD& a, OFFIELD& b, Compare comp)
{
    auto aHost = a.internalField().copyToHost();
    auto bSpan = std::span(b.primitiveFieldRef().data(), b.size());
    // nf a span might be shorter than bSpan for surface fields
    REQUIRE_THAT(aHost.span({0, bSpan.size()}), Catch::Matchers::RangeEquals(bSpan, comp));
}

TEST_CASE("Interpolation")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    NeoFOAM::Executor exec = GENERATE(
        NeoFOAM::Executor(NeoFOAM::SerialExecutor {}),
        NeoFOAM::Executor(NeoFOAM::CPUExecutor {}),
        NeoFOAM::Executor(NeoFOAM::GPUExecutor {})
    );

    std::string execName = std::visit([](auto e) { return e.print(); }, exec);

    auto meshPtr = Foam::createMesh(exec, runTime);
    Foam::fvccNeoMesh& mesh = *meshPtr;
    auto nfMesh = mesh.nfMesh();

    SECTION("scalar_" + execName)
    {
        Foam::IStringStream is("linear");
        Foam::tmp<Foam::surfaceInterpolationScheme<Foam::scalar>> foamInterPol =
            Foam::surfaceInterpolationScheme<Foam::scalar>::New(mesh, is);

        auto ofT = createRandomField(runTime, mesh);
        auto nfT = constructFrom(exec, nfMesh, ofT);
        nfT.correctBoundaryConditions();
        REQUIRE(nfT == ofT);

        Foam::surfaceScalarField ofSurfT(foamInterPol->interpolate(ofT));
        auto nfSurfT = constructSurfaceField(exec, nfMesh, ofSurfT);

        SECTION("linear")
        {
            auto linearKernel = fvcc::SurfaceInterpolationFactory::create("linear", exec, nfMesh);
            fvcc::SurfaceInterpolation interp(exec, nfMesh, std::move(linearKernel));
            interp.interpolate(nfT, nfSurfT);
            compare(nfSurfT, ofSurfT, ApproxScalar(1e-15));
        }
    }
}

TEST_CASE("GradOperator")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    NeoFOAM::Executor exec = GENERATE(
        NeoFOAM::Executor(NeoFOAM::CPUExecutor {}),
        NeoFOAM::Executor(NeoFOAM::SerialExecutor {}),
        NeoFOAM::Executor(NeoFOAM::GPUExecutor {})
    );

    std::string execName = std::visit([](auto e) { return e.print(); }, exec);

    std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec, runTime);
    Foam::fvccNeoMesh& mesh = *meshPtr;
    NeoFOAM::UnstructuredMesh& nfMesh = mesh.nfMesh();

    SECTION("GaussGrad on " + execName)
    {
        // linear interpolation hardcoded for now
        Foam::IStringStream is("linear");

        auto ofT = createRandomField(runTime, mesh);
        auto nfT = constructFrom(exec, nfMesh, ofT);
        nfT.correctBoundaryConditions();
        REQUIRE(nfT == ofT);

        Foam::fv::gaussGrad<Foam::scalar> foamGradScalar(mesh, is);
        Foam::volVectorField ofGradT("ofGradT", foamGradScalar.calcGrad(ofT, "test"));

        auto nfGradT = constructFrom(exec, nfMesh, ofGradT);
        NeoFOAM::fill(nfGradT.internalField(), NeoFOAM::Vector(0.0, 0.0, 0.0));
        NeoFOAM::fill(nfGradT.boundaryField().value(), NeoFOAM::Vector(0.0, 0.0, 0.0));
        fvcc::GaussGreenGrad(exec, nfMesh).grad(nfT, nfGradT);

        compare(nfGradT, ofGradT, ApproxVector(1e-12));
    }
}


TEST_CASE("DivOperator")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    NeoFOAM::Executor exec = GENERATE(
        NeoFOAM::Executor(NeoFOAM::CPUExecutor {}),
        NeoFOAM::Executor(NeoFOAM::SerialExecutor {}),
        NeoFOAM::Executor(NeoFOAM::GPUExecutor {})
    );
    std::string execName = std::visit([](auto e) { return e.print(); }, exec);

    auto meshPtr = Foam::createMesh(exec, runTime);
    Foam::fvccNeoMesh& mesh = *meshPtr;
    auto nfMesh = mesh.nfMesh();

    SECTION("gaussDiv_scalar_" + execName)
    {
        // linear interpolation hardcoded for now
        Foam::IStringStream is("linear");

        auto ofT = createRandomField(runTime, mesh);
        auto nfT = constructFrom(exec, nfMesh, ofT);
        nfT.correctBoundaryConditions();
        REQUIRE(nfT == ofT);

        Foam::surfaceScalarField ofPhi(
            Foam::IOobject(
                "phi", runTime.timeName(), mesh, Foam::IOobject::NO_READ, Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("phi", Foam::dimless, 0.0)
        );
        forAll(ofPhi, facei)
        {
            ofPhi[facei] = facei;
        }

        Foam::fv::gaussConvectionScheme<Foam::scalar> foamDivScalar(mesh, ofPhi, is);
        Foam::volScalarField ofDivT("ofDivT", foamDivScalar.fvcDiv(ofPhi, ofT));

        auto nfPhi = constructSurfaceField(exec, nfMesh, ofPhi);

        SECTION("Computes correct divT")
        {
            auto nfDivT = constructFrom(exec, nfMesh, ofDivT);
            NeoFOAM::fill(nfDivT.internalField(), 0.0);
            NeoFOAM::fill(nfDivT.boundaryField().value(), 0.0);
            fvcc::GaussGreenDiv(
                exec,
                nfMesh,
                fvcc::SurfaceInterpolation(
                    exec, nfMesh, std::make_unique<fvcc::Linear>(exec, nfMesh)
                )
            )
                .div(nfDivT, nfPhi, nfT);

            compare(nfDivT, ofDivT, ApproxScalar(1e-15));
        }
    }
}
