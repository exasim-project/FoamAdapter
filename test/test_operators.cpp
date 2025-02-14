// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenGrad.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenDiv.hpp"

#define namespaceFoam // Suppress <using namespace Foam;>
#include "gaussConvectionScheme.H"
#include "NeoFOAM/core/input.hpp"
#include "NeoFOAM/dsl/explicit.hpp"

#include "common.hpp"

namespace fvcc = NeoFOAM::finiteVolume::cellCentred;
namespace dsl = NeoFOAM::dsl;

extern Foam::Time* timePtr;    // A single time object
extern Foam::argList* argsPtr; // Some forks want argList access at createMesh.H
extern Foam::fvMesh* meshPtr;  // A single mesh object


TEST_CASE("Interpolation")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    NeoFOAM::Executor exec = GENERATE(
        NeoFOAM::Executor(NeoFOAM::SerialExecutor {}),
        NeoFOAM::Executor(NeoFOAM::CPUExecutor {}),
        NeoFOAM::Executor(NeoFOAM::GPUExecutor {})
    );

    std::string execName = std::visit([](auto e) { return e.name(); }, exec);

    auto meshPtr = Foam::createMesh(exec, runTime);
    Foam::MeshAdapter& mesh = *meshPtr;
    auto nfMesh = mesh.nfMesh();

    SECTION("scalar_" + execName)
    {
        Foam::IStringStream is("linear");
        Foam::tmp<Foam::surfaceInterpolationScheme<Foam::scalar>> foamInterPol =
            Foam::surfaceInterpolationScheme<Foam::scalar>::New(mesh, is);

        auto ofT = randomScalarField(runTime, mesh);
        auto nfT = constructFrom(exec, nfMesh, ofT);
        nfT.correctBoundaryConditions();
        REQUIRE(nfT == ofT);

        Foam::surfaceScalarField ofSurfT(foamInterPol->interpolate(ofT));
        auto nfSurfT = constructSurfaceField(exec, nfMesh, ofSurfT);

        SECTION("linear")
        {
            NeoFOAM::Input interpolationScheme = NeoFOAM::TokenList({std::string("linear")});
            auto linearKernel =
                fvcc::SurfaceInterpolationFactory::create(exec, nfMesh, interpolationScheme);
            fvcc::SurfaceInterpolation interp(exec, nfMesh, std::move(linearKernel));
            // TODO since it is constructed from ofField it is trivial
            // we should reset the field first
            interp.interpolate(nfT, nfSurfT);
            nfSurfT.correctBoundaryConditions();
            compare(nfSurfT, ofSurfT, ApproxScalar(1e-15), false);
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

    std::string execName = std::visit([](auto e) { return e.name(); }, exec);

    std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, runTime);
    Foam::MeshAdapter& mesh = *meshPtr;
    NeoFOAM::UnstructuredMesh& nfMesh = mesh.nfMesh();

    SECTION("GaussGrad on " + execName)
    {
        // linear interpolation hardcoded for now
        Foam::IStringStream is("linear");

        auto ofT = randomScalarField(runTime, mesh);
        auto nfT = constructFrom(exec, nfMesh, ofT);
        nfT.correctBoundaryConditions();
        REQUIRE(nfT == ofT);

        Foam::fv::gaussGrad<Foam::scalar> foamGradScalar(mesh, is);
        Foam::volVectorField ofGradT("ofGradT", foamGradScalar.calcGrad(ofT, "test"));

        auto nfGradT = constructFrom(exec, nfMesh, ofGradT);
        NeoFOAM::fill(nfGradT.internalField(), NeoFOAM::Vector(0.0, 0.0, 0.0));
        NeoFOAM::fill(nfGradT.boundaryField().value(), NeoFOAM::Vector(0.0, 0.0, 0.0));
        fvcc::GaussGreenGrad(exec, nfMesh).grad(nfT, nfGradT);
        nfGradT.correctBoundaryConditions();

        compare(nfGradT, ofGradT, ApproxVector(1e-12), false);
    }
}


TEST_CASE("DivOperator")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    NeoFOAM::Executor exec = GENERATE(
        NeoFOAM::Executor(NeoFOAM::SerialExecutor {}),
        NeoFOAM::Executor(NeoFOAM::CPUExecutor {}),
        NeoFOAM::Executor(NeoFOAM::GPUExecutor {})
    );
    std::string execName = std::visit([](auto e) { return e.name(); }, exec);

    auto meshPtr = Foam::createMesh(exec, runTime);
    Foam::MeshAdapter& mesh = *meshPtr;
    auto nfMesh = mesh.nfMesh();

    SECTION("gaussDiv_scalar_" + execName)
    {
        // linear interpolation hardcoded for now
        Foam::IStringStream is("linear");

        auto ofT = randomScalarField(runTime, mesh);
        auto nfT = constructFrom(exec, nfMesh, ofT);
        nfT.correctBoundaryConditions();
        REQUIRE(nfT == ofT);

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

        Foam::fv::gaussConvectionScheme<Foam::scalar> foamDivScalar(mesh, ofPhi, is);
        Foam::volScalarField ofDivT("ofDivT", foamDivScalar.fvcDiv(ofPhi, ofT));

        auto nfPhi = constructSurfaceField(exec, nfMesh, ofPhi);

        SECTION("Gauss-Green")
        {
            auto nfDivT = constructFrom(exec, nfMesh, ofDivT);
            NeoFOAM::TokenList scheme({std::string("linear")});
            // Reset
            NeoFOAM::fill(nfDivT.internalField(), 0.0);
            NeoFOAM::fill(nfDivT.boundaryField().value(), 0.0);
            fvcc::GaussGreenDiv(exec, nfMesh, scheme).div(nfDivT, nfPhi, nfT);
            nfDivT.correctBoundaryConditions();

            compare(nfDivT, ofDivT, ApproxScalar(1e-15), false);
        }

        SECTION("compute div from dsl::exp")
        {
            NeoFOAM::TokenList scheme =
                NeoFOAM::TokenList({std::string("Gauss"), std::string("linear")});

            auto nfDivT = constructFrom(exec, nfMesh, ofDivT);
            NeoFOAM::fill(nfDivT.internalField(), 0.0);
            NeoFOAM::fill(nfDivT.boundaryField().value(), 0.0);
            dsl::SpatialOperator divOp = dsl::exp::div(nfPhi, nfT);
            divOp.build(scheme);
            divOp.explicitOperation(nfDivT.internalField());

            nfDivT.correctBoundaryConditions();

            compare(nfT, ofT, ApproxScalar(1e-15), false);

            compare(nfDivT, ofDivT, ApproxScalar(1e-15), false);
        }
    }
}
