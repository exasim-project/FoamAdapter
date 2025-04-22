// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "common.hpp"

#include "fv.H"
#include "gaussConvectionScheme.H"


namespace fvcc = NeoN::finiteVolume::cellCentred;
namespace dsl = NeoN::dsl;


extern Foam::Time* timePtr;    // A single time object
extern Foam::argList* argsPtr; // Some forks want argList access at createMesh.H
extern Foam::fvMesh* meshPtr;  // A single mesh object

NeoN::TokenList interpolationScheme;

TEST_CASE("Interpolation")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    auto [execName, exec] = GENERATE(allAvailableExecutor());

    auto meshPtr = Foam::createMesh(exec, runTime);
    Foam::MeshAdapter& mesh = *meshPtr;
    auto nfMesh = mesh.nfMesh();

    auto ofT = randomScalarField(runTime, mesh, "T");
    auto nfT = constructFrom(exec, nfMesh, ofT);

    SECTION("Linear SurfaceInterpolation[scalar] on " + execName)
    {
        Foam::IStringStream is("linear");
        Foam::tmp<Foam::surfaceInterpolationScheme<Foam::scalar>> foamInterPol =
            Foam::surfaceInterpolationScheme<Foam::scalar>::New(mesh, is);
        Foam::surfaceScalarField ofSurfT(foamInterPol->interpolate(ofT));
        auto nfSurfT = constructSurfaceField(exec, nfMesh, ofSurfT);

        // Foam::surfaceInterpolationScheme<Foam::scalar>::New(mesh, is);
        // NeoN::TokenList interpolationScheme;
        interpolationScheme.insert(std::string("linear"));
        auto linearKernel = fvcc::SurfaceInterpolationFactory<NeoN::scalar>::create(
            exec,
            nfMesh,
            interpolationScheme
        );
        auto op = fvcc::SurfaceInterpolation(exec, nfMesh, std::move(linearKernel));
        // TODO since it is constructed from ofField it is trivial
        // we should reset the field first
        op.interpolate(nfT, nfSurfT);
        nfSurfT.correctBoundaryConditions();
        compare(nfSurfT, ofSurfT, ApproxScalar(1e-15), false);
    }

    SECTION("GaussGreenGrad[scalar] on " + execName)
    {
        // linear interpolation hardcoded for now
        Foam::IStringStream is("linear");

        Foam::fv::gaussGrad<Foam::scalar> foamGradScalar(mesh, is);
        Foam::volVectorField ofGradT("ofGradT", foamGradScalar.calcGrad(ofT, "test"));

        auto nfGradT = constructFrom(exec, nfMesh, ofGradT);
        // reset values
        NeoN::fill(nfGradT.internalVector(), NeoN::Vec3(0.0, 0.0, 0.0));
        NeoN::fill(nfGradT.boundaryData().value(), NeoN::Vec3(0.0, 0.0, 0.0));
        fvcc::GaussGreenGrad(exec, nfMesh).grad(nfT, nfGradT);
        nfGradT.correctBoundaryConditions();

        compare(nfGradT, ofGradT, ApproxVector(1e-12), false);
    }

    SECTION("Gauss-Green")
    {
        Foam::IStringStream is("linear");
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
        auto nfDivT = constructFrom(exec, nfMesh, ofDivT);
        auto nfPhi = constructSurfaceField(exec, nfMesh, ofPhi);
        // Reset
        NeoN::fill(nfDivT.internalVector(), 0.0);
        NeoN::fill(nfDivT.boundaryData().value(), 0.0);
        fvcc::GaussGreenDiv<NeoN::scalar>(exec, nfMesh, interpolationScheme)
            .div(nfDivT, nfPhi, nfT, dsl::Coeff(1.0));
        nfDivT.correctBoundaryConditions();

        compare(nfDivT, ofDivT, ApproxScalar(1e-15), false);
    }

    SECTION("compute div from dsl::exp")
    {
        Foam::IStringStream is("linear");
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
        NeoN::TokenList scheme = NeoN::TokenList({std::string("Gauss"), std::string("linear")});

        auto nfDivT = constructFrom(exec, nfMesh, ofDivT);
        auto nfPhi = constructSurfaceField(exec, nfMesh, ofPhi);
        NeoN::fill(nfDivT.internalVector(), 0.0);
        NeoN::fill(nfDivT.boundaryData().value(), 0.0);
        dsl::SpatialOperator divOp = dsl::exp::div(nfPhi, nfT);
        divOp.read(scheme);
        divOp.explicitOperation(nfDivT.internalVector());

        nfDivT.correctBoundaryConditions();

        compare(nfT, ofT, ApproxScalar(1e-15), false);

        compare(nfDivT, ofDivT, ApproxScalar(1e-15), false);
    }
}
