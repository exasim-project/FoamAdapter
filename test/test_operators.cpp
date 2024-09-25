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
#include "FoamAdapter/fvcc/mesh/fvccNeoMesh.hpp"

#define namespaceFoam // Suppress <using namespace Foam;>
#include "gaussConvectionScheme.H"

namespace fvcc = NeoFOAM::finiteVolume::cellCentred;

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

    std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec, runTime);
    Foam::fvccNeoMesh& mesh = *meshPtr;
    std::string execName = std::visit([](auto e) { return e.print(); }, exec);

    Foam::Info << "reading mesh with executor: " << execName << Foam::endl;
    NeoFOAM::UnstructuredMesh& uMesh = mesh.uMesh();

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);

    SECTION("scalar_" + execName)
    {
        Foam::IStringStream is("linear");
        Foam::tmp<Foam::surfaceInterpolationScheme<Foam::scalar>> foamInterPol =
            Foam::surfaceInterpolationScheme<Foam::scalar>::New(mesh, is);

        Foam::Info << "Reading field T\n" << Foam::endl;

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
        Foam::surfaceScalarField surfT(foamInterPol->interpolate(t));

        fvcc::VolumeField<NeoFOAM::scalar> neoT = constructFrom(exec, uMesh, t);
        neoT.correctBoundaryConditions();

        REQUIRE(neoT == t);

        fvcc::SurfaceField<NeoFOAM::scalar> neoSurfT = constructSurfaceField(exec, uMesh, surfT);

        SECTION("linear")
        {
            auto linearKernel = fvcc::SurfaceInterpolationFactory::create("linear", exec, uMesh);

            fvcc::SurfaceInterpolation interp(exec, uMesh, std::move(linearKernel));
            interp.interpolate(neoT, neoSurfT);
            auto sNeoSurfTHost = neoSurfT.internalField().copyToHost();
            std::span<Foam::scalar> surfTSpan(surfT.primitiveFieldRef().data(), surfT.size());
            REQUIRE_THAT(
                sNeoSurfTHost.span({0, surfT.size()}),
                Catch::Matchers::RangeEquals(surfTSpan, ApproxScalar(1e-15))
            );
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

    Foam::Info << "reading mesh with executor: " << execName << Foam::endl;
    std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec, runTime);
    Foam::fvccNeoMesh& mesh = *meshPtr;
    NeoFOAM::UnstructuredMesh& uMesh = mesh.uMesh();

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);

    SECTION("gauss_scalar_" + execName)
    {
        // linear interpolation hardcoded for now
        Foam::IStringStream is("linear");

        Foam::Info << "Reading field T\n" << Foam::endl;

        Foam::volScalarField t(
            Foam::IOobject(
                "T", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
            ),
            mesh
        );
        forAll(t, celli)
        {
            t[celli] = celli;
        }
        t.correctBoundaryConditions();
        t.write();
        std::span<Foam::scalar> sT(t.primitiveFieldRef().data(), t.size());

        Foam::fv::gaussGrad<Foam::scalar> foamGradScalar(mesh, is);
        Foam::volVectorField ofGradT("ofGradT", foamGradScalar.calcGrad(t, "test"));
        ofGradT.write();

        fvcc::VolumeField<NeoFOAM::scalar> neoT = constructFrom(exec, uMesh, t);
        neoT.correctBoundaryConditions();
        auto neoTHost = neoT.internalField().copyToHost();
        REQUIRE_THAT(neoTHost.span(), Catch::Matchers::RangeEquals(sT, ApproxScalar(1e-16)));

        fvcc::VolumeField<NeoFOAM::Vector> neoGradT = constructFrom(exec, uMesh, ofGradT);
        NeoFOAM::fill(neoGradT.internalField(), NeoFOAM::Vector(0.0, 0.0, 0.0));
        NeoFOAM::fill(neoGradT.boundaryField().value(), NeoFOAM::Vector(0.0, 0.0, 0.0));
        fvcc::GaussGreenGrad(exec, uMesh).grad(neoT, neoGradT);
        Foam::Info << "writing gradT field for exector: " << execName << Foam::endl;
        write(neoGradT.internalField(), mesh, "gradT_" + execName);

        std::span<Foam::vector> sOfGradT(ofGradT.primitiveFieldRef().data(), ofGradT.size());
        auto neoGradTHost = neoGradT.internalField().copyToHost();
        REQUIRE_THAT(
            neoGradTHost.span(), Catch::Matchers::RangeEquals(sOfGradT, ApproxVector(1e-12))
        );
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

    Foam::Info << "reading mesh with executor: " << execName << Foam::endl;
    std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec, runTime);
    Foam::fvccNeoMesh& mesh = *meshPtr;
    NeoFOAM::UnstructuredMesh& uMesh = mesh.uMesh();

    SECTION("gaussDiv_scalar_" + execName)
    {
        // linear interpolation hardcoded for now
        Foam::IStringStream is("linear");

        Foam::Info << "Reading field T\n" << Foam::endl;

        Foam::volScalarField t(
            Foam::IOobject(
                "T", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
            ),
            mesh
        );

        Foam::surfaceScalarField phi(
            Foam::IOobject(
                "phi", runTime.timeName(), mesh, Foam::IOobject::NO_READ, Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("phi", Foam::dimless, 0.0)
        );

        forAll(phi, facei)
        {
            phi[facei] = facei;
        }
        phi.write();

        forAll(t, celli)
        {
            t[celli] = celli;
        }
        t.correctBoundaryConditions();
        t.write();
        std::span<Foam::scalar> sT(t.primitiveFieldRef().data(), t.size());

        Foam::fv::gaussConvectionScheme<Foam::scalar> foamDivScalar(mesh, phi, is);
        Foam::volScalarField ofDivT("ofDivT", foamDivScalar.fvcDiv(phi, t));
        ofDivT.write();

        fvcc::VolumeField<NeoFOAM::scalar> neoT = constructFrom(exec, uMesh, t);


        fvcc::SurfaceField<NeoFOAM::scalar> neoPhi = constructSurfaceField(exec, uMesh, phi);
        std::span<Foam::scalar> sPhi(phi.primitiveFieldRef().data(), t.size());
        const auto sNeoPhiHost = neoPhi.internalField().copyToHost();

        REQUIRE_THAT(
            sNeoPhiHost.span({0, sPhi.size()}),
            Catch::Matchers::RangeEquals(sPhi, ApproxScalar(1e-15))
        );

        neoT.correctBoundaryConditions();
        auto neoTHost = neoT.internalField().copyToHost();
        REQUIRE_THAT(neoTHost.span(), Catch::Matchers::RangeEquals(sT, ApproxScalar(1e-15)));

        fvcc::VolumeField<NeoFOAM::scalar> neoDivT = constructFrom(exec, uMesh, ofDivT);
        NeoFOAM::fill(neoDivT.internalField(), 0.0);
        NeoFOAM::fill(neoDivT.boundaryField().value(), 0.0);
        fvcc::GaussGreenDiv(
            exec,
            uMesh,
            fvcc::SurfaceInterpolation(exec, uMesh, std::make_unique<fvcc::Linear>(exec, uMesh))
        )
            .div(neoDivT, neoPhi, neoT);
        Foam::Info << "writing divT field for exector: " << execName << Foam::endl;
        write(neoDivT.internalField(), mesh, "divT_" + execName);

        std::span<Foam::scalar> sOfDivT(ofDivT.primitiveFieldRef().data(), ofDivT.size());
        auto neoDivTHost = neoDivT.internalField().copyToHost();
        REQUIRE_THAT(
            neoDivTHost.span(), Catch::Matchers::RangeEquals(sOfDivT, ApproxScalar(1e-15))
        );
    }
}
