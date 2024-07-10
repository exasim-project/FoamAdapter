// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/catch_approx.hpp>

#include "NeoFOAM/fields/Field.hpp"
#include "NeoFOAM/fields/FieldOperations.hpp"
#include "NeoFOAM/fields/FieldTypeDefs.hpp"
#include "NeoFOAM/fields/comparisions/fieldComparision.hpp"

#include "NeoFOAM/fields/boundaryFields.hpp"
#include "NeoFOAM/fields/domainField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/fields/fvccVolField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/fields/fvccSurfaceField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/fvccBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/vol/scalar/fvccScalarFixedValueBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/vol/scalar/fvccScalarZeroGradientBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/surface/scalar/fvccSurfaceScalarCalculatedBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/surface/scalar/fvccSurfaceScalarEmptyBoundaryField.hpp"

#include "NeoFOAM/cellCentredFiniteVolume/grad/gaussGreenGrad.hpp"

#include "NeoFOAM/cellCentredFiniteVolume/surfaceInterpolation/linear.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/surfaceInterpolation/upwind.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/surfaceInterpolation/surfaceInterpolation.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/surfaceInterpolation/surfaceInterpolationSelector.hpp"

#include "NeoFOAM/cellCentredFiniteVolume/div/gaussGreenDiv.hpp"

#include "FoamAdapter/fvcc/surfaceInterpolation/surfaceInterpolationFactory.hpp"
#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/writers/writers.hpp"
#include "FoamAdapter/comparison/fieldComparison.hpp"
#include "FoamAdapter/setup/setup.hpp"
#include "FoamAdapter/fvcc/mesh/fvccNeoMesh.hpp"

#include <random>
#include <span>

#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"
#include "gaussConvectionScheme.H"

Foam::Time* timePtr;    // A single time object
Foam::argList* argsPtr; // Some forks want argList access at createMesh.H

int main(int argc, char* argv[])
{
    // Initialize Catch2
    Kokkos::initialize(argc, argv);
    Catch::Session session;

    // Specify command line options
    int returnCode = session.applyCommandLine(argc, argv);
    if (returnCode != 0) // Indicates a command line error
        return returnCode;

#include "setRootCase.H"
#include "createTime.H"
    argsPtr = &args;
    timePtr = &runTime;

    int result = session.run();

    // Run benchmarks if there are any
    Kokkos::finalize();

    return result;
}

struct ApproxScalar
{
    Foam::scalar margin;
    bool operator()(double rhs, double lhs) const
    {
        return Catch::Approx(rhs).margin(margin) == lhs;
    }
};

struct ApproxVector
{
    Foam::scalar margin;
    bool operator()(NeoFOAM::Vector rhs, Foam::vector lhs) const
    {
        NeoFOAM::Vector diff(rhs[0] - lhs[0], rhs[1] - lhs[1], rhs[2] - lhs[2]);

        return Catch::Approx(0).margin(margin) == mag(diff);
    }
};

TEST_CASE("Interpolation")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    NeoFOAM::executor exec = GENERATE(
        NeoFOAM::executor(NeoFOAM::CPUExecutor {}),
        NeoFOAM::executor(NeoFOAM::OMPExecutor {}),
        NeoFOAM::executor(NeoFOAM::GPUExecutor {})
    );

    std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec, runTime);
    Foam::fvccNeoMesh& mesh = *meshPtr;
    std::string exec_name = std::visit([](auto e) { return e.print(); }, exec);

    Foam::Info << "reading mesh with executor: " << exec_name << Foam::endl;
    NeoFOAM::unstructuredMesh& uMesh = mesh.uMesh();

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);

    SECTION("scalar_" + exec_name)
    {
        Foam::IStringStream is("linear");
        Foam::tmp<Foam::surfaceInterpolationScheme<Foam::scalar>> foamInterPol =
            Foam::surfaceInterpolationScheme<Foam::scalar>::New(mesh, is);

        Foam::Info << "Reading field T\n" << Foam::endl;

        Foam::volScalarField T(
            Foam::IOobject(
                "T", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
            ),
            mesh
        );

        forAll(T, celli)
        {
            T[celli] = dis(gen);
        }

        T.correctBoundaryConditions();
        Foam::surfaceScalarField surfT(foamInterPol->interpolate(T));

        NeoFOAM::fvccVolField<NeoFOAM::scalar> neoT = constructFrom(exec, uMesh, T);
        neoT.correctBoundaryConditions();

        REQUIRE(neoT == T);
        std::vector<std::unique_ptr<NeoFOAM::fvccSurfaceBoundaryField<NeoFOAM::scalar>>> bcs;
        bcs.push_back(std::make_unique<NeoFOAM::fvccSurfaceScalarCalculatedBoundaryField>(uMesh, 0)
        );
        bcs.push_back(std::make_unique<NeoFOAM::fvccSurfaceScalarEmptyBoundaryField>(uMesh, 1));
        NeoFOAM::fvccSurfaceField<NeoFOAM::scalar> neoSurfT(exec, uMesh, std::move(bcs));

        SECTION("linear")
        {
            // std::unique_ptr<NeoFOAM::surfaceInterpolationKernel> linearKernel(new
            // NeoFOAM::linear(exec, uMesh));

            NeoFOAM::surfaceInterpolation interp(
                NeoFOAM::surfaceInterpolationSelector(std::string("linear"), exec, mesh.uMesh())
            );
            interp.interpolate(neoSurfT, neoT);
            auto s_neoSurfT = neoSurfT.internalField().copyToHost().field();
            std::span<Foam::scalar> surfT_span(surfT.primitiveFieldRef().data(), surfT.size());
            REQUIRE_THAT(
                s_neoSurfT.subspan(0, surfT.size()),
                Catch::Matchers::RangeEquals(surfT_span, ApproxScalar(1e-15))
            );
        }
    }
}

TEST_CASE("GradOperator")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    NeoFOAM::executor exec = GENERATE(
        NeoFOAM::executor(NeoFOAM::CPUExecutor {}),
        NeoFOAM::executor(NeoFOAM::OMPExecutor {}),
        NeoFOAM::executor(NeoFOAM::GPUExecutor {})
    );

    // NeoFOAM::executor exec = NeoFOAM::CPUExecutor{};

    std::string exec_name = std::visit([](auto e) { return e.print(); }, exec);

    Foam::Info << "reading mesh with executor: " << exec_name << Foam::endl;
    std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec, runTime);
    Foam::fvccNeoMesh& mesh = *meshPtr;
    NeoFOAM::unstructuredMesh& uMesh = mesh.uMesh();

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);

    SECTION("gauss_scalar_" + exec_name)
    {
        // linear interpolation hardcoded for now
        Foam::IStringStream is("linear");

        Foam::Info << "Reading field T\n" << Foam::endl;

        Foam::volScalarField T(
            Foam::IOobject(
                "T", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
            ),
            mesh
        );
        forAll(T, celli)
        {
            // T[celli] = dis(gen);
            T[celli] = celli;
        }
        T.correctBoundaryConditions();
        T.write();
        std::span<Foam::scalar> s_T(T.primitiveFieldRef().data(), T.size());

        Foam::fv::gaussGrad<Foam::scalar> foamGradScalar(mesh, is);
        Foam::volVectorField ofGradT("ofGradT", foamGradScalar.calcGrad(T, "test"));
        ofGradT.write();

        NeoFOAM::fvccVolField<NeoFOAM::scalar> neoT = constructFrom(exec, uMesh, T);
        neoT.correctBoundaryConditions();
        REQUIRE_THAT(
            neoT.internalField().copyToHost().field(),
            Catch::Matchers::RangeEquals(s_T, ApproxScalar(1e-16))
        );

        NeoFOAM::fvccVolField<NeoFOAM::Vector> neoGradT = constructFrom(exec, uMesh, ofGradT);
        NeoFOAM::fill(neoGradT.internalField(), NeoFOAM::Vector(0.0, 0.0, 0.0));
        NeoFOAM::fill(neoGradT.boundaryField().value(), NeoFOAM::Vector(0.0, 0.0, 0.0));
        NeoFOAM::gaussGreenGrad(exec, uMesh).grad(neoGradT, neoT);
        Foam::Info << "writing gradT field for exector: " << exec_name << Foam::endl;
        write(neoGradT.internalField(), mesh, "gradT_" + exec_name);

        std::span<Foam::vector> s_ofGradT(ofGradT.primitiveFieldRef().data(), ofGradT.size());
        REQUIRE_THAT(
            neoGradT.internalField().copyToHost().field(),
            Catch::Matchers::RangeEquals(s_ofGradT, ApproxVector(1e-12))
        );
    }
}


TEST_CASE("DivOperator")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    NeoFOAM::executor exec = GENERATE(
        NeoFOAM::executor(NeoFOAM::CPUExecutor {}),
        NeoFOAM::executor(NeoFOAM::OMPExecutor {}),
        NeoFOAM::executor(NeoFOAM::GPUExecutor {})
    );

    std::string exec_name = std::visit([](auto e) { return e.print(); }, exec);

    Foam::Info << "reading mesh with executor: " << exec_name << Foam::endl;
    std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec, runTime);
    Foam::fvccNeoMesh& mesh = *meshPtr;
    NeoFOAM::unstructuredMesh& uMesh = mesh.uMesh();

    SECTION("gaussDiv_scalar_" + exec_name)
    {
        // linear interpolation hardcoded for now
        Foam::IStringStream is("linear");

        Foam::Info << "Reading field T\n" << Foam::endl;

        Foam::volScalarField T(
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
            // T[celli] = dis(gen);
            phi[facei] = facei;
        }
        phi.write();

        forAll(T, celli)
        {
            // T[celli] = dis(gen);
            T[celli] = celli;
        }
        T.correctBoundaryConditions();
        T.write();
        std::span<Foam::scalar> s_T(T.primitiveFieldRef().data(), T.size());

        Foam::fv::gaussConvectionScheme<Foam::scalar> foamDivScalar(mesh, phi, is);
        Foam::volScalarField ofDivT("ofDivT", foamDivScalar.fvcDiv(phi, T));
        ofDivT.write();

        NeoFOAM::fvccVolField<NeoFOAM::scalar> neoT = constructFrom(exec, uMesh, T);


        NeoFOAM::fvccSurfaceField<NeoFOAM::scalar> neoPhi = constructSurfaceField(exec, uMesh, phi);
        std::span<Foam::scalar> s_phi(phi.primitiveFieldRef().data(), T.size());
        const auto s_neoPhi_host = neoPhi.internalField().copyToHost().field();
        REQUIRE_THAT(
            s_neoPhi_host.subspan(0, s_phi.size()),
            Catch::Matchers::RangeEquals(s_phi, ApproxScalar(1e-15))
        );

        neoT.correctBoundaryConditions();
        REQUIRE_THAT(
            neoT.internalField().copyToHost().field(),
            Catch::Matchers::RangeEquals(s_T, ApproxScalar(1e-15))
        );

        NeoFOAM::fvccVolField<NeoFOAM::scalar> neoDivT = constructFrom(exec, uMesh, ofDivT);
        NeoFOAM::fill(neoDivT.internalField(), 0.0);
        NeoFOAM::fill(neoDivT.boundaryField().value(), 0.0);
        NeoFOAM::gaussGreenDiv(
            exec, uMesh, NeoFOAM::surfaceInterpolationSelector(std::string("linear"), exec, uMesh)
        )
            .div(neoDivT, neoPhi, neoT);
        Foam::Info << "writing divT field for exector: " << exec_name << Foam::endl;
        write(neoDivT.internalField(), mesh, "divT_" + exec_name);

        std::span<Foam::scalar> s_ofDivT(ofDivT.primitiveFieldRef().data(), ofDivT.size());
        REQUIRE_THAT(
            neoDivT.internalField().copyToHost().field(),
            Catch::Matchers::RangeEquals(s_ofDivT, ApproxScalar(1e-15))
        );
    }
}
