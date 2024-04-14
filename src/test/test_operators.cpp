// SPDX-License-Identifier: GPL-3.0
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
#include "NeoFOAM/cellCentredFiniteVolume/surfaceInterpolation/surfaceInterpolationFactory.hpp"

#include "FoamAdapter/fvcc/surfaceInterpolation/surfaceInterpolationFactory.hpp"
#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/writers/writers.hpp"
#include "FoamAdapter/comparision/fieldComparision.hpp"

#include <random>
#include <span>

#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"

Foam::Time *timePtr;    // A single time object
Foam::argList *argsPtr; // Some forks want argList access at createMesh.H
Foam::fvMesh *meshPtr;  // A single mesh object




int main(int argc, char *argv[])
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
#include "createMesh.H"
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
    bool operator()(double rhs,double lhs) const
    {
        return Catch::Approx(rhs).margin(margin) == lhs;
    }
};

TEST_CASE("Interpolation")
{
    Foam::Time &runTime = *timePtr;
    Foam::argList &args = *argsPtr;
#include "createMesh.H"

    NeoFOAM::executor exec = GENERATE(
        NeoFOAM::executor(NeoFOAM::CPUExecutor{}),
        NeoFOAM::executor(NeoFOAM::OMPExecutor{}),
        NeoFOAM::executor(NeoFOAM::GPUExecutor{})
        );
    std::string exec_name = std::visit([](auto e)
                                       { return e.print(); },
                                       exec);

    Foam::Info << "reading mesh" << Foam::endl;
    NeoFOAM::unstructuredMesh uMesh = readOpenFOAMMesh(exec, mesh);

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);

    SECTION("scalar_" + exec_name)
    {
        Foam::IStringStream is("linear");
        Foam::tmp<Foam::surfaceInterpolationScheme<Foam::scalar>> foamInterPol = Foam::surfaceInterpolationScheme<Foam::scalar>::New(mesh, is);

        Foam::Info
            << "Reading field T\n"
            << Foam::endl;

        Foam::volScalarField T(
            Foam::IOobject(
                "T",
                runTime.timeName(),
                mesh,
                Foam::IOobject::MUST_READ,
                Foam::IOobject::AUTO_WRITE),
            mesh);

        forAll(T, celli)
        {
            T[celli] = dis(gen);
        }
        T.correctBoundaryConditions();
        Foam::surfaceScalarField surfT = foamInterPol->interpolate(T);

        NeoFOAM::fvccVolField<NeoFOAM::scalar> neoT = constructFrom(exec, uMesh, T);
        neoT.correctBoundaryConditions();

        REQUIRE(neoT == T);
        std::vector<std::unique_ptr<NeoFOAM::fvccSurfaceBoundaryField<NeoFOAM::scalar>>> bcs;
        bcs.push_back(std::make_unique<NeoFOAM::fvccSurfaceScalarCalculatedBoundaryField>(uMesh, 0));
        bcs.push_back(std::make_unique<NeoFOAM::fvccSurfaceScalarEmptyBoundaryField>(uMesh, 1));
        NeoFOAM::fvccSurfaceField<NeoFOAM::scalar> neoSurfT(
            exec,
            uMesh,
            std::move(bcs));

        Foam::Info << "Creating linear kernel" << Foam::endl;
        if (Foam::surfaceInterpolationFactory::dictionaryConstructorTablePtr_ != nullptr) {
            Foam::Info << Foam::surfaceInterpolationFactory::dictionaryConstructorTablePtr_->sortedToc() << Foam::endl;
        }
        else {
            Foam::Info << "dictionaryConstructorTablePtr_ is nullptr" << Foam::endl;
        }
        // Foam::surfaceInterpolationFactory::New(exec, uMesh);
        // NeoFOAM::surfaceInterpolationFactory surfInterFactory{};
        // surfInterFactory.registerClass("linear", [](const NeoFOAM::executor& exec, const NeoFOAM::unstructuredMesh& mesh) {
        //     return std::make_unique<NeoFOAM::linear>(exec, mesh);
        // });
        // std::cout << "Creating linear kernel" << std::endl;
        
        // std::cout << "linear::registered: " << NeoFOAM::linear::registered << std::endl;
        std::cout << "map size: " << NeoFOAM::surfaceInterpolationFactory::number_of_instances() << std::endl;
        // surfInterFactory.registerClass("upwind", [](const NeoFOAM::executor& exec, const NeoFOAM::unstructuredMesh& mesh) {
        //     return std::make_unique<NeoFOAM::upwind>(exec, mesh);
        // });

        // std::unique_ptr<NeoFOAM::surfaceInterpolationKernel> linearKernel = NeoFOAM::surfaceInterpolationFactory::New("linear", exec, uMesh);

        std::unique_ptr<NeoFOAM::surfaceInterpolationKernel> linearKernel(new NeoFOAM::linear(exec, uMesh));

        NeoFOAM::surfaceInterpolation interp(exec, uMesh, std::move(linearKernel));
        interp.interpolate(neoSurfT, neoT);
        auto s_neoSurfT = neoSurfT.internalField().copyToHost().field();
        std::span<Foam::scalar> surfT_span(surfT.primitiveFieldRef().data(), surfT.size());
        REQUIRE_THAT(s_neoSurfT.subspan(0,surfT.size()), Catch::Matchers::RangeEquals(surfT_span, ApproxScalar(1e-12)));
    }
}

TEST_CASE("GradOperator")
{
    Foam::Time &runTime = *timePtr;
    Foam::argList &args = *argsPtr;
#include "createMesh.H"
    NeoFOAM::executor exec = NeoFOAM::CPUExecutor();
    Foam::Info << "reading mesh" << Foam::endl;
    NeoFOAM::unstructuredMesh uMesh = readOpenFOAMMesh(exec, mesh);

    Foam::Info << "Reading field T\n"
               << Foam::endl;

    Foam::volScalarField T(
        Foam::IOobject(
            "T",
            runTime.timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::AUTO_WRITE),
        mesh);

    Foam::IStringStream is("linear");
    Foam::fv::gaussGrad<Foam::scalar> foamGradScalar(mesh, is);
    Foam::volVectorField ofGradT = foamGradScalar.calcGrad(T, "test");

    NeoFOAM::fvccVolField<NeoFOAM::scalar> neoT = constructFrom(exec, uMesh, T);
    NeoFOAM::fill(neoT.internalField(), 1.0);
    Foam::scalar pi = Foam::constant::mathematical::pi;
    const NeoFOAM::Field<NeoFOAM::Vector> &cc = uMesh.cellCentres();
    Foam::scalar spread = 0.05;
    auto s_cc = cc.field();

    neoT.internalField().apply(KOKKOS_LAMBDA(int i) {
        return std::exp(-0.5 * (std::pow((s_cc[i][0] - 0.05) / spread, 2.0) + std::pow((s_cc[i][1] - 0.075) / spread, 2.0)));
    });
    neoT.correctBoundaryConditions();

    NeoFOAM::vectorField nofGradT = NeoFOAM::gaussGreenGrad(exec, uMesh).grad(neoT.internalField());
}