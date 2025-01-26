// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "NeoFOAM/core/primitives/scalar.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/fields/volumeField.hpp"
#include <cstddef>
#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenGrad.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenDiv.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/sourceTerm.hpp"

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


TEST_CASE("matrix multiplication")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    NeoFOAM::Database db;
    fvcc::FieldCollection& fieldCol = fvcc::FieldCollection::instance(db, "fieldCollection");

    // NeoFOAM::Executor exec = GENERATE(
    //     NeoFOAM::Executor(NeoFOAM::SerialExecutor {}),
    //     NeoFOAM::Executor(NeoFOAM::CPUExecutor {}),
    //     NeoFOAM::Executor(NeoFOAM::GPUExecutor {})
    // );
    NeoFOAM::Executor exec = NeoFOAM::SerialExecutor {};

    std::string execName = std::visit([](auto e) { return e.name(); }, exec);

    auto meshPtr = Foam::createMesh(exec, runTime);
    Foam::MeshAdapter& mesh = *meshPtr;
    auto nfMesh = mesh.nfMesh();
    runTime.setDeltaT(1);

    SECTION("ddt_" + execName)
    {

        auto ofT = randomScalarField(runTime, mesh);
        ofT.correctBoundaryConditions();


        fvcc::VolumeField<NeoFOAM::scalar>& nfT =
            fieldCol.registerField<fvcc::VolumeField<NeoFOAM::scalar>>(
                Foam::CreateFromFoamField<Foam::volScalarField> {
                    .exec = exec,
                    .nfMesh = nfMesh,
                    .foamField = ofT,
                    .name = "nfT"
                }
            );
        auto& nfTOld = fvcc::oldTime(nfT);
        const auto nfTOldSpan = nfTOld.internalField().span();
        NeoFOAM::map(
            nfTOld.internalField(),
            KOKKOS_LAMBDA(const std::size_t celli) { return nfTOldSpan[celli] - 1.0; }
        );

        ofT.oldTime() -= Foam::dimensionedScalar("value", Foam::dimTemperature, 1);
        ofT.oldTime().correctBoundaryConditions();

        Foam::fvScalarMatrix matrix(Foam::fvm::ddt(ofT));
        Foam::volScalarField ddt(
            "ddt",
            matrix & ofT
        ); // we should get a uniform field with a value of 1
        ddt.write();
    }

    SECTION("sourceterm_" + execName)
    {
        NeoFOAM::scalar coeff = 2.0;
        auto ofT = randomScalarField(runTime, mesh);
        fvcc::VolumeField<NeoFOAM::scalar> nfT = constructFrom(exec, nfMesh, ofT);

        NeoFOAM::map(
            nfT.internalField(),
            KOKKOS_LAMBDA(const std::size_t celli) { return celli; }
        );
        auto coefficients = nfT;
        NeoFOAM::fill(coefficients.internalField(), coeff);
        fvcc::SourceTerm sourceTerm(dsl::Operator::Type::Implicit, coefficients, nfT);
        NeoFOAM::Field<NeoFOAM::scalar> source(nfT.exec(), nfT.internalField().size(), 0.0);
        sourceTerm.explicitOperation(source);

        auto sourceHost = source.copyToHost();
        auto nftHost = nfT.internalField().copyToHost();
        for (size_t i = 0; i < sourceHost.size(); i++)
        {
            REQUIRE(sourceHost[i] == coeff * nftHost[i]);
        }

        auto ls = sourceTerm.createEmptyLinearSystem();
        sourceTerm.implicitOperation(ls);
        auto implicitHost = NeoFOAM::la::SpMV(ls, nfT.internalField()).copyToHost();
        for (size_t i = 0; i < implicitHost.size(); i++)
        {
            REQUIRE(implicitHost[i] == coeff * nftHost[i]);
        }
    }

    // SECTION("laplacian_" + execName)
    // {
    //     auto ofT = randomScalarField(runTime, mesh);
    //     ofT.correctBoundaryConditions();

    //     Foam::fvScalarMatrix matrix(Foam::fvm::laplacian(ofT));

    //     Foam::volScalarField imp_laplacian("imp_laplacian", matrix & ofT);
    //     Foam::Info << "imp_laplacian: \n" << imp_laplacian.primitiveField() << Foam::endl;

    //     Foam::volScalarField exp_laplacian("exp_laplacian", Foam::fvc::laplacian(ofT));
    //     Foam::Info << "exp_laplacian: \n" << exp_laplacian.primitiveField() << Foam::endl;

    //     auto diff = imp_laplacian.primitiveField() - exp_laplacian.primitiveField();
    //     Foam::Info << "diff: \n" << diff << Foam::endl;

    // }
    // REQUIRE(false);
}
