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
#include "NeoFOAM/finiteVolume/cellCentred/operators/fvccSparsityPattern.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/fvccLinearSystem.hpp"

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

    // SECTION("ddt_" + execName)
    // {

    //     auto ofT = randomScalarField(runTime, mesh);
    //     ofT.correctBoundaryConditions();


    //     fvcc::VolumeField<NeoFOAM::scalar>& nfT =
    //         fieldCol.registerField<fvcc::VolumeField<NeoFOAM::scalar>>(
    //             Foam::CreateFromFoamField<Foam::volScalarField> {
    //                 .exec = exec,
    //                 .nfMesh = nfMesh,
    //                 .foamField = ofT,
    //                 .name = "nfT"
    //             }
    //         );
    //     auto& nfTOld = fvcc::oldTime(nfT);
    //     const auto nfTOldSpan = nfTOld.internalField().span();
    //     NeoFOAM::map(
    //         nfTOld.internalField(),
    //         KOKKOS_LAMBDA(const std::size_t celli) { return nfTOldSpan[celli] - 1.0; }
    //     );

    //     ofT.oldTime() -= Foam::dimensionedScalar("value", Foam::dimTemperature, 1);
    //     ofT.oldTime().correctBoundaryConditions();

    //     Foam::fvScalarMatrix matrix(Foam::fvm::ddt(ofT));
    //     Foam::volScalarField ddt(
    //         "ddt",
    //         matrix & ofT
    //     ); // we should get a uniform field with a value of 1
    //     ddt.write();
    // }

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
        fvcc::LinearSystem<NeoFOAM::scalar> ls2(
            nfT,
            ls,
            fvcc::SparsityPattern::readOrCreate(nfMesh)
        );

        // check diag
        NeoFOAM::Field<NeoFOAM::scalar> diag(nfT.exec(), nfT.internalField().size(), 0.0);
        ls2.diag(diag);
        auto diagHost = diag.copyToHost();
        for (size_t i = 0; i < diagHost.size(); i++)
        {
            REQUIRE(diagHost[i] == coeff);
        }
        auto result = ls2 & nfT;
        auto implicitHost = result.internalField().copyToHost();
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
    //     Foam::Info << "matrix.diag() " << matrix.diag() << Foam::endl;
    //     Foam::Info << "matrix.upper() " << matrix.upper() << Foam::endl;
    //     Foam::Info << "matrix.lower() " << matrix.lower() << Foam::endl;

    //     Foam::volScalarField imp_laplacian("imp_laplacian", matrix & ofT);
    //     Foam::Info << "imp_laplacian: \n" << imp_laplacian.primitiveField() << Foam::endl;

    //     Foam::volScalarField exp_laplacian("exp_laplacian", Foam::fvc::laplacian(ofT));
    //     Foam::Info << "exp_laplacian: \n" << exp_laplacian.primitiveField() << Foam::endl;

    //     auto diff = imp_laplacian.primitiveField() - exp_laplacian.primitiveField();
    //     Foam::Info << "diff: \n" << diff << Foam::endl;

    // }

    SECTION("div_" + execName)
    {
        auto ofT = randomScalarField(runTime, mesh);
        forAll(ofT, celli)
        {
            ofT[celli] = celli;
        }
        ofT.correctBoundaryConditions();
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
            // ofPhi[facei] = facei;
            ofPhi[facei] = 1;
        }
        // FIXME
        // the boundary conditions ofPhi are zero

        auto nfPhi = constructSurfaceField(exec, nfMesh, ofPhi);
        // the boundary conditions nfPhi are not zero


        Foam::fvScalarMatrix matrix(Foam::fvm::div(ofPhi, ofT));
        Foam::volScalarField imp_div("imp_div", matrix & ofT);

        NeoFOAM::TokenList scheme({std::string("linear")});
        fvcc::SparsityPattern pattern(nfMesh);
        la::LinearSystem<NeoFOAM::scalar, NeoFOAM::localIdx> ls = pattern.linearSystem();
        fvcc::GaussGreenDiv(exec, nfMesh, scheme).div(ls, nfPhi, nfT);

        fvcc::LinearSystem<NeoFOAM::scalar> fvLs(
            nfT,
            ls,
            fvcc::SparsityPattern::readOrCreate(nfMesh)
        );

        auto values = ls.matrix().values();
        auto colIdx = ls.matrix().colIdxs();
        auto rowPtrs = ls.matrix().rowPtrs();

        auto result = fvLs & nfT;
        auto implicitHost = result.internalField().copyToHost();
        const auto vol = nfMesh.cellVolumes().span();

        Foam::Info << "matrix.lower() " << matrix.lower() << Foam::endl;

        Foam::scalarField Diag = matrix.diag();
        Diag *= 0;
        std::span<Foam::scalar> ds(Diag.data(), Diag.size());

        Foam::labelUList& l = const_cast<Foam::labelUList&>(matrix.lduAddr().lowerAddr());
        Foam::labelUList& u = const_cast<Foam::labelUList&>(matrix.lduAddr().upperAddr());
        std::span<Foam::label> ls2(l.data(), l.size());
        std::span<Foam::label> us(u.data(), u.size());
        auto ownerSpan = nfMesh.faceOwner().span();
        auto neighbourSpan = nfMesh.faceNeighbour().span();
        std::span<Foam::scalar> lowerSpan(matrix.lower().data(), matrix.lower().size());
        std::span<Foam::scalar> upperSpan(matrix.upper().data(), matrix.upper().size());
        std::span<Foam::scalar> diagSpan(matrix.diag().data(), matrix.diag().size());

        for (Foam::label face = 0; face < l.size(); face++)
        {
            auto lower = matrix.lower()[face];
            auto upper = matrix.upper()[face];
            auto lf = l[face];
            auto uf = u[face];
            Diag[l[face]] -= matrix.lower()[face];
            Diag[u[face]] -= matrix.upper()[face];
            lf = 0;
        }

        for (size_t celli = 0; celli < implicitHost.size(); celli++)
        {
            // will fail no boundary conditions
            std::cout << "celli: " << celli << " " << implicitHost[celli] << " "
                      << imp_div[celli] * vol[celli] << std::endl;
            Foam::Info << "   cellCells[celli]: " << mesh.cellCells()[celli] << Foam::endl;
            Foam::Info << "   matrix.diag celli: " << matrix.diag()[celli] << Foam::endl;

            for (size_t coli = rowPtrs[celli]; coli < rowPtrs[celli + 1]; coli++)
            {
                Foam::Info << "    colIdx: " << colIdx[coli] << " value: " << values[coli]
                           << " x: " << Foam::endl;
            }
        }
        // REQUIRE(implicitHost[13] / vol[13] == imp_div[13]); // <-- only with no boundary
        // conditions
    }
    REQUIRE(false);
}
