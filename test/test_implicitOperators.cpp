// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include <cstddef>
#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenGrad.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenDiv.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/ddtOperator.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/sourceTerm.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/sparsityPattern.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/linearSystem.hpp"

#include "gaussConvectionScheme.H"
#include "NeoFOAM/dsl.hpp"

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
        // ddt.write();
        fvcc::DdtOperator ddtOp(dsl::Operator::Type::Implicit, nfT);

        auto ls = ddtOp.createEmptyLinearSystem();
        ddtOp.implicitOperation(ls, runTime.value(), runTime.deltaTValue());
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
            REQUIRE(diagHost[i] == 1.0);
        }
        auto result = ls2 & nfT;
        auto implicitHost = result.internalField().copyToHost();
        for (size_t i = 0; i < implicitHost.size(); i++)
        {
            REQUIRE(implicitHost[i] == Catch::Approx(ddt[i]).margin(1e-16));
        }
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
        fvcc::LinearSystem<NeoFOAM::scalar> ls2(
            nfT,
            ls,
            fvcc::SparsityPattern::readOrCreate(nfMesh)
        );

        // check diag
        NeoFOAM::Field<NeoFOAM::scalar> diag(nfT.exec(), nfT.internalField().size(), 0.0);
        ls2.diag(diag);
        auto diagHost = diag.copyToHost();
        auto cellVolumes = nfMesh.cellVolumes().copyToHost();
        for (size_t i = 0; i < diagHost.size(); i++)
        {
            REQUIRE(diagHost[i] == coeff * cellVolumes[i]);
        }
        auto result = ls2 & nfT;
        auto implicitHost = result.internalField().copyToHost();
        for (size_t i = 0; i < implicitHost.size(); i++)
        {
            REQUIRE(implicitHost[i] == coeff * nftHost[i] * cellVolumes[i]);
        }
    }

    SECTION("solve sourceterm_" + execName)
    {
        auto ofT = randomScalarField(runTime, mesh);
        ofT.primitiveFieldRef() = 1.0;
        ofT.correctBoundaryConditions();

        fvcc::VolumeField<NeoFOAM::scalar> nfT = constructFrom(exec, nfMesh, ofT);
        NeoFOAM::fill(nfT.internalField(), 1.0);

        auto nfCoeff1 = nfT;
        NeoFOAM::fill(nfCoeff1.internalField(), 1.0);

        auto nfCoeff2 = nfT;
        NeoFOAM::fill(nfCoeff2.internalField(), 2.0);

        Foam::dimensionedScalar coeff1("coeff", Foam::dimless, 1.0);
        Foam::dimensionedScalar coeff2("coeff", Foam::dimless, 2.0);
        Foam::fvScalarMatrix matrix(Foam::fvm::Sp(coeff1, ofT) - coeff2 * ofT);


        dsl::Expression eqnSys(dsl::imp::Source(nfCoeff1, nfT) - dsl::exp::Source(nfCoeff2, nfT));
        NeoFOAM::scalar t = 0;
        NeoFOAM::scalar dt = 1;
        NeoFOAM::Dictionary fvSchemesDict {};
        NeoFOAM::Dictionary fvSolutionDict {};
        fvSolutionDict.insert("maxIters", 100);
        fvSolutionDict.insert("relTol", float(1e-7));


        dsl::solve(eqnSys, nfT, t, dt, fvSchemesDict, fvSolutionDict);
        matrix.solve();

        auto nfTHost = nfT.internalField().copyToHost();
        for (size_t celli = 0; celli < nfTHost.size(); celli++)
        {
            REQUIRE(nfTHost[celli] == Catch::Approx(ofT[celli]).margin(1e-16));
        }
    }


    SECTION("solve div" + execName)
    {
        auto ofT = randomScalarField(runTime, mesh);
        forAll(ofT, celli)
        {
            ofT[celli] = celli;
        }
        ofT.correctBoundaryConditions();

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
            ofPhi[facei] = 1;
        }

        auto nfPhi = constructSurfaceField(exec, nfMesh, ofPhi);
        auto nfCoeff1 = nfT;
        NeoFOAM::scalar coeff = 1000;
        Foam::dimensionedScalar coeff1("coeff", Foam::dimensionSet(0, -3, 0, 0, 0), coeff);
        Foam::dimensionedScalar coeff2("coeff", Foam::dimensionSet(0, -3, 0, 0, 0), coeff);
        NeoFOAM::fill(nfCoeff1.internalField(), coeff);

        auto nfCoeff2 = nfT;
        NeoFOAM::fill(nfCoeff2.internalField(), coeff);

        Foam::volScalarField testfvcDiv(Foam::fvc::div(ofPhi, ofT));
        std::span<Foam::scalar> testfvcDivSpan(
            testfvcDiv.primitiveFieldRef().data(),
            testfvcDiv.primitiveFieldRef().size()
        );

        Foam::fvScalarMatrix matrix(
            Foam::fvm::Sp(coeff1, ofT) + Foam::fvm::div(ofPhi, ofT) - coeff2 * ofT
        );

        Foam::fvScalarMatrix lap(Foam::fvm::laplacian(ofT));

        dsl::Expression eqnSys(
            dsl::imp::Source(nfCoeff1, nfT) + dsl::imp::div(nfPhi, nfT)
            - dsl::exp::Source(nfCoeff2, nfT)
        );

        NeoFOAM::scalar t = 0;
        NeoFOAM::scalar dt = 1;
        NeoFOAM::Dictionary fvSchemesDict {};
        NeoFOAM::Dictionary divSchemes {};
        divSchemes.insert(
            "div(phi,T)",
            NeoFOAM::TokenList {std::string("Gauss"), std::string("linear")}
        );
        fvSchemesDict.insert("divSchemes", divSchemes);

        NeoFOAM::Dictionary fvSolutionDict {};
        fvSolutionDict.insert("maxIters", 100);
        fvSolutionDict.insert("relTol", float(1e-8));

        matrix.solve();
        dsl::solve(eqnSys, nfT, t, dt, fvSchemesDict, fvSolutionDict);

        std::span<Foam::scalar> ofTSpan(ofT.data(), ofT.size());
        auto nfTHost = nfT.internalField().copyToHost();
        std::span<NeoFOAM::scalar> nfTHostSpan(nfTHost.data(), nfTHost.size());
        for (size_t celli = 0; celli < nfTHost.size(); celli++)
        {
            REQUIRE(nfTHost[celli] == Catch::Approx(ofT[celli]).margin(1e-16));
        }
    }
}
