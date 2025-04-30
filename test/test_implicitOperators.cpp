// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023-2025 NeoFOAM authors

#include <cstddef>
#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "common.hpp"

#include "gaussConvectionScheme.H"


namespace fvcc = NeoN::finiteVolume::cellCentred;
namespace dsl = NeoN::dsl;

extern Foam::Time* timePtr;    // A single time object
extern Foam::argList* argsPtr; // Some forks want argList access at createMesh.H
extern Foam::fvMesh* meshPtr;  // A single mesh object


TEST_CASE("Implicit operators")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    NeoN::Database db;
    fvcc::VectorCollection& fieldCol = fvcc::VectorCollection::instance(db, "VectorCollection");

    auto [execName, exec] = GENERATE(allAvailableExecutor());

    auto meshPtr = Foam::createMesh(exec, runTime);
    Foam::MeshAdapter& mesh = *meshPtr;
    auto nfMesh = mesh.nfMesh();

    runTime.setDeltaT(1);

    auto ofT = randomScalarField(runTime, mesh, "T");
    ofT.correctBoundaryConditions();

    fvcc::VolumeField<NeoN::scalar>& nfT = fieldCol.registerVector<fvcc::VolumeField<NeoN::scalar>>(
        Foam::CreateFromFoamField<Foam::volScalarField> {
            .exec = exec,
            .nfMesh = nfMesh,
            .foamField = ofT,
            .name = "nfT"

        }
    );

    NeoN::Dictionary fvSchemesDict {};
    NeoN::Dictionary fvSolutionDict {};

    /* this tests computes the residual ie Ax-b for OF and NF and compares the results
     *
     * Since b is not initialised it is zero. A is computed from u^1 and u^1 with a  difference
     * of 1, thus ddt(U) = 1
     *
     * res Ax -b
     */
    SECTION("ddt operator on " + execName)
    {
        // setup OpenFOAM
        ofT.oldTime() -= Foam::dimensionedScalar("value", Foam::dimTemperature, 1);
        ofT.oldTime().correctBoundaryConditions();
        Foam::fvScalarMatrix matrix(Foam::fvm::ddt(ofT));
        Foam::volScalarField ofRes("ofRes", matrix & ofT);

        // setup NeoFOAM
        auto& nfTOld = fvcc::oldTime(nfT);
        nfTOld -= 1.0;
        nfTOld.correctBoundaryConditions();

        auto expr = dsl::Expression<NeoN::scalar>(exec);
        expr.addOperator(NeoN::dsl::imp::ddt(nfT));

        auto rhs = fvcc::VolumeField<NeoN::scalar>(nfT);
        NeoN::fill(rhs.internalVector(), 0.0);

        auto ls =
            NeoN::finiteVolume::cellCentred::assembleLinearSystem(expr, rhs, fvSchemesDict, 0, 1.0);

        auto res = ls & nfT;
        auto resHost = res.internalVector().copyToHost();
        for (size_t i = 0; i < resHost.size(); i++)
        {
            REQUIRE(resHost.view()[i] == Catch::Approx(ofRes[i]).margin(1e-16));
        }
    }

    SECTION("sourceterm_" + execName)
    {
        auto coefficients = nfT;

        NeoN::scalar coeff = 2.0;
        NeoN::fill(coefficients.internalVector(), coeff);

        NeoN::Vector<NeoN::scalar> source(exec, nfT.internalVector().size(), 0.0);
        fvcc::SourceTerm sourceTerm(dsl::Operator::Type::Implicit, coefficients, nfT);
        sourceTerm.explicitOperation(source);

        auto sourceHost = source.copyToHost();
        const auto sourceView = sourceHost.view();
        auto nftHost = nfT.internalVector().copyToHost();
        const auto hostnfTView = nftHost.view();
        for (size_t i = 0; i < sourceHost.size(); i++)
        {
            REQUIRE(sourceView[i] == coeff * hostnfTView[i]);
        }

        // TODO finish the sourceterm operator
        // auto ls = sourceTerm.createEmptyLinearSystem();
        // sourceTerm.implicitOperation(ls);
        // fvcc::Expression<NeoN::scalar> ls2(
        //     nfT,
        //     ls,
        //     fvcc::SparsityPattern::readOrCreate(nfMesh)
        // );

        // // check diag
        // NeoN::Vector<NeoN::scalar> diag(nfT.exec(), nfT.internalVector().size(), 0.0);
        // ls2.diag(diag);
        // auto diagHost = diag.copyToHost();
        // auto cellVolumes = nfMesh.cellVolumes().copyToHost();
        // for (size_t i = 0; i < diagHost.size(); i++)
        // {
        //     REQUIRE(diagHost[i] == coeff * cellVolumes[i]);
        // }
        // auto result = ls2 & nfT;
        // auto implicitHost = result.internalVector().copyToHost();
        // for (size_t i = 0; i < implicitHost.size(); i++)
        // {
        //     REQUIRE(implicitHost[i] == coeff * nftHost[i] * cellVolumes[i]);
        // }
    }

    // SECTION("solve sourceterm_" + execName)
    // {
    //     auto ofT = randomScalarField(runTime, mesh, "T");
    //     ofT.primitiveFieldRef() = 1.0;
    //     ofT.correctBoundaryConditions();
    //
    //     fvcc::VolumeField<NeoN::scalar> nfT = constructFrom(exec, nfMesh, ofT);
    //     NeoN::fill(nfT.internalVector(), 1.0);
    //
    //     auto nfCoeff1 = nfT;
    //     NeoN::fill(nfCoeff1.internalVector(), 1.0);
    //
    //     auto nfCoeff2 = nfT;
    //     NeoN::fill(nfCoeff2.internalVector(), 2.0);
    //
    //     Foam::dimensionedScalar coeff1("coeff", Foam::dimless, 1.0);
    //     Foam::dimensionedScalar coeff2("coeff", Foam::dimless, 2.0);
    //     Foam::fvScalarMatrix matrix(Foam::fvm::Sp(coeff1, ofT) - coeff2 * ofT);
    //
    //
    //     dsl::Expression eqnSys(dsl::imp::source(nfCoeff1, nfT) - dsl::exp::source(nfCoeff2,
    //     nfT)); NeoN::scalar t = 0; NeoN::scalar dt = 1; NeoN::Dictionary fvSchemesDict {};
    //     NeoN::Dictionary fvSolutionDict {};
    //     fvSolutionDict.insert("maxIters", 100);
    //     fvSolutionDict.insert("relTol", float(1e-7));
    //
    //
    //     dsl::solve(eqnSys, nfT, t, dt, fvSchemesDict, fvSolutionDict);
    //     matrix.solve();
    //
    //     auto nfTHost = nfT.internalVector().copyToHost();
    //     const auto nfTHostView = nfTHost.view();
    //     for (size_t celli = 0; celli < nfTHost.size(); celli++)
    //     {
    //         REQUIRE(nfTHostView[celli] == Catch::Approx(ofT[celli]).margin(1e-16));
    //     }
    // }


    SECTION("solve div" + execName)
    {
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
        NeoN::scalar coeff = 1000;
        Foam::dimensionedScalar coeff1("coeff", Foam::dimensionSet(0, -3, 0, 0, 0), coeff);
        Foam::dimensionedScalar coeff2("coeff", Foam::dimensionSet(0, -3, 0, 0, 0), coeff);
        NeoN::fill(nfCoeff1.internalVector(), coeff);

        auto nfCoeff2 = nfT;
        NeoN::fill(nfCoeff2.internalVector(), coeff);

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
            dsl::imp::source(nfCoeff1, nfT) + dsl::imp::div(nfPhi, nfT)
            - dsl::exp::source(nfCoeff2, nfT)
        );

        NeoN::scalar t = 0;
        NeoN::scalar dt = 1;
        NeoN::Dictionary fvSchemesDict {};
        NeoN::Dictionary divSchemes {};
        divSchemes.insert(
            "div(phi,T)",
            NeoN::TokenList {std::string("Gauss"), std::string("linear")}
        );
        fvSchemesDict.insert("divSchemes", divSchemes);

        NeoN::Dictionary fvSolutionDict {
            {{"solver", std::string {"Ginkgo"}},
             {"type", "solver::Cg"},
             {"criteria", NeoN::Dictionary {{{"iteration", 100}, {"relative_residual_norm", 1e-8}}}}
            }
        };

        matrix.solve();
        dsl::solve(eqnSys, nfT, t, dt, fvSchemesDict, fvSolutionDict);

        std::span<Foam::scalar> ofTSpan(ofT.data(), ofT.size());
        auto nfTHost = nfT.internalVector().copyToHost();
        std::span<NeoN::scalar> nfTHostSpan(nfTHost.data(), nfTHost.size());
        // for (size_t celli = 0; celli < nfTHost.size(); celli++)
        // {
        //     REQUIRE(nfTHost.view()[celli] == Catch::Approx(ofT[celli]).margin(1e-16));
        // }
    }
}
