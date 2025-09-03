// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 FoamAdapter authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "common.hpp"

using Foam::Info;
using Foam::endl;
using Foam::nl;
namespace fvc = Foam::fvc;
namespace fvm = Foam::fvm;

namespace dsl = NeoN::dsl;
namespace fvcc = NeoN::finiteVolume::cellCentred;

extern Foam::Time* timePtr; // A single time object


void initFields(Foam::volScalarField& T, Foam::volVectorField& U, Foam::surfaceScalarField& phi)
{

    const Foam::fvMesh& mesh = T.mesh();


    Foam::scalar spread = 0.05;
    Foam::scalar pi = Foam::constant::mathematical::pi;
    forAll(U, celli)
    {
        // initialize U
        Foam::scalar x = mesh.C()[celli].x();
        Foam::scalar y = mesh.C()[celli].y();

        U[celli].x() = -Foam::sin(2.0 * pi * y) * Foam::pow(Foam::sin(pi * x), 2.0);
        U[celli].y() = Foam::sin(2.0 * pi * x) * Foam::pow(Foam::sin(pi * y), 2.0);
        U[celli].z() = 0.0;

        // initialize T
        T[celli] = std::exp(
            -0.5
            * (std::pow((mesh.C()[celli].x() - 0.5) / spread, 2.0)
               + std::pow((mesh.C()[celli].y() - 0.75) / spread, 2.0))
        );
    }
    T.correctBoundaryConditions();

    phi = Foam::linearInterpolate(U) & mesh.Sf();
}

TEST_CASE("Advection Equation")
{
    Foam::Time& runTime = *timePtr;

    NeoN::Database db;
    fvcc::VectorCollection& vectorCollection =
        fvcc::VectorCollection::instance(db, "VectorCollection");

    auto [execName, exec] = GENERATE(allAvailableExecutor());

    Foam::scalar startTime = 0.0;
    Foam::label startTimeIndex = 0;
    runTime.setTime(startTime, startTimeIndex);

    // create mesh
    std::unique_ptr<FoamAdapter::MeshAdapter> meshAdapterPtr =
        FoamAdapter::createMesh(exec, runTime);
    FoamAdapter::MeshAdapter& mesh = *meshAdapterPtr;
    NeoN::UnstructuredMesh& nfMesh = mesh.nfMesh();

    Foam::volScalarField T(
        Foam::IOobject(
            "T",
            runTime.timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::AUTO_WRITE
        ),
        mesh
    );

    Foam::volVectorField U(
        Foam::IOobject(
            "U",
            runTime.timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::AUTO_WRITE
        ),
        mesh
    );

    Foam::surfaceScalarField phi("phi", Foam::linearInterpolate(U) & mesh.Sf());

    initFields(T, U, phi);

    // Copies of initial U and phi for use when flow is periodic
    Foam::surfaceScalarField phi0 = phi;
    Foam::volVectorField U0 = U;

    fvcc::VolumeField<NeoN::scalar>& nfT =
        vectorCollection.registerVector<fvcc::VolumeField<NeoN::scalar>>(
            FoamAdapter::CreateFromFoamField<Foam::volScalarField> {
                .exec = exec,
                .nfMesh = nfMesh,
                .foamField = T,
                .name = "nfT"
            }
        );
    auto nfPhi0 = FoamAdapter::constructSurfaceField(exec, nfMesh, phi0);
    auto nfPhi = FoamAdapter::constructSurfaceField(exec, nfMesh, phi);

    NeoN::Dictionary controlDict = FoamAdapter::convert(runTime.controlDict());
    NeoN::Dictionary fvSchemesDict = FoamAdapter::convert(mesh.schemesDict());
    Foam::scalar endTime = controlDict.get<Foam::scalar>("endTime");


    SECTION("Scalar advection with " + execName + " and " + "forwardEuler")
    {
        std::string timeIntegration = "forwardEuler";
        fvSchemesDict.get<NeoN::Dictionary>("ddtSchemes").insert("type", timeIntegration);

        while (runTime.run())
        {
            NeoN::Dictionary fvSolutionDict {
                {{"solver", std::string {"Ginkgo"}},
                 {"type", "solver::Bicgstab"},
                 {"criteria",
                  NeoN::Dictionary {{{"iteration", 500}, {"relative_residual_norm", 1e-14}}}}}
            };

            Foam::scalar t = runTime.time().value();
            Foam::scalar dt = runTime.deltaT().value();

            auto& nfOldT = fvcc::oldTime(nfT);
            nfOldT.internalVector() = nfT.internalVector();

            Foam::scalar pi = Foam::constant::mathematical::pi;
            U = U0 * Foam::cos(pi * (t + 0.5 * dt) / endTime);
            phi = phi0 * Foam::cos(pi * (t + 0.5 * dt) / endTime);

            nfPhi.internalVector() =
                nfPhi0.internalVector() * std::cos(pi * (t + 0.5 * dt) / endTime);

            runTime++;

            // advance Foam fields in time
            {
                Foam::fvScalarMatrix TEqn(fvm::ddt(T) + fvc::div(phi, T));
                TEqn.solve();
            }

            // advance FoamAdapter fields in time
            {
                NeoN::dsl::Expression eqnSys(
                    NeoN::dsl::imp::ddt(nfT) + NeoN::dsl::exp::div(nfPhi, nfT)
                );
                NeoN::dsl::solve(eqnSys, nfT, t, dt, fvSchemesDict, fvSolutionDict);
            }

            // for debugging with paraview
            if (runTime.outputTime())
            {
                write(nfT.internalVector(), mesh, "nfTExp_" + execName);
                T.write(); // for some reason T was not written
            }

            // runTime.write();
            // runTime.printExecutionTime(Info);
        }
        FoamAdapter::compare(nfT, T, ApproxScalar(1e-10), false);
    }

    std::string timeIntegration = "backwardEuler";
    NeoN::Dictionary fvSolutionDict {
        {{"solver", std::string {"Ginkgo"}},
         {"type", "solver::Bicgstab"},
         {"preconditioner",
          NeoN::Dictionary {{{"type", "preconditioner::Jacobi"}, {"max_block_size", 1}}}},
         {"criteria", NeoN::Dictionary {{{"iteration", 20}, {"relative_residual_norm", 1e-14}}}}}
    };

    SECTION("Scalar advection with " + execName + " and " + timeIntegration)
    {
        fvSchemesDict.get<NeoN::Dictionary>("ddtSchemes").insert("type", timeIntegration);

        while (runTime.run())
        {
            Foam::scalar t = runTime.time().value();
            Foam::scalar dt = runTime.deltaT().value();
            runTime++;

            auto& nfOldT = fvcc::oldTime(nfT);
            nfOldT.internalVector() = nfT.internalVector();

            Foam::scalar pi = Foam::constant::mathematical::pi;
            U = U0 * Foam::cos(pi * (t + 0.5 * dt) / endTime);
            phi = phi0 * Foam::cos(pi * (t + 0.5 * dt) / endTime);

            nfPhi.internalVector() =
                nfPhi0.internalVector() * std::cos(pi * (t + 0.5 * dt) / endTime);

            FoamAdapter::compare(nfT, T, ApproxScalar(1e-04), false);

            // advance Foam fields in time
            {
                Foam::fvScalarMatrix TEqn(fvm::ddt(T) + fvm::div(phi, T));
                TEqn.solve();
            }

            // advance FoamAdapter fields in time
            {
                NeoN::dsl::Expression eqnSys(
                    NeoN::dsl::imp::ddt(nfT) + NeoN::dsl::imp::div(nfPhi, nfT)
                );
                NeoN::dsl::solve(eqnSys, nfT, t, dt, fvSchemesDict, fvSolutionDict);
            }

            // for debugging with paraview
            if (runTime.outputTime())
            {
                write(nfT.internalVector(), mesh, "nfTImp_" + execName);
                T.write(); // for some reason T was not written
            }

            // runTime.write();
            // runTime.printExecutionTime(Info);
        }

        FoamAdapter::compare(nfT, T, ApproxScalar(1e-8), false);
    }
}
