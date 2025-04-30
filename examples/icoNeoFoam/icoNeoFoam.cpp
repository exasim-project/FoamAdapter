// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 FoamAdapter authors

#include "NeoN/NeoN.hpp"

#include "FoamAdapter/FoamAdapter.hpp"


#include "fvCFD.H"
#include "pisoControl.H"

using Foam::Info;
using Foam::endl;
using Foam::nl;
namespace fvc = Foam::fvc;
namespace fvm = Foam::fvm;

namespace dsl = NeoN::dsl;
namespace nnfvcc = NeoN::finiteVolume::cellCentred;
namespace nffvcc = FoamAdapter;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    Kokkos::initialize(argc, argv);
    {
#include "addCheckCaseOptions.H"
#include "setRootCase.H"
#include "createTime.H"

        NeoN::Database db;

        fvcc::VectorCollection& vectorCollection =
            fvcc::VectorCollection::instance(db, "VectorCollection");


        NeoN::Dictionary controlDict = FoamAdapter::convert(runTime.controlDict());
        NeoN::Executor exec = FoamAdapter::createExecutor(runTime.controlDict());

        std::unique_ptr<FoamAdapter::MeshAdapter> meshPtr = FoamAdapter::createMesh(exec, runTime);
        FoamAdapter::MeshAdapter& mesh = *meshPtr;

        Foam::pisoControl piso(mesh);


#include "createFields.H"

        auto [adjustTimeStep, maxCo, maxDeltaT] = FoamAdapter::timeControls(runTime);

        // Foam::scalar coNum = Foam::calculateCoNum(phi);
        // if (adjustTimeStep)
        // {
        //     Foam::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
        // }


        NeoN::Dictionary fvSchemesDict = FoamAdapter::convert(mesh.schemesDict());
        NeoN::Dictionary fvSolutionDict = FoamAdapter::convert(mesh.solutionDict());
        auto& solverDict = fvSolutionDict.get<NeoN::Dictionary>("solvers");

        Info << "creating FoamAdapter mesh" << endl;
        NeoN::UnstructuredMesh& nfMesh = mesh.nfMesh();

        Info << "creating FoamAdapter pressure fields" << endl;
        fvcc::VolumeField<NeoN::scalar>& nfp =
            vectorCollection.registerVector<fvcc::VolumeField<NeoN::scalar>>(
                FoamAdapter::CreateFromFoamField<Foam::volScalarField> {
                    .exec = exec,
                    .nfMesh = nfMesh,
                    .foamField = p,
                    .name = "nfp"
                }
            );

        Info << "creating FoamAdapter velocity fields" << endl;
        fvcc::VolumeField<NeoN::Vec3>& nfU =
            vectorCollection.registerVector<fvcc::VolumeField<NeoN::Vec3>>(
                FoamAdapter::CreateFromFoamField<Foam::volVectorField> {
                    .exec = exec,
                    .nfMesh = nfMesh,
                    .foamField = U,
                    .name = "nfU"
                }
            );

        auto nuBCs = fvcc::createCalculatedBCs<fvcc::SurfaceBoundary<NeoN::scalar>>(nfMesh);
        fvcc::SurfaceField<NeoN::scalar> nfNu(exec, "nfNu", nfMesh, nuBCs);
        NeoN::fill(nfNu.internalVector(), nu.value());
        NeoN::fill(nfNu.boundaryData().value(), nu.value());

        NeoN::scalar endTime = controlDict.get<NeoN::scalar>("endTime");

        Foam::surfaceScalarField flux("flux", phi);
        auto nfPhi = FoamAdapter::constructSurfaceField(exec, nfMesh, phi);

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        Info << "\nStarting time loop\n" << endl;

        while (runTime.loop())
        {
            Foam::scalar t = runTime.time().value();
            Foam::scalar dt = runTime.deltaT().value();

            auto& nfOldU = fvcc::oldTime(nfU);
            nfOldU.internalVector() = nfU.internalVector();
            auto [adjustTimeStep, maxCo, maxDeltaT] = FoamAdapter::timeControls(runTime);
            auto coNum = fvcc::computeCoNum(nfPhi, dt);
            Foam::Info << "max(phi) : " << max(phi).value() << Foam::endl;
            Foam::Info << "max(U) : " << max(U).value() << Foam::endl;
            if (adjustTimeStep)
            {
                FoamAdapter::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
            }

            Foam::Info << "Time = " << runTime.timeName() << Foam::nl << Foam::endl;

            // Momentum predictor
            nffvcc::Expression<NeoN::Vec3> UEqn2(
                NeoN::dsl::imp::ddt(nfU) + NeoN::dsl::imp::div(nfPhi, nfU)
                    - NeoN::dsl::imp::laplacian(nfNu, nfU),
                nfU,
                fvSchemesDict,
                solverDict.get<NeoN::Dictionary>("nfU")
            );

            UEqn2.assemble(t, dt);

            // --- PISO loop
            while (piso.correct())
            {
                Info << "PISO loop" << endl;
                auto [nfrAU, nfHbyA] = nffvcc::discreteMomentumFields(UEqn2);
                nffvcc::constrainHbyA(nfHbyA, nfU, nfp);

                nnfvcc::SurfaceField<NeoN::scalar> nfrAUf =
                    nnfvcc::SurfaceInterpolation<NeoN::scalar>(
                        exec,
                        nfMesh,
                        NeoN::TokenList({std::string("linear")})
                    )
                        .interpolate(nfrAU);
                nfrAUf.name = "nfrAUf";

                // TODO: + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)
                auto nfPhiHbyA = nffvcc::flux(nfHbyA);

                // TODO:
                // Foam::adjustPhi(phiHbyA, U, p);

                // Update the pressure BCs to ensure flux consistency
                // Foam::constrainPressure(p, U, phiHbyA, rAU);

                // Non-orthogonal pressure corrector loop
                while (piso.correctNonOrthogonal())
                {
                    // Pressure corrector

                    nffvcc::Expression<NeoN::scalar> pEqn2(
                        NeoN::dsl::imp::laplacian(nfrAUf, nfp) - NeoN::dsl::exp::div(nfPhiHbyA),
                        nfp,
                        fvSchemesDict,
                        solverDict.get<NeoN::Dictionary>("nfP")
                    );

                    pEqn2.assemble(t, dt);

                    if (p.needReference() && pRefCell >= 0)
                    {
                        pEqn2.setReference(pRefCell, pRefValue);
                    }
                    pEqn2.solve(t, dt);
                    nfp.correctBoundaryConditions();

                    if (piso.finalNonOrthogonalIter())
                    {
                        nffvcc::updateFaceVelocity(nfPhi, nfPhiHbyA, pEqn2);
                    }
                }
                // TODO:
                // #include "continuityErrs.H"

                nffvcc::updateVelocity(nfU, nfHbyA, nfrAU, nfp);
                nfU.correctBoundaryConditions();
            }

            runTime.write();
            if (runTime.outputTime())
            {
                Info << "writing nfp field" << endl;
                write(nfp.internalVector(), mesh, "nfp");
                write(nfU.internalVector(), mesh, "nfU");
            }

            runTime.printExecutionTime(Info);
        }

        Info << "End\n" << endl;
    }
    Kokkos::finalize();

    return 0;
}

// ************************************************************************* //
