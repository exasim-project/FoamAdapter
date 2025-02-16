// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 NeoFOAM authors

#include "NeoFOAM/NeoFOAM.hpp"

#include "FoamAdapter/FoamAdapter.hpp"
#include "FoamAdapter/readers/foamDictionary.hpp"


#include "fvCFD.H"
#include "pisoControl.H"

using Foam::Info;
using Foam::endl;
using Foam::nl;
namespace fvc = Foam::fvc;
namespace fvm = Foam::fvm;

namespace dsl = NeoFOAM::dsl;
namespace fvcc = NeoFOAM::finiteVolume::cellCentred;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    Kokkos::initialize(argc, argv);
    {
#include "addCheckCaseOptions.H"
#include "setRootCase.H"
#include "createTime.H"

        NeoFOAM::Database db;

        fvcc::FieldCollection& fieldCollection =
            fvcc::FieldCollection::instance(db, "fieldCollection");


        NeoFOAM::Dictionary controlDict = Foam::readFoamDictionary(runTime.controlDict());
        NeoFOAM::Executor exec = createExecutor(runTime.controlDict());

        std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, runTime);
        Foam::MeshAdapter& mesh = *meshPtr;

        Foam::pisoControl piso(mesh);

        auto [adjustTimeStep, maxCo, maxDeltaT] = Foam::timeControls(runTime);

#include "createFields.H"


        std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);

        Foam::scalar coNum = Foam::calculateCoNum(phi);
        if (adjustTimeStep)
        {
            Foam::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
        }


        NeoFOAM::Dictionary fvSchemesDict = Foam::readFoamDictionary(mesh.schemesDict());
        NeoFOAM::Dictionary fvSolutionDict = Foam::readFoamDictionary(mesh.solutionDict());

        Info << "creating NeoFOAM mesh" << endl;
        NeoFOAM::UnstructuredMesh& nfMesh = mesh.nfMesh();

        Info << "creating NeoFOAM fields" << endl;
        fvcc::VolumeField<NeoFOAM::scalar>& nfp =
            fieldCollection.registerField<fvcc::VolumeField<NeoFOAM::scalar>>(
                Foam::CreateFromFoamField<Foam::volScalarField> {
                    .exec = exec,
                    .nfMesh = nfMesh,
                    .foamField = p,
                    .name = "nfp"
                }
            );


        Foam::scalar endTime = controlDict.get<Foam::scalar>("endTime");


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        Info << "\nStarting time loop\n" << endl;

        while (runTime.loop())
        {

            std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);
            coNum = calculateCoNum(phi);
            Foam::Info << "max(phi) : " << max(phi).value() << Foam::endl;
            Foam::Info << "max(U) : " << max(U).value() << Foam::endl;
            if (adjustTimeStep)
            {
                Foam::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
            }
            runTime++;

            Foam::scalar t = runTime.time().value();
            Foam::scalar dt = runTime.deltaT().value();

            Info << "Time = " << runTime.timeName() << nl << endl;

            // #include "CourantNo.H"

            // Momentum predictor

            Foam::fvVectorMatrix UEqn(fvm::ddt(U) + fvm::div(phi, U) - fvm::laplacian(nu, U));

            if (piso.momentumPredictor())
            {
                solve(UEqn == -fvc::grad(p));
            }

            // --- PISO loop
            while (piso.correct())
            {
                Foam::volScalarField rAU(1.0 / UEqn.A());
                Foam::volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));
                Foam::surfaceScalarField phiHbyA(
                    "phiHbyA",
                    fvc::flux(HbyA) + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)
                );

                Foam::adjustPhi(phiHbyA, U, p);

                // Update the pressure BCs to ensure flux consistency
                Foam::constrainPressure(p, U, phiHbyA, rAU);

                // Non-orthogonal pressure corrector loop
                while (piso.correctNonOrthogonal())
                {
                    // Pressure corrector
                    auto nfPhiHbyA = Foam::constructSurfaceField(exec, nfMesh, phiHbyA);
                    Foam::surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
                    auto nfrAUf = Foam::constructSurfaceField(exec, nfMesh, rAUf);
                    Foam::fvScalarMatrix pEqn(fvm::laplacian(rAUf, p) == fvc::div(phiHbyA));

                    // dsl::Expression pEqn2(dsl::imp::laplacian(nfrAUf,nfp) -
                    // dsl::exp::div(nfPhiHbyA));
                    dsl::Expression pEqn2(
                        dsl::imp::laplacian(nfrAUf, nfp) - dsl::exp::div(nfPhiHbyA, nfp)
                    );

                    dsl::solve(pEqn2, nfp, t, dt, fvSchemesDict, fvSolutionDict);

                    pEqn.setReference(pRefCell, pRefValue);

                    pEqn.solve(p.select(piso.finalInnerIter()));

                    if (piso.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                // #include "continuityErrs.H"

                U = HbyA - rAU * fvc::grad(p);
                U.correctBoundaryConditions();
            }

            runTime.write();

            runTime.printExecutionTime(Info);
        }

        Info << "End\n" << endl;
    }
    Kokkos::finalize();

    return 0;
}

// ************************************************************************* //
