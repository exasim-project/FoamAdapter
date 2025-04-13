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
        auto& solverDict = fvSolutionDict.get<NeoFOAM::Dictionary>("solvers");

        Info << "creating NeoFOAM mesh" << endl;
        NeoFOAM::UnstructuredMesh& nfMesh = mesh.nfMesh();

        Info << "creating NeoFOAM pressure fields" << endl;
        fvcc::VolumeField<NeoFOAM::scalar>& nfp =
            fieldCollection.registerField<fvcc::VolumeField<NeoFOAM::scalar>>(
                Foam::CreateFromFoamField<Foam::volScalarField> {
                    .exec = exec,
                    .nfMesh = nfMesh,
                    .foamField = p,
                    .name = "nfp"
                }
            );

        Info << "creating NeoFOAM velocity fields" << endl;
        fvcc::VolumeField<NeoFOAM::Vector>& nfU =
            fieldCollection.registerField<fvcc::VolumeField<NeoFOAM::Vector>>(
                Foam::CreateFromFoamField<Foam::volVectorField> {
                    .exec = exec,
                    .nfMesh = nfMesh,
                    .foamField = U,
                    .name = "nfU"
                }
            );

        auto nuBCs = fvcc::createCalculatedBCs<fvcc::SurfaceBoundary<NeoFOAM::scalar>>(nfMesh);
        fvcc::SurfaceField<NeoFOAM::scalar> nfNu(exec, "nfNu", nfMesh, nuBCs);
        fill(nfNu.internalField(), nu.value());
        fill(nfNu.boundaryField().value(), nu.value());

        NeoFOAM::scalar endTime = controlDict.get<NeoFOAM::scalar>("endTime");

        Foam::surfaceScalarField flux("flux", phi);
        auto nfPhi = Foam::constructSurfaceField(exec, nfMesh, phi);

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        Info << "\nStarting time loop\n" << endl;

        while (runTime.loop())
        {
            auto& nfOldU = fvcc::oldTime(nfU);
            nfOldU.internalField() = nfU.internalField();
            std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);
            coNum = calculateCoNum(phi);
            Foam::Info << "max(phi) : " << max(phi).value() << Foam::endl;
            Foam::Info << "max(U) : " << max(U).value() << Foam::endl;
            if (adjustTimeStep)
            {
                Foam::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
            }

            Foam::scalar t = runTime.time().value();
            Foam::scalar dt = runTime.deltaT().value();

            Info << "Time = " << runTime.timeName() << nl << endl;

            // Momentum predictor
            fvcc::Expression<NeoFOAM::Vector> UEqn2(
                dsl::imp::ddt(nfU) + dsl::imp::div(nfPhi, nfU) - dsl::imp::laplacian(nfNu, nfU),
                nfU,
                fvSchemesDict,
                solverDict.get<NeoFOAM::Dictionary>("nfU")
            );

            UEqn2.assemble(t, dt);

            // --- PISO loop
            while (piso.correct())
            {
                Info << "PISO loop" << endl;
                auto [nfrAU, nfHbyA] = fvcc::discreteMomentumFields(UEqn2);
                fvcc::constrainHbyA(nfHbyA, nfU, nfp);

                fvcc::SurfaceField<NeoFOAM::scalar> nfrAUf =
                    fvcc::SurfaceInterpolation<NeoFOAM::scalar>(
                        exec,
                        nfMesh,
                        NeoFOAM::TokenList({std::string("linear")})
                    )
                        .interpolate(nfrAU);
                nfrAUf.name = "nfrAUf";

                // TODO: + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)
                auto nfPhiHbyA = fvcc::flux(nfHbyA);

                // TODO:
                // Foam::adjustPhi(phiHbyA, U, p);

                // Update the pressure BCs to ensure flux consistency
                // Foam::constrainPressure(p, U, phiHbyA, rAU);

                // Non-orthogonal pressure corrector loop
                while (piso.correctNonOrthogonal())
                {
                    // Pressure corrector

                    fvcc::Expression<NeoFOAM::scalar> pEqn2(
                        dsl::imp::laplacian(nfrAUf, nfp) - dsl::exp::div(nfPhiHbyA),
                        nfp,
                        fvSchemesDict,
                        solverDict.get<NeoFOAM::Dictionary>("nfP")
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
                        fvcc::updateFaceVelocity(nfPhi, nfPhiHbyA, pEqn2);
                    }
                }
                // TODO:
                // #include "continuityErrs.H"

                fvcc::updateVelocity(nfU, nfHbyA, nfrAU, nfp);
                nfU.correctBoundaryConditions();
            }

            runTime.write();
            if (runTime.outputTime())
            {
                Info << "writing nfp field" << endl;
                write(nfp.internalField(), mesh, "nfp");
                write(nfU.internalField(), mesh, "nfU");
            }

            runTime.printExecutionTime(Info);
        }

        Info << "End\n" << endl;
    }
    Kokkos::finalize();

    return 0;
}

// ************************************************************************* //
