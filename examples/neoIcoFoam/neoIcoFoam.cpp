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

        NeoN::Dictionary fvSchemesDict = FoamAdapter::convert(mesh.schemesDict());
        NeoN::Dictionary fvSolutionDict = FoamAdapter::convert(mesh.solutionDict());
        auto& solverDict = fvSolutionDict.get<NeoN::Dictionary>("solvers");

        Info << "creating FoamAdapter mesh" << endl;
        NeoN::UnstructuredMesh& nfMesh = mesh.nfMesh();

        Info << "creating FoamAdapter pressure fields" << endl;
        fvcc::VolumeField<NeoN::scalar>& p =
            vectorCollection.registerVector<fvcc::VolumeField<NeoN::scalar>>(
                FoamAdapter::CreateFromFoamField<Foam::volScalarField> {
                    .exec = exec,
                    .nfMesh = nfMesh,
                    .foamField = ofp,
                    .name = "p"
                }
            );

        Info << "creating FoamAdapter velocity fields" << endl;
        fvcc::VolumeField<NeoN::Vec3>& U =
            vectorCollection.registerVector<fvcc::VolumeField<NeoN::Vec3>>(
                FoamAdapter::CreateFromFoamField<Foam::volVectorField> {
                    .exec = exec,
                    .nfMesh = nfMesh,
                    .foamField = ofU,
                    .name = "U"
                }
            );

        auto nuBCs = fvcc::createCalculatedBCs<fvcc::SurfaceBoundary<NeoN::scalar>>(nfMesh);
        fvcc::SurfaceField<NeoN::scalar> nu(exec, "nu", nfMesh, nuBCs);
        NeoN::fill(nu.internalVector(), viscosity.value());
        NeoN::fill(nu.boundaryData().value(), viscosity.value());

        NeoN::scalar endTime = controlDict.get<NeoN::scalar>("endTime");

        auto phi = FoamAdapter::constructSurfaceField(exec, nfMesh, ofphi);

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        Info << "\nStarting time loop\n" << endl;

        while (runTime.loop())
        {
            Foam::scalar t = runTime.time().value();
            Foam::scalar dt = runTime.deltaT().value();

            auto& oldU = fvcc::oldTime(U);
            oldU.internalVector() = U.internalVector();
            auto [adjustTimeStep, maxCo, maxDeltaT] = FoamAdapter::timeControls(runTime);
            auto coNum = fvcc::computeCoNum(phi, dt);
            // Foam::Info << "max(phi) : " << max(phi).value() << Foam::endl;
            // Foam::Info << "max(U) : " << max(U).value() << Foam::endl;
            if (adjustTimeStep)
            {
                FoamAdapter::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
            }

            Foam::Info << "Time = " << runTime.timeName() << Foam::nl << Foam::endl;

            // Momentum predictor
            nffvcc::Expression<NeoN::Vec3> UEqn(
                NeoN::dsl::imp::ddt(U) + NeoN::dsl::imp::div(phi, U)
                    - NeoN::dsl::imp::laplacian(nu, U),
                U,
                fvSchemesDict,
                solverDict.get<NeoN::Dictionary>("U")
            );

            UEqn.assemble(t, dt);

            // --- PISO loop
            while (piso.correct())
            {
                Info << "PISO loop" << endl;
                auto [rAU, HbyA] = nffvcc::discreteMomentumFields(UEqn);
                nffvcc::constrainHbyA(HbyA, U, p);

                nnfvcc::SurfaceField<NeoN::scalar> nfrAUf =
                    nnfvcc::SurfaceInterpolation<NeoN::scalar>(
                        exec,
                        nfMesh,
                        NeoN::TokenList({std::string("linear")})
                    )
                        .interpolate(rAU);
                nfrAUf.name = "rAUf";

                // TODO: + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)
                auto phiHbyA = nffvcc::flux(HbyA);

                // TODO:
                // Foam::adjustPhi(phiHbyA, U, p);

                // Update the pressure BCs to ensure flux consistency
                // Foam::constrainPressure(p, U, phiHbyA, rAU);

                // Non-orthogonal pressure corrector loop
                while (piso.correctNonOrthogonal())
                {
                    // Pressure corrector

                    nffvcc::Expression<NeoN::scalar> pEqn(
                        NeoN::dsl::imp::laplacian(nfrAUf, p) - NeoN::dsl::exp::div(phiHbyA),
                        p,
                        fvSchemesDict,
                        solverDict.get<NeoN::Dictionary>("nfP")
                    );

                    pEqn.assemble(t, dt);

                    if (ofp.needReference() && pRefCell >= 0)
                    {
                        pEqn.setReference(pRefCell, pRefValue);
                    }
                    pEqn.solve(t, dt);
                    p.correctBoundaryConditions();

                    if (piso.finalNonOrthogonalIter())
                    {
                        nffvcc::updateFaceVelocity(phi, phiHbyA, pEqn);
                    }
                }
                // TODO:
                // #include "continuityErrs.H"

                nffvcc::updateVelocity(U, HbyA, rAU, p);
                U.correctBoundaryConditions();
            }

            runTime.write();
            if (runTime.outputTime())
            {
                Info << "writing p field" << endl;
                write(p.internalVector(), mesh, "p");
                write(U.internalVector(), mesh, "U");
            }

            runTime.printExecutionTime(Info);
        }

        Info << "End\n" << endl;
    }
    Kokkos::finalize();

    return 0;
}

// ************************************************************************* //
