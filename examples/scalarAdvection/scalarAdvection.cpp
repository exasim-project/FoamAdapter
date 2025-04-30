// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 FoamAdapter authors

#include "NeoN/NeoN.hpp"

#include "FoamAdapter/FoamAdapter.hpp"

#include "fvCFD.H"

using Foam::Info;
using Foam::endl;
using Foam::nl;
namespace fvc = Foam::fvc;

namespace dsl = NeoN::dsl;
namespace fvcc = NeoN::finiteVolume::cellCentred;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    Kokkos::initialize(argc, argv);
    {
#include "addCheckCaseOptions.H"
#include "setRootCase.H"
#include "createTime.H"

        NeoN::Database db;

        fvcc::VectorCollection& VectorCollection =
            fvcc::VectorCollection::instance(db, "VectorCollection");

        NeoN::Dictionary controlDict = FoamAdapter::convert(runTime.controlDict());
        NeoN::Executor exec = FoamAdapter::createExecutor(runTime.controlDict());

        std::unique_ptr<FoamAdapter::MeshAdapter> meshPtr = FoamAdapter::createMesh(exec, runTime);
        FoamAdapter::MeshAdapter& mesh = *meshPtr;

#include "createControl.H"

        auto [adjustTimeStep, maxCo, maxDeltaT] = FoamAdapter::timeControls(runTime);

#include "createFields.H"

        // Foam::scalar coNum = fvcc::computeCoNum(phi);
        // if (adjustTimeStep)
        // {
        //     FoamAdapter::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
        // }


        NeoN::Dictionary fvSchemesDict = FoamAdapter::convert(mesh.schemesDict());
        NeoN::Dictionary fvSolutionDict = FoamAdapter::convert(mesh.solutionDict());

        Info << "creating FoamAdapter mesh" << endl;
        NeoN::UnstructuredMesh& nfMesh = mesh.nfMesh();

        Info << "creating FoamAdapter fields" << endl;
        fvcc::VolumeField<NeoN::scalar>& nfT =
            VectorCollection.registerVector<fvcc::VolumeField<NeoN::scalar>>(
                FoamAdapter::CreateFromFoamField<Foam::volScalarField> {
                    .exec = exec,
                    .nfMesh = nfMesh,
                    .foamField = T,
                    .name = "nfT"
                }
            );
        auto nfPhi0 = FoamAdapter::constructSurfaceField(exec, nfMesh, phi0);
        auto nfPhi = FoamAdapter::constructSurfaceField(exec, nfMesh, phi);

        Foam::scalar endTime = controlDict.get<Foam::scalar>("endTime");

        while (runTime.run())
        {
            Foam::scalar t = runTime.time().value();
            Foam::scalar dt = runTime.deltaT().value();

            if (controlDict.get<int>("setFields"))
            {
                Foam::scalar pi = Foam::constant::mathematical::pi;
                U = U0 * Foam::cos(pi * (t + 0.5 * dt) / endTime);
                phi = phi0 * Foam::cos(pi * (t + 0.5 * dt) / endTime);

                nfPhi.internalVector() =
                    nfPhi0.internalVector() * std::cos(pi * (t + 0.5 * dt) / endTime);
            }


            std::tie(adjustTimeStep, maxCo, maxDeltaT) = FoamAdapter::timeControls(runTime);
            auto coNum = fvcc::computeCoNum(nfPhi, dt);
            Foam::Info << "max(phi) : " << max(phi).value() << Foam::endl;
            Foam::Info << "max(U) : " << max(U).value() << Foam::endl;
            if (adjustTimeStep)
            {
                FoamAdapter::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
            }
            runTime++;

            Info << "Time = " << runTime.timeName() << endl;

            {
                NeoN::dsl::Expression eqnSys(
                    NeoN::dsl::imp::ddt(nfT) + NeoN::dsl::imp::div(nfPhi, nfT)
                );

                NeoN::dsl::solve(eqnSys, nfT, t, dt, fvSchemesDict, fvSolutionDict);
            }

            if (runTime.outputTime())
            {
                Foam::Info << "writing nfT field" << Foam::endl;
                FoamAdapter::write(nfT.internalVector(), mesh, "nfT");
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
