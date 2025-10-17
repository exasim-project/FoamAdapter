// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023-2025 FoamAdapter authors

#include "NeoN/NeoN.hpp"

#include "FoamAdapter/FoamAdapter.hpp"

#include "fvCFD.H"

using Foam::Info;
using Foam::endl;
using Foam::nl;

namespace dsl = NeoN::dsl;
namespace fvcc = NeoN::finiteVolume::cellCentred;
namespace nf = FoamAdapter;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    Kokkos::initialize(argc, argv);
    {
#include "addCheckCaseOptions.H"
#include "setRootCase.H"
#include "createTime.H"

        auto rt = nf::createAdapterRunTime(runTime);
        auto& mesh = rt.mesh;

#include "createControl.H"
#include "createFields.H"

        Info << "creating FoamAdapter fields" << endl;
        fvcc::VectorCollection& VectorCollection =
            fvcc::VectorCollection::instance(rt.db, "VectorCollection");
        fvcc::VolumeField<NeoN::scalar>& nfT =
            VectorCollection.registerVector<fvcc::VolumeField<NeoN::scalar>>(
                FoamAdapter::CreateFromFoamField<Foam::volScalarField> {
                    .exec = rt.exec,
                    .nfMesh = rt.nfMesh,
                    .foamField = T,
                    .name = "nfT"
                }
            );
        auto nfPhi0 = FoamAdapter::constructSurfaceField(rt.exec, rt.nfMesh, phi0);
        auto nfPhi = FoamAdapter::constructSurfaceField(rt.exec, rt.nfMesh, phi);

        Foam::scalar endTime = rt.controlDict.get<Foam::scalar>("endTime");

        while (runTime.run())
        {
            Foam::scalar t = runTime.time().value();
            Foam::scalar dt = runTime.deltaT().value();

            auto& nfOldT = fvcc::oldTime(nfT);
            nfOldT.internalVector() = nfT.internalVector();

            if (rt.controlDict.get<int>("setFields"))
            {
                Foam::scalar pi = Foam::constant::mathematical::pi;
                U = U0 * Foam::cos(pi * (t + 0.5 * dt) / endTime);
                phi = phi0 * Foam::cos(pi * (t + 0.5 * dt) / endTime);

                nfPhi.internalVector() =
                    nfPhi0.internalVector() * std::cos(pi * (t + 0.5 * dt) / endTime);
            }

            auto coNum = fvcc::computeCoNum(nfPhi, dt);
            Foam::Info << "max(phi) : " << max(phi).value() << Foam::endl;
            Foam::Info << "max(U) : " << max(U).value() << Foam::endl;
            if (rt.adjustTimeStep)
            {
                FoamAdapter::setDeltaT(runTime, rt, coNum);
            }
            runTime++;

            Info << "Time = " << runTime.timeName() << endl;
            {
                dsl::Expression eqnSys(dsl::imp::ddt(nfT) + dsl::imp::div(nfPhi, nfT));

                dsl::solve(eqnSys, nfT, t, dt, rt.fvSchemesDict, rt.fvSolutionDict);
            }

            if (runTime.outputTime())
            {
                Foam::Info << "writing nfT field" << Foam::endl;
                FoamAdapter::write(nfT.internalVector(), rt.mesh, "nfT");
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
