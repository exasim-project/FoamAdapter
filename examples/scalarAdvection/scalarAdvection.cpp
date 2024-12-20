// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "FoamAdapter/FoamAdapter.hpp"
#include "NeoFOAM/dsl/expression.hpp"
#include "NeoFOAM/dsl/solver.hpp"
#include "NeoFOAM/dsl/ddt.hpp"
#include "FoamAdapter/readers/foamDictionary.hpp"

#include "NeoFOAM/dsl/implicit.hpp"
#include "NeoFOAM/dsl/explicit.hpp"


#define namespaceFoam
#include "fvCFD.H"

using Foam::Info;
using Foam::endl;
using Foam::nl;
namespace fvc = Foam::fvc;

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

        NeoFOAM::Dictionary controlDict = Foam::readFoamDictionary(runTime.controlDict());
        NeoFOAM::Executor exec = createExecutor(runTime.controlDict());

        std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, runTime);
        Foam::MeshAdapter& mesh = *meshPtr;

#include "createControl.H"

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
        auto nfT = Foam::constructFrom(exec, nfMesh, T);
        auto nfPhi0 = Foam::constructSurfaceField(exec, nfMesh, phi0);
        auto nfPhi = Foam::constructSurfaceField(exec, nfMesh, phi);

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

                nfPhi.internalField() =
                    nfPhi0.internalField() * std::cos(pi * (t + 0.5 * dt) / endTime);
            }


            std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);
            coNum = calculateCoNum(phi);
            Foam::Info << "max(phi) : " << max(phi).value() << Foam::endl;
            Foam::Info << "max(U) : " << max(U).value() << Foam::endl;
            if (adjustTimeStep)
            {
                Foam::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
            }
            runTime++;

            Info << "Time = " << runTime.timeName() << endl;


            {
                dsl::Expression eqnSys(dsl::imp::ddt(nfT) + dsl::exp::div(nfPhi, nfT));

                dsl::solve(eqnSys, nfT, t, dt, fvSchemesDict, fvSolutionDict);
            }

            if (runTime.outputTime())
            {
                Info << "writing nfT field" << endl;
                write(nfT.internalField(), mesh, "nfT");
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
