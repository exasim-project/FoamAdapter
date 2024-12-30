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
namespace fvm = Foam::fvm;

namespace dsl = NeoFOAM::dsl;
namespace fvcc = NeoFOAM::finiteVolume::cellCentred;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{

#include "addCheckCaseOptions.H"
#include "setRootCase.H"
#include "createTime.H"
    // #include "createMesh.H" requires using Foam which is incompatible with NVCC

    std::unique_ptr<Foam::fvMesh> meshPtr = Foam::createMesh(runTime);
    Foam::fvMesh& mesh = *meshPtr;

#include "createControl.H"

    auto [adjustTimeStep, maxCo, maxDeltaT] = Foam::timeControls(runTime);

#include "createFields.H"

    std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);

    Foam::scalar coNum = Foam::calculateCoNum(phi);
    if (adjustTimeStep)
    {
        Foam::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
    }

    Info << "\nStarting time loop\n" << endl;

    auto endTime = runTime.endTime().value();

    while (runTime.run())
    {
        std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);


        runTime++;
        Info << "Time = " << runTime.timeName() << endl;

        Foam::scalar t = runTime.time().value();
        Foam::scalar dt = runTime.deltaT().value();

        Foam::scalar pi = Foam::constant::mathematical::pi;
        phi = phi0 * Foam::cos(pi * (t + 0.5 * dt) / endTime);

        if (adjustTimeStep)
        {
            coNum = calculateCoNum(phi);
            Foam::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
        }

        // Info << "max(phi) : " << max(phi).value() << endl;
        // Info << "max(U) : " << max(U).value() << endl;

        // advance Foam fields in time
        Foam::fvScalarMatrix TEqn(fvm::ddt(T) + fvc::div(phi, T));

        TEqn.solve();

        if (runTime.outputTime())
        {
            U = U0 * Foam::cos(pi * (t + 0.5 * dt) / endTime);
        }


        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info << "End\n" << endl;


    return 0;
}

// ************************************************************************* //
