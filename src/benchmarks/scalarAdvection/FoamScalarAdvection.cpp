// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"

#include "FoamAdapter/setup/setup.hpp"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{

#include "addProfilingOption.H"
#include "addCheckCaseOptions.H"
#include "setRootCase.H"
#include "createTime.H"
    std::unique_ptr<Foam::fvMesh> meshPtr = Foam::createMesh(runTime);
    Foam::fvMesh& mesh = *meshPtr;
#include "createControl.H"
    // #include "createTimeControls.H"
    auto [adjustTimeStep, maxCo, maxDeltaT] = Foam::timeControls(runTime);


#include "createFields.H"
    // set temperature
    Foam::scalar spread = 0.05;
    forAll(T, celli)
    {
        T[celli] = std::exp(
            -0.5
            * (std::pow((mesh.C()[celli].x() - 0.5) / spread, 2.0)
               + std::pow((mesh.C()[celli].y() - 0.75) / spread, 2.0))
        );
    }
    T.correctBoundaryConditions();
    T.write();

    std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);

    Foam::scalar CoNum = Foam::calculateCoNum(phi);
    if (adjustTimeStep)
    {
        Foam::setDeltaT(runTime, maxCo, CoNum, maxDeltaT);
    }

    Foam::scalar pi = Foam::constant::mathematical::pi;
    {
        Foam::scalarField X(mesh.C().component(0));
        Foam::scalarField Y(mesh.C().component(1));
        Foam::scalarField u(-Foam::sin(2.0 * pi * Y) * Foam::pow(Foam::sin(pi * X), 2.0));
        Foam::scalarField w(Foam::sin(2.0 * pi * X) * Foam::pow(Foam::sin(pi * Y), 2.0));
        forAll(U0, celli)
        {
            U0[celli].x() = u[celli];
            U0[celli].y() = w[celli];
            U0[celli].z() = 0.0;
        }
    }
    phi0 = Foam::linearInterpolate(U0) & mesh.Sf();


    while (runTime.run())
    {
        // #include "readTimeControls.H"
        std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);
        CoNum = calculateCoNum(phi);
        Foam::Info << "max(phi) : " << max(phi) << Foam::endl;
        Foam::Info << "max(U) : " << max(U) << Foam::endl;
        if (adjustTimeStep)
        {
            Foam::setDeltaT(runTime, maxCo, CoNum, maxDeltaT);
        }

        runTime++;

        Foam::Info << "Time = " << runTime.timeName() << Foam::nl << Foam::endl;

        if (spirallingFlow > 0)
        {
            Foam::Info << "Spiralling flow: " << spirallingFlow << Foam::endl;
            Foam::scalar t = runTime.time().value();
            Foam::scalar dt = runTime.deltaT().value();
            U = U0 * Foam::cos(pi * (t + 0.5 * dt) / spirallingFlow);
            phi = phi0 * Foam::cos(pi * (t + 0.5 * dt) / spirallingFlow);
        }

        {
            addProfiling(foamAdvection, "foamAdvection");
            Foam::fvScalarMatrix TEqn(Foam::fvm::ddt(T) + Foam::fvc::div(phi, T));

            TEqn.solve();
        }

        runTime.write();

        runTime.printExecutionTime(Foam::Info);
    }

    Foam::Info << "End\n" << Foam::endl;

    return 0;
}

// ************************************************************************* //
