/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    gradOperator


Description


\*---------------------------------------------------------------------------*/



// #include "Kokkos_Core.hpp"

#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/fields/FieldTypeDefs.hpp"
#include "NeoFOAM/mesh/unstructuredMesh/unstructuredMesh.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/grad/gaussGreenGrad.hpp"

#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"

#include "NeoFOAM_GPL/readers/foamMesh.hpp"
#include "NeoFOAM_GPL/writers/writers.hpp"
#include "NeoFOAM_GPL/setup/setup.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{   
    Kokkos::initialize(argc, argv);
    {
        #include "addProfilingOption.H"
        #include "addCheckCaseOptions.H"
        #include "setRootCase.H"
        #include "createTime.H"
        #include "createMesh.H"
        #include "createControl.H"
        // #include "createTimeControls.H"
        auto [adjustTimeStep, maxCo, maxDeltaT] = timeControls(runTime);


        #include "createFields.H"
        // set temperature
        Foam::scalar spread = 0.05;
        forAll(T,celli)
        {
            // T[celli] = Foam::exp(10*(1 - (Foam::mag(mesh.C()[celli] - Foam::vector(0.5,0.75,0.0)))))/Foam::constant::mathematical::e;
            T[celli] = std::exp(-0.5 * (std::pow((mesh.C()[celli].x() - 0.5) / spread, 2.0) + std::pow((mesh.C()[celli].y() - 0.75) / spread, 2.0)));
        }
        T.correctBoundaryConditions();
        T.write();
        // #include "readTimeControls.H"
        // [adjustTimeStep, maxCo, maxDeltaT] = createTimeControls(runTime);
        // updateTimeControls(runTime, adjustTimeStep, maxCo, maxDeltaT);
        std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);
        // #include "createUfIfPresent.H"
        // #include "CourantNo.H"
        Foam::scalar CoNum = calculateCoNum(phi);
        if (adjustTimeStep)
        {
            setDeltaT(runTime, maxCo,CoNum, maxDeltaT);
        }

        Foam::scalar pi = Foam::constant::mathematical::pi;
        {
            Foam::scalarField X = mesh.C().component(0);
            Foam::scalarField Y = mesh.C().component(1);
            Foam::scalarField u = -Foam::sin(2.0*pi*Y)*Foam::pow(Foam::sin(pi*X),2.0);
            Foam::scalarField w = Foam::sin(2.0*pi*X)*Foam::pow(Foam::sin(pi*Y),2.0);
            forAll(U0,celli)
            {
                U0[celli].x() = u[celli];
                U0[celli].y() = w[celli];
                U0[celli].z() = 0.0;
            }
        }
        phi0 = linearInterpolate(U0) & mesh.Sf();

        while (runTime.run())
        {
            // #include "readTimeControls.H"
            std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);
            CoNum = calculateCoNum(phi);
            Foam::Info << "max(phi) : " << max(phi) << Foam::endl;
            Foam::Info << "max(U) : " << max(U) << Foam::endl;
            if (adjustTimeStep)
            {
                setDeltaT(runTime, maxCo,CoNum, maxDeltaT);
            }

            runTime++;

            if(spirallingFlow > 0)
            {
                Foam::Info << "Spiralling flow: " << spirallingFlow << Foam::endl;
                Foam::scalar t = runTime.time().value();
                Foam::scalar dt = runTime.deltaT().value();
                U = U0*Foam::cos(pi*(t+ 0.5*dt)/spirallingFlow);
                phi = phi0*Foam::cos(pi*(t+ 0.5*dt)/spirallingFlow);
                // Foam::volScalarField Test("Test", Foam::fvc::div(phi, T));
                // Test.write();
            }

            Foam::fvScalarMatrix TEqn
            (
                Foam::fvm::ddt(T)
                + Foam::fvc::div(phi, T)
            );

            TEqn.solve();

            runTime.write();
        }

        Foam::Info << "End\n" << Foam::endl;
    }
    Kokkos::finalize();

    return 0;
}

// ************************************************************************* //

