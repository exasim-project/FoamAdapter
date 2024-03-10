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
#include "NeoFOAM/fields/field.hpp"
#include "NeoFOAM/mesh/unstructuredMesh/unstructuredMesh.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/grad/gaussGreenGrad.hpp"

#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"

#include "NeoFOAM_GPL/readers/foamMesh.hpp"
#include "NeoFOAM_GPL/writers/writers.hpp"



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
        #include "createTimeControls.H"

        #include "createFields.H"
        #include "readTimeControls.H"
        #include "createUfIfPresent.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"


        scalar pi = constant::mathematical::pi;
        {
            scalarField X = mesh.C().component(0);
            scalarField Z = mesh.C().component(2);
            scalarField u = -Foam::sin(2.0*pi*Z)*Foam::pow(Foam::sin(pi*X),2.0);
            scalarField w = Foam::sin(2.0*pi*X)*Foam::pow(Foam::sin(pi*Z),2.0);
            forAll(U0,celli)
            {
                U0[celli].x() = u[celli];
                U0[celli].y() = 0.0;
                U0[celli].z() = w[celli];
            }
        }

        phi = linearInterpolate(U) & mesh.Sf();

        while (runTime.run())
        {
            #include "readTimeControls.H"
            #include "CourantNo.H"
            #include "setDeltaT.H"

            runTime++;

            if(spirallingFlow > 0)
            {
                scalar t = runTime.time().value();
                scalar dt = runTime.deltaT().value();
                U = U0*Foam::cos(pi*(t+ 0.5*dt)/spirallingFlow);
                phi = phi0*Foam::cos(pi*(t+ 0.5*dt)/spirallingFlow);
            }

            fvScalarMatrix UEqn
            (
                fvm::ddt(T)
                + fvc::div(phi, T)
            );


            
            runTime.write();
        }

        Info << "End\n" << Foam::endl;
    }
    Kokkos::finalize();

    return 0;
}

// ************************************************************************* //

