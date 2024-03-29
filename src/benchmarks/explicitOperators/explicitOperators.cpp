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
    matrixAssembly


Description


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "profiling.H"
// #include "fields.H"

// void (const fvMesh& mesh,asdf)
// {

// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "addProfilingOption.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"


    #include "createFields.H"
    runTime.setDeltaT(1.0);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    ++runTime;

    Info<< "Time = " << runTime.timeName() << nl << endl;
    

    volVectorField gradT(fvc::grad(T));
    for (int i; i<10; i++)
    {
        addProfiling(grad, "OpenFOAM");
        gradT = fvc::grad(T);
    }

    {
        addProfiling(grad, "NeoFOAM");
        
    }


    profiling::writeNow();
    
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
