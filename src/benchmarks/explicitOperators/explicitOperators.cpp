// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors


#include "fvCFD.H"
#include "profiling.H"
// #include "fields.H"

// void (const fvMesh& mesh,asdf)
// {

// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
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

    Info << "Time = " << runTime.timeName() << nl << endl;


    volVectorField gradT(fvc::grad(T));
    for (int i; i < 10; i++)
    {
        addProfiling(grad, "OpenFOAM");
        gradT = fvc::grad(T);
    }

    {
        addProfiling(grad, "NeoFOAM");
    }


    profiling::writeNow();

    runTime.printExecutionTime(Info);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
