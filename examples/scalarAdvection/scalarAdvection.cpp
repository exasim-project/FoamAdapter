// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "FoamAdapter/FoamAdapter.hpp"
#include "NeoFOAM/dsl/expression.hpp"
#include "NeoFOAM/dsl/solver.hpp"
#include "NeoFOAM/dsl/ddt.hpp"
#include "FoamAdapter/readers/foamDictionary.hpp"

// #include "NeoFOAM/dsl/implicit.hpp"
// #include "NeoFOAM/dsl/explicit.hpp"


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
#include "createNFTime.H"
#include "createControl.H"
#include "createFields.H"
        NeoFOAM::Dictionary fvSchemesDict = Foam::readFoamDictionary(mesh.schemesDict());
        NeoFOAM::Dictionary fvSolutionDict = Foam::readFoamDictionary(mesh.solutionDict());

        Info << "creating NeoFOAM mesh" << endl;
        auto nfMesh = mesh.nfMesh();

        Info << "creating NeoFOAM fields" << endl;
        auto nfT = Foam::constructFrom(exec, nfMesh, T);
        auto nfPhi = Foam::constructSurfaceField(exec, nfMesh, phi);

        Foam::volScalarField ofDivT("ofDivT", fvc::div(phi, T));
        auto nfDivT = constructFrom(exec, nfMesh, ofDivT);

        while (runTime.run())
        {
            runTime++;

            Info << "Time = " << runTime.timeName() << nl << max(phi) << nl << max(U) << endl;



            {
                // dsl::Expression eqnSys(
                //     dsl::Implicit::ddt(nfT)
                //   + dsl::Implicit::ddt(nfT)
                // );

                // NeoFOAM::scalar dt = runTime.deltaT().value();
                // dsl::solve(eqnSys, nfT, dt, fvSchemesDict, fvSolutionDict);
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
