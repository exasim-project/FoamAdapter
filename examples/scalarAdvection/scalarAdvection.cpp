// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/fields/field.hpp"
#include "NeoFOAM/mesh/unstructured.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenDiv.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenGrad.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/interpolation/surfaceInterpolation.hpp"

#include "fvCFD.H"

#include "FoamAdapter/meshAdapter.hpp"
#include "FoamAdapter/readers/foamFields.hpp"

#include "FoamAdapter/writers.hpp"
#include "FoamAdapter/setup.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace fvcc = NeoFOAM::finiteVolume::cellCentred;

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

        Info << "creating NeoFOAM mesh" << endl;
        auto nfMesh = mesh.nfMesh();

        Info << "creating NeoFOAM fields" << endl;
        auto nfT = constructFrom(exec, nfMesh, T);
        auto nfPhi = constructSurfaceField(exec, nfMesh, phi);

        volScalarField ofDivT("ofDivT", fvc::div(phi, T));
        auto nfDivT = constructFrom(exec, nfMesh, ofDivT);

        while (runTime.run())
        {
            runTime++;

            Info << "Time = " << runTime.timeName() << nl << max(phi) << nl << max(U) << endl;

            // NeoFOAM Euler
            // NOTE for now hardcoded
            // this will soon be replaced by the NeoFOAM DSL
            {
                fvcc::GaussGreenDiv(
                    exec,
                    nfMesh,
                    fvcc::SurfaceInterpolation(
                        exec,
                        nfMesh,
                        fvcc::SurfaceInterpolationFactory::create("upwind", exec, nfMesh)
                    )
                )
                    .div(nfDivT, nfPhi, nfT);

                nfT.internalField() =
                    nfT.internalField() - nfDivT.internalField() * runTime.deltaT().value();
                nfT.correctBoundaryConditions();
                Kokkos::fence();
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
