// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/fields/field.hpp"
#include "NeoFOAM/mesh/unstructured.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenDiv.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenGrad.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/interpolation/surfaceInterpolation.hpp"

#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"

#include "FoamAdapter/meshAdapter.hpp"
#include "FoamAdapter/readers/foamFields.hpp"

#include "FoamAdapter/writers/writers.hpp"
#include "FoamAdapter/setup.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace fvcc = NeoFOAM::finiteVolume::cellCentred;

int main(int argc, char* argv[])
{
    Kokkos::initialize(argc, argv);

    {
#include "addProfilingOption.H"
#include "addCheckCaseOptions.H"
#include "setRootCase.H"
#include "createTime.H"
        NeoFOAM::Executor exec = Foam::createExecutor(runTime.controlDict());

        std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec, runTime);
        Foam::fvccNeoMesh& mesh = *meshPtr;
#include "createControl.H"
#include "createFields.H"

        auto [adjustTimeStep, maxCo, maxDeltaT] = Foam::timeControls(runTime);

        Foam::Info << "creating NeoFOAM mesh" << Foam::endl;
        NeoFOAM::UnstructuredMesh uMesh = Foam::readOpenFOAMMesh(exec, mesh);

        Foam::Info << "creating NeoFOAM fields" << Foam::endl;
        auto neoT = Foam::constructFrom(exec, uMesh, T);
        neoT.correctBoundaryConditions();
        auto neoPhi = constructSurfaceField(exec, uMesh, phi);

        Foam::Info << "writing neoT field" << Foam::endl;
        write(neoT.internalField(), mesh, "neoT");

        std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);
        Foam::scalar coNum = Foam::calculateCoNum(phi);
        if (adjustTimeStep)
        {
            Foam::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
        }

        Foam::volScalarField ofDivT("ofDivT", Foam::fvc::div(phi, T));

        auto neoDivT = constructFrom(exec, uMesh, ofDivT);
        NeoFOAM::fill(neoDivT.internalField(), 0.0);
        NeoFOAM::fill(neoDivT.boundaryField().value(), 0.0);

        while (runTime.run())
        {
            std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);
            coNum = calculateCoNum(phi);

            if (adjustTimeStep)
            {
                Foam::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
            }

            runTime++;

            Foam::Info << "Time = " << runTime.timeName() << Foam::nl << max(phi) << Foam::nl
                       << max(U) << Foam::endl;

            // NeoFOAM Euler
            // NOTE for now hardcoded
            // this will soon be replaced by the NeoFOAM DSL
            {
                NeoFOAM::fill(neoDivT.internalField(), 0.0);
                NeoFOAM::fill(neoDivT.boundaryField().value(), 0.0);

                fvcc::GaussGreenDiv(
                    exec,
                    uMesh,
                    fvcc::SurfaceInterpolation(
                        exec,
                        uMesh,
                        fvcc::SurfaceInterpolationFactory::create("upwind", exec, uMesh)
                    )
                )
                    .div(neoDivT, neoPhi, neoT);

                neoT.internalField() =
                    neoT.internalField() - neoDivT.internalField() * runTime.deltaT().value();
                neoT.correctBoundaryConditions();
                Kokkos::fence();
            }

            if (runTime.outputTime())
            {
                Foam::Info << "writing neoT field" << Foam::endl;
                write(neoT.internalField(), mesh, "neoT");
            }

            runTime.write();
            runTime.printExecutionTime(Foam::Info);
        }

        Foam::Info << "End\n" << Foam::endl;
    }
    Kokkos::finalize();

    return 0;
}

// ************************************************************************* //
