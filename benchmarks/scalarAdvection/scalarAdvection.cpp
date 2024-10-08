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

#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/readers/foamFields.hpp"

#include "FoamAdapter/writers/writers.hpp"
#include "FoamAdapter/setup/setup.hpp"

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
        // creating neofoam fields
        Foam::Info << "creating neofoam mesh" << Foam::endl;
        NeoFOAM::UnstructuredMesh uMesh = Foam::readOpenFOAMMesh(exec, mesh);
        fvcc::VolumeField<NeoFOAM::scalar> neoT = Foam::constructFrom(exec, uMesh, T);
        neoT.correctBoundaryConditions();
        // fvcc::VolumeField<NeoFOAM::Vector> neoU = constructFrom(exec, uMesh, U);

        //     std::pow((s_cc[celli][1] - 0.75) / spread, 2.0)));
        // });
        // neoT.correctBoundaryConditions();
        fvcc::SurfaceField<NeoFOAM::scalar> neoPhi = constructSurfaceField(exec, uMesh, phi);

        Foam::Info << "writing neoT field" << Foam::endl;
        write(neoT.internalField(), mesh, "neoT");

        // #include "readTimeControls.H"
        // [adjustTimeStep, maxCo, maxDeltaT] = createTimeControls(runTime);
        // updateTimeControls(runTime, adjustTimeStep, maxCo, maxDeltaT);
        std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);
        // #include "createUfIfPresent.H"
        // #include "CourantNo.H"
        Foam::scalar coNum = Foam::calculateCoNum(phi);
        if (adjustTimeStep)
        {
            Foam::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
        }

        Foam::scalar pi = Foam::constant::mathematical::pi;
        {
            Foam::scalarField x(mesh.C().component(0));
            Foam::scalarField y(mesh.C().component(1));
            Foam::scalarField u(-Foam::sin(2.0 * pi * y) * Foam::pow(Foam::sin(pi * x), 2.0));
            Foam::scalarField w(Foam::sin(2.0 * pi * x) * Foam::pow(Foam::sin(pi * y), 2.0));
            forAll(u, celli)
            {
                U0[celli].x() = u[celli];
                U0[celli].y() = w[celli];
                U0[celli].z() = 0.0;
            }
        }
        phi0 = Foam::linearInterpolate(U0) & mesh.Sf();
        fvcc::SurfaceField<NeoFOAM::scalar> neoPhi0 = constructSurfaceField(exec, uMesh, phi0);

        Foam::volScalarField ofDivT("ofDivT", Foam::fvc::div(phi, T));

        fvcc::VolumeField<NeoFOAM::scalar> neoDivT = constructFrom(exec, uMesh, ofDivT);
        NeoFOAM::fill(neoDivT.internalField(), 0.0);
        NeoFOAM::fill(neoDivT.boundaryField().value(), 0.0);


        while (runTime.run())
        {
            // #include "readTimeControls.H"
            std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);
            coNum = calculateCoNum(phi);
            Foam::Info << "max(phi) : " << max(phi) << Foam::endl;
            Foam::Info << "max(U) : " << max(U) << Foam::endl;
            if (adjustTimeStep)
            {
                Foam::setDeltaT(runTime, maxCo, coNum, maxDeltaT);
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
                neoPhi.internalField() =
                    neoPhi0.internalField() * std::cos(pi * (t + 0.5 * dt) / spirallingFlow);
            }

            {
                addProfiling(foamAdvection, "foamAdvection");
                Foam::fvScalarMatrix tEqn(Foam::fvm::ddt(T) + Foam::fvc::div(phi, T));

                tEqn.solve();
            }

            // NeoFOAM Euler hardcoded
            {
                addProfiling(neoFoamAdvection, "neoFoamAdvection");
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
