// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/fields/field.hpp"
#include "NeoFOAM/mesh/unstructured.hpp"
#include "NeoFOAM/finiteVolume/operators/gaussGreenDiv.hpp"
#include "NeoFOAM/finiteVolume/operators/gaussGreenGrad.hpp"
#include "NeoFOAM/finiteVolume/interpolation/surfaceInterpolation.hpp"

#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"

#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/readers/foamFields.hpp"

#include "FoamAdapter/writers/writers.hpp"
#include "FoamAdapter/setup/setup.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

        NeoFOAM::fvccSurfaceField<NeoFOAM::scalar> neoPhi = constructSurfaceField(exec, uMesh, phi);

        Foam::Info << "writing neoT field" << Foam::endl;
        write(neoT.internalField(), mesh, "neoT");

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
        NeoFOAM::fvccSurfaceField<NeoFOAM::scalar> neoPhi0 =
            constructSurfaceField(exec, uMesh, phi0);

        Foam::volScalarField ofDivT("ofDivT", Foam::fvc::div(phi, T));

        fvcc::VolumeField<NeoFOAM::scalar> neoDivT = constructFrom(exec, uMesh, ofDivT);
        NeoFOAM::fill(neoDivT.internalField(), 0.0);
        NeoFOAM::fill(neoDivT.boundaryField().value(), 0.0);


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
                neoPhi.internalField() =
                    neoPhi0.internalField() * std::cos(pi * (t + 0.5 * dt) / spirallingFlow);
            }

            // NeoFOAM Euler hardcoded
            {
                addProfiling(neoFoamAdvection, "neoFoamAdvection");
                NeoFOAM::fill(neoDivT.internalField(), 0.0);
                NeoFOAM::fill(neoDivT.boundaryField().value(), 0.0);
                NeoFOAM::gaussGreenDiv(
                    exec,
                    uMesh,
                    NeoFOAM::surfaceInterpolationSelector(std::string("upwind"), exec, mesh.uMesh())
                )
                    .div(neoDivT, neoPhi, neoT);
                // neoT.internalField() = neoT.internalField() - neoDivT.internalField() *
                // runTime.deltaT().value();
                double dt = runTime.deltaT().value();
                auto s_neoT = neoT.internalField().span();
                auto s_neoDivT = neoDivT.internalField().span();
                neoT.internalField().apply(KOKKOS_LAMBDA(const int celli) {
                    return s_neoT[celli] + dt * s_neoDivT[celli];
                });
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
