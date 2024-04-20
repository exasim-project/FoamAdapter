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
#include "NeoFOAM/cellCentredFiniteVolume/div/gaussGreenDiv.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/surfaceInterpolation/surfaceInterpolationSelector.hpp"

#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"

#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/readers/foamFields.hpp"

#include "FoamAdapter/writers/writers.hpp"
#include "FoamAdapter/setup/setup.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{   
    Kokkos::initialize(argc, argv);
    
    {
        #include "addProfilingOption.H"
        #include "addCheckCaseOptions.H"
        #include "setRootCase.H"
        #include "createTime.H"
        NeoFOAM::executor exec = Foam::createExecutor(runTime.controlDict());
        
        std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec,runTime);
        Foam::fvccNeoMesh& mesh = *meshPtr;
        #include "createControl.H"
        // #include "createTimeControls.H"
        auto [adjustTimeStep, maxCo, maxDeltaT] = Foam::timeControls(runTime);


        #include "createFields.H"
        // set temperature
        Foam::scalar spread = 0.05;
        forAll(T,celli)
        {
            T[celli] = std::exp(-0.5 * (std::pow((mesh.C()[celli].x() - 0.5) / spread, 2.0) + std::pow((mesh.C()[celli].y() - 0.75) / spread, 2.0)));
        }
        T.correctBoundaryConditions();
        T.write();
        // creating neofoam fields
        Foam::Info << "creating neofoam mesh" << Foam::endl;
        NeoFOAM::unstructuredMesh uMesh = Foam::readOpenFOAMMesh(exec, mesh);
        NeoFOAM::fvccVolField<NeoFOAM::scalar> neoT = Foam::constructFrom(exec, uMesh, T);
        neoT.correctBoundaryConditions();
        // NeoFOAM::fvccVolField<NeoFOAM::Vector> neoU = constructFrom(exec, uMesh, U);

        // auto s_cc = uMesh.cellCentres().field();
        // neoT.internalField().apply(KOKKOS_LAMBDA(int celli) 
        // {
        //     return std::exp(-0.5 * (std::pow((s_cc[celli][0] - 0.5) / spread, 2.0) + std::pow((s_cc[celli][1] - 0.75) / spread, 2.0)));
        // });
        // neoT.correctBoundaryConditions();
        NeoFOAM::fvccSurfaceField<NeoFOAM::scalar> neoPhi = constructSurfaceField(exec, uMesh, phi);
        
        Foam::Info << "writing neoT field" << Foam::endl;
        write(neoT.internalField(), mesh, "neoT");

        // #include "readTimeControls.H"
        // [adjustTimeStep, maxCo, maxDeltaT] = createTimeControls(runTime);
        // updateTimeControls(runTime, adjustTimeStep, maxCo, maxDeltaT);
        std::tie(adjustTimeStep, maxCo, maxDeltaT) = timeControls(runTime);
        // #include "createUfIfPresent.H"
        // #include "CourantNo.H"
        Foam::scalar CoNum = Foam::calculateCoNum(phi);
        if (adjustTimeStep)
        {
            Foam::setDeltaT(runTime, maxCo,CoNum, maxDeltaT);
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
        phi0 = Foam::linearInterpolate(U0) & mesh.Sf();
        NeoFOAM::fvccSurfaceField<NeoFOAM::scalar> neoPhi0 = constructSurfaceField(exec, uMesh, phi0);

        Foam::volScalarField ofDivT("ofDivT",Foam::fvc::div(phi, T));

        NeoFOAM::fvccVolField<NeoFOAM::scalar> neoDivT = constructFrom(exec, uMesh, ofDivT);
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
                Foam::setDeltaT(runTime, maxCo,CoNum, maxDeltaT);
            }

            runTime++;

            Foam::Info << "Time = " << runTime.timeName() << Foam::nl << Foam::endl;

            if(spirallingFlow > 0)
            {
                Foam::Info << "Spiralling flow: " << spirallingFlow << Foam::endl;
                Foam::scalar t = runTime.time().value();
                Foam::scalar dt = runTime.deltaT().value();
                U = U0*Foam::cos(pi*(t+ 0.5*dt)/spirallingFlow);
                phi = phi0*Foam::cos(pi*(t+ 0.5*dt)/spirallingFlow);
                neoPhi.internalField() = neoPhi0.internalField()*std::cos(pi*(t+ 0.5*dt)/spirallingFlow);
            }

            {
                addProfiling(foamAdvection, "foamAdvection");
                Foam::fvScalarMatrix TEqn
                (
                    Foam::fvm::ddt(T)
                    + Foam::fvc::div(phi, T)
                );

                TEqn.solve();
            }

            // NeoFOAM Euler hardcoded
            {
                addProfiling(neoFoamAdvection, "neoFoamAdvection");
                NeoFOAM::fill(neoDivT.internalField(), 0.0);
                NeoFOAM::fill(neoDivT.boundaryField().value(), 0.0);
                NeoFOAM::gaussGreenDiv(
                    exec,
                    uMesh,
                    NeoFOAM::surfaceInterpolationSelector(std::string("upwind"),exec ,mesh.uMesh())
                ).div(neoDivT, neoPhi ,neoT);
                neoT.internalField() = neoT.internalField() - neoDivT.internalField() * runTime.deltaT().value();
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

