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

#include "Time.H"
#include "fvMesh.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "fvm.H"
#include "linear.H"
#include "uniformDimensionedFields.H"
#include "calculatedFvPatchFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "findRefCell.H"
#include "IOMRFZoneList.H"
#include "constants.H"
#include "gravityMeshObject.H"

#include "columnFvMesh.H"

#include "OSspecific.H"
#include "argList.H"
#include "timeSelector.H"

#include "profiling.H"

#include "Kokkos_Core.hpp"

#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/fields/FieldTypeDefs.hpp"
#include "NeoFOAM/mesh/unstructuredMesh/unstructuredMesh.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/grad/gaussGreenGrad.hpp"

#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/writers/writers.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <typename T>
void print_field(NeoFOAM::Field<T> a)
{
    std::cout << "a has a size of: " << a.size() << std::endl;
    auto tmp_view = a.copyToHost().field();
    for (int i = 0; i < a.size(); i++)
    {
        std::cout << "tmp_view: " << tmp_view[i] << " at: " << i << std::endl;
    }
}
template <>
void print_field(NeoFOAM::Field<NeoFOAM::Vector> a)
{
    std::cout << "a has a size of: " << a.size() << std::endl;
    auto tmp_view = a.copyToHost().field();
    for (int i = 0; i < a.size(); i++)
    {
        std::cout << "tmp_view: " << tmp_view[i](0) << " " << tmp_view[i](1) << " " << tmp_view[i](2) << " at: " << i << std::endl;
    }
}

int main(int argc, char *argv[])
{
    Kokkos::initialize(argc, argv);
    {
#include "addProfilingOption.H"
#include "addCheckCaseOptions.H"
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"

        runTime++;
        int N = 1000;

        for (int celli = 0; celli < mesh.nCells(); celli++)
        {
            T[celli] = celli;
        }
        T.write();
        // for (int i = 0; i < N; i++)
        {
            addProfiling(foamGradOperator, "foamGradOperator");
            Foam::volVectorField test = Foam::fvc::grad(T);
            test.write();
        }


        {
            NeoFOAM::executor exec = NeoFOAM::GPUExecutor();
            Foam::Info << "reading mesh" << Foam::endl;
            NeoFOAM::unstructuredMesh uMesh = readOpenFOAMMesh(exec, mesh);

            NeoFOAM::scalarField Temperature(exec, T.internalField().size());
            Temperature.apply(KOKKOS_LAMBDA(int i) { return i; });

            Foam::Info << "writing temperature field" << Foam::endl;
            write(Temperature, mesh, "Temperature");

            {
                addProfiling(neofoamGradOperator, "neofoamGradOperator");
                // NeoFOAM::vectorField gradT = NeoFOAM::gaussGreenGrad(exec, uMesh).grad(Temperature);
                // Kokkos::fence();
                // write(gradT, mesh, "gradTGPU");
            }
        }

        {
            NeoFOAM::executor exec = NeoFOAM::CPUExecutor();
            Foam::Info << "reading mesh" << Foam::endl;
            NeoFOAM::unstructuredMesh uMesh = readOpenFOAMMesh(exec, mesh);

            NeoFOAM::scalarField Temperature(exec, T.internalField().size());
            Temperature.apply(KOKKOS_LAMBDA(int i) { return i; });

            Foam::Info << "writing temperature field" << Foam::endl;
            write(Temperature, mesh, "Temperature");

            {
                addProfiling(neofoamGradOperator, "neofoamGradOperator");
                // NeoFOAM::vectorField gradT = NeoFOAM::gaussGreenGrad(exec, uMesh).grad(Temperature);
                // Kokkos::fence();
                // write(gradT, mesh, "gradTCPU");
            }
        }

        {
            NeoFOAM::executor exec = NeoFOAM::OMPExecutor();
            Foam::Info << "reading mesh" << Foam::endl;
            NeoFOAM::unstructuredMesh uMesh = readOpenFOAMMesh(exec, mesh);

            NeoFOAM::scalarField Temperature(exec, T.internalField().size());
            Temperature.apply(KOKKOS_LAMBDA(int i) { return i; });

            Foam::Info << "writing temperature field" << Foam::endl;
            write(Temperature, mesh, "Temperature");

            {
                addProfiling(neofoamGradOperator, "neofoamGradOperator");
            //     NeoFOAM::vectorField gradT = NeoFOAM::gaussGreenGrad(exec, uMesh).grad(Temperature);
            //     Kokkos::fence();
            //     write(gradT, mesh, "gradTOMP");
            }
        }

        Foam::profiling::print(Foam::Info);
        runTime.write();

        Foam::Info << "End\n"
                   << Foam::endl;
    }

    Kokkos::finalize();

    return 0;
}

// ************************************************************************* //
