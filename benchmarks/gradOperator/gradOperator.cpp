// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

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
#include "NeoFOAM/fields/field.hpp"
#include "NeoFOAM/mesh/unstructured.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenGrad.hpp"

#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/writers/writers.hpp"
#include "FoamAdapter/setup/setup.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<typename T>
void printField(NeoFOAM::Field<T> a)
{
    std::cout << "a has a size of: " << a.size() << std::endl;
    auto tmpViewHost = a.copyToHost();
    auto tmpView = tmpViewHost.span();
    for (int i = 0; i < a.size(); i++)
    {
        std::cout << "tmp_view: " << tmpView[i] << " at: " << i << std::endl;
    }
}

template<>
void printField(NeoFOAM::Field<NeoFOAM::Vector> a)
{
    std::cout << "a has a size of: " << a.size() << std::endl;
    auto tmpViewHost = a.copyToHost();
    auto tmpView = tmpViewHost.span();
    for (int i = 0; i < a.size(); i++)
    {
        std::cout << "tmp_view: " << tmpView[i](0) << " " << tmpView[i](1) << " " << tmpView[i](2)
                  << " at: " << i << std::endl;
    }
}

int main(int argc, char* argv[])
{
    Kokkos::initialize(argc, argv);
    {
#include "addProfilingOption.H"
#include "addCheckCaseOptions.H"
#include "setRootCase.H"
#include "createTime.H"
        std::unique_ptr<Foam::fvMesh> meshPtr = Foam::createMesh(runTime);
        Foam::fvMesh& mesh = *meshPtr;
#include "createFields.H"

        runTime++;

        for (int celli = 0; celli < mesh.nCells(); celli++)
        {
            T[celli] = celli;
        }
        T.write();
        {
            addProfiling(foamGradOperator, "foamGradOperator");
            Foam::volVectorField test(Foam::fvc::grad(T));
            test.write();
        }

        {
            NeoFOAM::Executor exec = NeoFOAM::GPUExecutor();
            Foam::Info << "reading mesh" << Foam::endl;
            NeoFOAM::UnstructuredMesh uMesh = readOpenFOAMMesh(exec, mesh);

            NeoFOAM::scalarField temperature(exec, T.internalField().size());
            temperature.apply(KOKKOS_LAMBDA(const NeoFOAM::size_t i) { return NeoFOAM::scalar(i); }
            );

            Foam::Info << "writing temperature field" << Foam::endl;
            write(temperature, mesh, "temperature");

            {
                addProfiling(neofoamGradOperator, "neofoamGradOperator");
                // NeoFOAM::vectorField gradT = NeoFOAM::gaussGreenGrad(exec,
                // uMesh).grad(temperature); Kokkos::fence(); write(gradT, mesh, "gradTGPU");
            }
        }

        {
            NeoFOAM::Executor exec = NeoFOAM::SerialExecutor();
            Foam::Info << "reading mesh" << Foam::endl;
            NeoFOAM::UnstructuredMesh uMesh = readOpenFOAMMesh(exec, mesh);

            NeoFOAM::scalarField temperature(exec, T.internalField().size());
            temperature.apply(KOKKOS_LAMBDA(const NeoFOAM::size_t i) { return NeoFOAM::scalar(i); }
            );

            Foam::Info << "writing temperature field" << Foam::endl;
            write(temperature, mesh, "temperature");

            {
                addProfiling(neofoamGradOperator, "neofoamGradOperator");
                // NeoFOAM::vectorField gradT = NeoFOAM::gaussGreenGrad(exec,
                // uMesh).grad(temperature); Kokkos::fence(); write(gradT, mesh, "gradTCPU");
            }
        }

        {
            NeoFOAM::Executor exec = NeoFOAM::CPUExecutor();
            Foam::Info << "reading mesh" << Foam::endl;
            NeoFOAM::UnstructuredMesh uMesh = readOpenFOAMMesh(exec, mesh);

            NeoFOAM::scalarField temperature(exec, T.internalField().size());
            temperature.apply(KOKKOS_LAMBDA(const NeoFOAM::size_t i) { return NeoFOAM::scalar(i); }
            );

            Foam::Info << "writing temperature field" << Foam::endl;
            write(temperature, mesh, "temperature");

            {
                addProfiling(neofoamGradOperator, "neofoamGradOperator");
                //     NeoFOAM::vectorField gradT = NeoFOAM::gaussGreenGrad(exec,
                //     uMesh).grad(temperature); Kokkos::fence(); write(gradT, mesh, "gradTOMP");
            }
        }

        Foam::profiling::print(Foam::Info);
        runTime.write();

        Foam::Info << "End\n" << Foam::endl;
    }

    Kokkos::finalize();

    return 0;
}

// ************************************************************************* //
