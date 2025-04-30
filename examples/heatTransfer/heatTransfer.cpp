// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 FoamAdapter authors

#include "NeoN/NeoN.hpp"

#include "FoamAdapter/FoamAdapter.hpp"


#include "fvCFD.H"
#include "pisoControl.H"

using Foam::Info;
using Foam::endl;
using Foam::nl;
namespace fvc = Foam::fvc;
namespace fvm = Foam::fvm;

namespace dsl = NeoN::dsl;
namespace fvcc = NeoN::finiteVolume::cellCentred;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    Kokkos::initialize(argc, argv);
    {
#include "addCheckCaseOptions.H"
#include "setRootCase.H"
#include "createTime.H"

        NeoN::Database db;

        fvcc::VectorCollection& vectorCollection =
            fvcc::VectorCollection::instance(db, "VectorCollection");


        NeoN::Dictionary controlDict = FoamAdapter::convert(runTime.controlDict());
        NeoN::Executor exec = FoamAdapter::createExecutor(runTime.controlDict());

        std::unique_ptr<FoamAdapter::MeshAdapter> meshPtr = FoamAdapter::createMesh(exec, runTime);
        FoamAdapter::MeshAdapter& mesh = *meshPtr;

        Foam::pisoControl piso(mesh);

        auto [adjustTimeStep, maxCo, maxDeltaT] = FoamAdapter::timeControls(runTime);

#include "createFields.H"

        NeoN::Dictionary fvSchemesDict = FoamAdapter::convert(mesh.schemesDict());
        NeoN::Dictionary fvSolutionDict = FoamAdapter::convert(mesh.solutionDict());

        Info << "creating FoamAdapter mesh" << endl;
        NeoN::UnstructuredMesh& nfMesh = mesh.nfMesh();

        Info << "creating FoamAdapter fields" << endl;
        fvcc::VolumeField<NeoN::scalar>& nfT =
            vectorCollection.registerVector<fvcc::VolumeField<NeoN::scalar>>(
                FoamAdapter::CreateFromFoamField<Foam::volScalarField> {
                    .exec = exec,
                    .nfMesh = nfMesh,
                    .foamField = T,
                    .name = "nfT"
                }
            );
        auto& nfTOld = oldTime(nfT);
        nfTOld.internalVector() = nfT.internalVector();
        nfT.correctBoundaryConditions();

        auto nfKappa = FoamAdapter::constructSurfaceField(exec, nfMesh, kappa);

        Foam::scalar endTime = controlDict.get<Foam::scalar>("endTime");

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        Foam::Info << "\nStarting time loop\n" << Foam::endl;

        while (runTime.loop())
        {

            runTime++;

            Foam::scalar t = runTime.time().value();
            Foam::scalar dt = runTime.deltaT().value();

            Foam::Info << "Time = " << runTime.timeName() << Foam::nl << Foam::endl;

            Foam::fvScalarMatrix TEqn(Foam::fvm::ddt(T) - Foam::fvm::laplacian(kappa, T));

            TEqn.solve();

            NeoN::dsl::Expression nfTEqn(
                NeoN::dsl::imp::ddt(nfT) - NeoN::dsl::imp::laplacian(nfKappa, nfT)
            );

            NeoN::dsl::solve(
                nfTEqn,
                nfT,
                t,
                dt,
                fvSchemesDict,
                fvSolutionDict.get<NeoN::Dictionary>("solvers").get<NeoN::Dictionary>("nfT")
            );

            runTime.write();
            if (runTime.outputTime())
            {
                Foam::Info << "writing nfT field" << Foam::endl;
                write(nfT.internalVector(), mesh, "nfT");
            }

            runTime.printExecutionTime(Info);
        }

        Foam::Info << "End\n" << Foam::endl;
    }
    Kokkos::finalize();

    return 0;
}

// ************************************************************************* //
