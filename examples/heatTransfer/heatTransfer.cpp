// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 NeoFOAM authors

#include "NeoN/NeoN.hpp"

#include "FoamAdapter/FoamAdapter.hpp"
#include "FoamAdapter/readers/foamDictionary.hpp"


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


        NeoN::Dictionary controlDict = Foam::readFoamDictionary(runTime.controlDict());
        NeoN::Executor exec = createExecutor(runTime.controlDict());

        std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, runTime);
        Foam::MeshAdapter& mesh = *meshPtr;

        Foam::pisoControl piso(mesh);

        auto [adjustTimeStep, maxCo, maxDeltaT] = Foam::timeControls(runTime);

#include "createFields.H"


        NeoN::Dictionary fvSchemesDict = Foam::readFoamDictionary(mesh.schemesDict());
        NeoN::Dictionary fvSolutionDict = Foam::readFoamDictionary(mesh.solutionDict());

        Info << "creating NeoFOAM mesh" << endl;
        NeoN::UnstructuredMesh& nfMesh = mesh.nfMesh();

        Info << "creating NeoFOAM fields" << endl;
        fvcc::VolumeField<NeoN::scalar>& nfT =
            vectorCollection.registerVector<fvcc::VolumeField<NeoN::scalar>>(
                Foam::CreateFromFoamField<Foam::volScalarField> {
                    .exec = exec,
                    .nfMesh = nfMesh,
                    .foamField = T,
                    .name = "nfT"
                }
            );
        auto& nfTOld = oldTime(nfT);
        nfTOld.internalVector() = nfT.internalVector();
        nfT.correctBoundaryConditions();

        auto nfKappa = Foam::constructSurfaceField(exec, nfMesh, kappa);


        Foam::scalar endTime = controlDict.get<Foam::scalar>("endTime");


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        Info << "\nStarting time loop\n" << endl;

        while (runTime.loop())
        {

            runTime++;

            Foam::scalar t = runTime.time().value();
            Foam::scalar dt = runTime.deltaT().value();

            Info << "Time = " << runTime.timeName() << nl << endl;

            Foam::fvScalarMatrix TEqn(fvm::ddt(T) - fvm::laplacian(kappa, T));

            TEqn.solve();

            dsl::Expression nfTEqn(dsl::imp::ddt(nfT) - dsl::imp::laplacian(nfKappa, nfT));

            dsl::solve(
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
                Info << "writing nfT field" << endl;
                write(nfT.internalVector(), mesh, "nfT");
            }

            runTime.printExecutionTime(Info);
        }

        Info << "End\n" << endl;
    }
    Kokkos::finalize();

    return 0;
}

// ************************************************************************* //
