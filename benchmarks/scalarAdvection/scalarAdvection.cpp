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

#include "NeoFOAM/DSL/eqnTerm.hpp"
#include "NeoFOAM/DSL/eqnSystem.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/timeIntegration/timeIntegration.hpp"


namespace dsl = NeoFOAM::DSL;
namespace fvcc = NeoFOAM::finiteVolume::cellCentred;

class Temporal
{

public:

    Temporal(fvcc::VolumeField<NeoFOAM::scalar>& Phi)
        : termType_(dsl::EqnTerm::Type::Temporal), Phi_(Phi), exec_(Phi.exec()),
          nCells_(Phi.mesh().nCells())
    {}

    std::string display() const { return "Temporal"; }

    void temporalOperation(NeoFOAM::Field<NeoFOAM::scalar>& field, NeoFOAM::scalar scale) {}

    dsl::EqnTerm::Type getType() const { return termType_; }

    fvcc::VolumeField<NeoFOAM::scalar>* volumeField() { return &Phi_; }

    const NeoFOAM::Executor& exec() const { return exec_; }

    std::size_t nCells() const { return nCells_; }

    dsl::EqnTerm::Type termType_;


    fvcc::VolumeField<NeoFOAM::scalar>& Phi_;
    const NeoFOAM::Executor exec_;
    const std::size_t nCells_;
};

class Divergence
{

public:

    Divergence(
        const fvcc::SurfaceField<NeoFOAM::scalar>& faceFlux, fvcc::VolumeField<NeoFOAM::scalar>& Phi
    )
        : termType_(dsl::EqnTerm::Type::Explicit), exec_(Phi.exec()), nCells_(Phi.mesh().nCells()),
          faceFlux_(faceFlux), Phi_(Phi),
          div_(
              Phi.exec(),
              Phi.mesh(),
              fvcc::SurfaceInterpolation(
                  Phi.exec(),
                  Phi.mesh(),
                  fvcc::SurfaceInterpolationFactory::create("upwind", Phi.exec(), Phi.mesh())
              )
          )
    {}

    std::string display() const { return "Divergence"; }

    void explicitOperation(NeoFOAM::Field<NeoFOAM::scalar>& source, NeoFOAM::scalar scale)
    {
        div_.div(source, faceFlux_, Phi_);
    }

    dsl::EqnTerm::Type getType() const { return termType_; }

    fvcc::VolumeField<NeoFOAM::scalar>* volumeField() { return nullptr; }


    const NeoFOAM::Executor& exec() const { return exec_; }

    std::size_t nCells() const { return nCells_; }

    dsl::EqnTerm::Type termType_;

    const NeoFOAM::Executor exec_;
    const std::size_t nCells_;
    const fvcc::SurfaceField<NeoFOAM::scalar>& faceFlux_;
    fvcc::VolumeField<NeoFOAM::scalar>& Phi_;
    fvcc::GaussGreenDiv div_;
};


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
        fvcc::SurfaceField<NeoFOAM::scalar> neoPhi0 = constructSurfaceField(exec, uMesh, phi0);

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

            {
                addProfiling(foamAdvection, "foamAdvection");
                Foam::fvScalarMatrix TEqn(Foam::fvm::ddt(T) + Foam::fvc::div(phi, T));

                TEqn.solve();
            }

            // NeoFOAM Euler hardcoded
            {
                addProfiling(neoFoamAdvection, "neoFoamAdvection");
                
                dsl::EqnTerm neoTimeTerm = Temporal(neoT);
                dsl::EqnTerm neoDivTerm = Divergence(neoPhi, neoT);
                dsl::EqnSystem eqnSys = neoTimeTerm + neoDivTerm;
                eqnSys.dt = runTime.deltaT().value();

                NeoFOAM::Dictionary dict;
                dict.insert("type", std::string("forwardEuler"));
                
                fvcc::TimeIntegration timeIntergrator(eqnSys, dict);
                timeIntergrator.solve();
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
