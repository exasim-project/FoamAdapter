// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023-2025 FoamAdapter authors

#include "FoamAdapter/auxiliary/setup.hpp"
#include "FoamAdapter/datastructures/meshAdapter.hpp"
#include "FoamAdapter/auxiliary/readers.hpp"

#include "fvc.H"

namespace FoamAdapter
{

std::tuple<bool, Foam::scalar, Foam::scalar> timeControls(const Foam::Time& runTime)
{
    bool adjustTimeStep = runTime.controlDict().getOrDefault("adjustTimeStep", false);

    Foam::scalar maxCo = runTime.controlDict().getOrDefault<Foam::scalar>("maxCo", 1);

    Foam::scalar maxDeltaT =
        runTime.controlDict().getOrDefault<Foam::scalar>("maxDeltaT", Foam::GREAT);

    return std::make_tuple(adjustTimeStep, maxCo, maxDeltaT);
}


void setDeltaT(Foam::Time& runTime, Foam::scalar maxCo, Foam::scalar coNum, Foam::scalar maxDeltaT)
{
    Foam::scalar maxDeltaTFact = maxCo / (coNum + Foam::SMALL);
    Foam::scalar deltaTFact = Foam::min(Foam::min(maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact), 1.2);

    runTime.setDeltaT(Foam::min(deltaTFact * runTime.deltaTValue(), maxDeltaT));

    Foam::Info << "deltaT = " << runTime.deltaTValue() << Foam::endl;
}


std::unique_ptr<MeshAdapter> createMesh(const NeoN::Executor& exec, const Foam::Time& runTime)
{
    Foam::word regionName(Foam::polyMesh::defaultRegion);
    Foam::IOobject io(regionName, runTime.timeName(), runTime, Foam::IOobject::MUST_READ);
    return std::make_unique<MeshAdapter>(exec, io);
}

std::unique_ptr<Foam::fvMesh> createMesh(const Foam::Time& runTime)
{
    std::unique_ptr<Foam::fvMesh> meshPtr;
    Foam::Info << "Create mesh";
    Foam::word regionName(Foam::polyMesh::defaultRegion);
    Foam::Info << " for time = " << runTime.timeName() << Foam::nl;

    meshPtr.reset(new Foam::fvMesh(
        Foam::IOobject(regionName, runTime.timeName(), runTime, Foam::IOobject::MUST_READ),
        false
    ));
    meshPtr->init(true); // initialise all (lower levels and current)

    return meshPtr;
}

/* @brief create a NeoN executor from a name
 * @return the Neon::Executor
 */
NeoN::Executor createExecutor(const Foam::word& execName)
{
    Foam::Info << "Creating Executor: " << execName << Foam::endl;
    if (execName == "Serial")
    {
        Foam::Info << "Serial Executor" << Foam::endl;
        return NeoN::SerialExecutor();
    }
    if (execName == "CPU")
    {
        Foam::Info << "CPU Executor" << Foam::endl;
        return NeoN::CPUExecutor();
    }
    if (execName == "GPU")
    {
        Foam::Info << "GPU Executor" << Foam::endl;
        return NeoN::GPUExecutor();
    }
    Foam::FatalError << "unknown Executor: " << execName << Foam::nl
                     << "Available executors: Serial, CPU, GPU" << Foam::nl
                     << Foam::abort(Foam::FatalError);

    return NeoN::SerialExecutor();
}

/* @brief create a NeoN executor from a dictionary
 * @return the Neon::Executor
 */
NeoN::Executor createExecutor(const Foam::dictionary& dict)
{
    auto execName = dict.get<Foam::word>("executor");
    return createExecutor(execName);
}

FoamAdapter::RunTime createAdapterRunTime(const Foam::Time& in)
{
    std::cout << __FILE__ << ":"
              << "Creating FoamAdapter runTime\n";
    auto [adjustTimeStep, maxCo, maxDeltaT] = timeControls(in);
    auto exec = createExecutor(in.controlDict());
    std::unique_ptr<MeshAdapter> meshPtr = createMesh(exec, in);
    MeshAdapter& mesh = *meshPtr;

    auto& nfMesh = mesh.nfMesh();
    return FoamAdapter::RunTime {
        .db = NeoN::Database(),
        .meshPtr = std::move(meshPtr),
        .mesh = mesh,
        .nfMesh = mesh.nfMesh(),
        .exec = exec,
        .t = in.time().value(),
        .dt = in.deltaT().value(),
        .adjustTimeStep = adjustTimeStep,
        .maxCo = maxCo,
        .maxDeltaT = maxDeltaT,
        .controlDict = convert(in.controlDict()),
        .fvSolutionDict = convert(mesh.solutionDict()),
        .fvSchemesDict = convert(mesh.schemesDict())
    };
}

}
