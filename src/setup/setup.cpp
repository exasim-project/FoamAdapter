// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "FoamAdapter/setup/setup.hpp"
#include "error.H"


std::tuple<bool, Foam::scalar, Foam::scalar> Foam::timeControls(const Foam::Time& runTime)
{
    bool adjustTimeStep = runTime.controlDict().getOrDefault("adjustTimeStep", false);

    Foam::scalar maxCo = runTime.controlDict().getOrDefault<Foam::scalar>("maxCo", 1);

    Foam::scalar maxDeltaT =
        runTime.controlDict().getOrDefault<Foam::scalar>("maxDeltaT", Foam::GREAT);

    return std::make_tuple(adjustTimeStep, maxCo, maxDeltaT);
}


Foam::scalar Foam::calculateCoNum(const Foam::surfaceScalarField& phi)
{
    const Foam::fvMesh& mesh = phi.mesh();
    const Foam::Time& runTime = mesh.time();
    Foam::scalarField sumPhi(Foam::fvc::surfaceSum(mag(phi))().primitiveField());
    Foam::scalar CoNum = 0.5 * Foam::gMax(sumPhi / mesh.V().field()) * runTime.deltaTValue();
    Foam::scalar meanCoNum = 0.5 * (gSum(sumPhi) / gSum(mesh.V().field())) * runTime.deltaTValue();

    Foam::Info << "Courant Number mean: " << meanCoNum << " max: " << CoNum << Foam::endl;
    return CoNum;
}

void Foam::setDeltaT(
    Foam::Time& runTime, Foam::scalar maxCo, Foam::scalar CoNum, Foam::scalar maxDeltaT
)
{
    Foam::scalar maxDeltaTFact = maxCo / (CoNum + Foam::SMALL);
    Foam::scalar deltaTFact = Foam::min(Foam::min(maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact), 1.2);

    runTime.setDeltaT(Foam::min(deltaTFact * runTime.deltaTValue(), maxDeltaT));

    Foam::Info << "deltaT = " << runTime.deltaTValue() << Foam::endl;
}


std::unique_ptr<Foam::fvccNeoMesh>
Foam::createMesh(const NeoFOAM::Executor& exec, const Foam::Time& runTime)
{
    Foam::word regionName(Foam::polyMesh::defaultRegion);
    Foam::IOobject io(regionName, runTime.timeName(), runTime, Foam::IOobject::MUST_READ);
    return std::make_unique<Foam::fvccNeoMesh>(exec, io);
}

std::unique_ptr<Foam::fvMesh> Foam::createMesh(const Foam::Time& runTime)
{
    std::unique_ptr<Foam::fvMesh> meshPtr;
    Foam::Info << "Create mesh";
    Foam::word regionName(Foam::polyMesh::defaultRegion);
    Foam::Info << " for time = " << runTime.timeName() << Foam::nl;

    meshPtr.reset(new Foam::fvMesh(
        Foam::IOobject(regionName, runTime.timeName(), runTime, Foam::IOobject::MUST_READ), false
    ));
    meshPtr->init(true); // initialise all (lower levels and current)

    return meshPtr;
}

NeoFOAM::Executor Foam::createExecutor(const Foam::dictionary& dict)
{
    auto exec_name = dict.get<Foam::word>("executor");
    Foam::Info << "Creating Executor: " << exec_name << Foam::endl;
    if (exec_name == "Serial")
    {
        Foam::Info << "Serial Executor" << Foam::endl;
        return NeoFOAM::SerialExecutor();
    }
    if (exec_name == "CPU")
    {
        Foam::Info << "CPU Executor" << Foam::endl;
        return NeoFOAM::CPUExecutor();
    }
    if (exec_name == "GPU")
    {
        Foam::Info << "GPU Executor" << Foam::endl;
        return NeoFOAM::GPUExecutor();
    }
    Foam::FatalError << "unknown Executor: " << exec_name << Foam::nl
                     << "Available executors: Serial, CPU, GPU" << Foam::nl
                     << Foam::abort(Foam::FatalError);

    return NeoFOAM::SerialExecutor();
}
