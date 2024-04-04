#include "FoamAdapter/setup/setup.hpp"


std::tuple<bool, Foam::scalar, Foam::scalar> Foam::timeControls(const Foam::Time& runTime) {
    bool adjustTimeStep =
    runTime.controlDict().getOrDefault("adjustTimeStep", false);

    Foam::scalar maxCo =
        runTime.controlDict().getOrDefault<Foam::scalar>("maxCo", 1);

    Foam::scalar maxDeltaT =
        runTime.controlDict().getOrDefault<Foam::scalar>("maxDeltaT",Foam::GREAT);

    return std::make_tuple(adjustTimeStep, maxCo, maxDeltaT);
}



Foam::scalar Foam::calculateCoNum(const Foam::surfaceScalarField& phi) {
    const Foam::fvMesh& mesh = phi.mesh();
    const Foam::Time& runTime = mesh.time();
    Foam::scalarField sumPhi(Foam::fvc::surfaceSum(mag(phi))().primitiveField());
    Foam::scalar CoNum = 0.5 * Foam::gMax(sumPhi / mesh.V().field()) * runTime.deltaTValue();
    Foam::scalar meanCoNum = 0.5 * (gSum(sumPhi) / gSum(mesh.V().field())) * runTime.deltaTValue();
    
    Foam::Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum << Foam::endl;
    return CoNum;
}

void Foam::setDeltaT(Foam::Time& runTime, Foam::scalar maxCo, Foam::scalar CoNum, Foam::scalar maxDeltaT)
{
    Foam::scalar maxDeltaTFact = maxCo/(CoNum + Foam::SMALL);
    Foam::scalar deltaTFact = Foam::min(Foam::min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    runTime.setDeltaT
    (
        Foam::min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Foam::Info<< "deltaT = " <<  runTime.deltaTValue() << Foam::endl;   
}


std::unique_ptr<Foam::fvccNeoMesh> Foam::createMesh(const NeoFOAM::executor& exec, const Foam::IOobject& io)
{
    return std::make_unique<Foam::fvccNeoMesh>(exec, io);
}
