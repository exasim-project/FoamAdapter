#include <tuple>
#include "DiagonalMatrix.H"
#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"


std::tuple<bool, Foam::scalar, Foam::scalar> createTimeControls(const Foam::Time& runTime);

void updateTimeControls(const Foam::Time& runTime, bool& adjustTimeStep, Foam::scalar& maxCo, Foam::scalar& maxDeltaT);

Foam::scalar calculateCoNum(const Foam::surfaceScalarField& phi);

void setDeltaT(Foam::Time& runTime, Foam::scalar maxCo, Foam::scalar CoNum, Foam::scalar maxDeltaT);