#include <tuple>
#include "DiagonalMatrix.H"
#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"


std::tuple<bool, Foam::scalar, Foam::scalar> timeControls(const Foam::Time& runTime);


Foam::scalar calculateCoNum(const Foam::surfaceScalarField& phi);

void setDeltaT(Foam::Time& runTime, Foam::scalar maxCo, Foam::scalar CoNum, Foam::scalar maxDeltaT);