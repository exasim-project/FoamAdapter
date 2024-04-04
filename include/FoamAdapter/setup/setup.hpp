// SPDX-License-Identifier: GPL-3.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include <tuple>
#include "DiagonalMatrix.H"
#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"
#include "FoamAdapter/fvcc/mesh/fvccNeoMesh.hpp"
#include <memory>


namespace Foam
{
std::tuple<bool, Foam::scalar, Foam::scalar> timeControls(const Foam::Time& runTime);


Foam::scalar calculateCoNum(const Foam::surfaceScalarField& phi);

void setDeltaT(Foam::Time& runTime, Foam::scalar maxCo, Foam::scalar CoNum, Foam::scalar maxDeltaT);

std::unique_ptr<Foam::fvccNeoMesh> createMesh(const NeoFOAM::executor& exec, const Foam::IOobject& io);

} // namespace Foam
