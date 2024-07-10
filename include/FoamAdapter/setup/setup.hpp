// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "DiagonalMatrix.H"
#include <tuple>
#define namespaceFoam // Suppress <using namespace Foam;>
#include "FoamAdapter/fvcc/mesh/fvccNeoMesh.hpp"
#include "fvCFD.H"
#include <memory>
#include "NeoFOAM/core/executor/executor.hpp"

namespace Foam
{

std::tuple<bool, Foam::scalar, Foam::scalar> timeControls(const Foam::Time& runTime);

Foam::scalar calculateCoNum(const Foam::surfaceScalarField& phi);

void setDeltaT(Foam::Time& runTime, Foam::scalar maxCo, Foam::scalar CoNum, Foam::scalar maxDeltaT);

std::unique_ptr<Foam::fvccNeoMesh>
createMesh(const NeoFOAM::executor& exec, const Foam::Time& runTime);

std::unique_ptr<Foam::fvMesh> createMesh(const Foam::Time& runTime);

NeoFOAM::executor createExecutor(const Foam::dictionary& dict);

} // namespace Foam
