// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include <tuple>
#include <memory>

#define namespaceFoam // Suppress <using namespace Foam;>

#include "FoamAdapter/fvcc/mesh/fvccNeoMesh.hpp"

#include "NeoFOAM/core/executor/executor.hpp"

#include "fvCFD.H" // include after NeoFOAM to avoid ambiguous sqrt error


namespace Foam
{

std::tuple<bool, Foam::scalar, Foam::scalar> timeControls(const Foam::Time& runTime);

Foam::scalar calculateCoNum(const Foam::surfaceScalarField& phi);

void setDeltaT(Foam::Time& runTime, Foam::scalar maxCo, Foam::scalar CoNum, Foam::scalar maxDeltaT);

std::unique_ptr<Foam::fvccNeoMesh>
createMesh(const NeoFOAM::Executor& exec, const Foam::Time& runTime);

std::unique_ptr<Foam::fvMesh> createMesh(const Foam::Time& runTime);

NeoFOAM::Executor createExecutor(const Foam::dictionary& dict);

} // namespace Foam
