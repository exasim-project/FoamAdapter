// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include <tuple>
#include <memory>

#define namespaceFoam // Suppress <using namespace Foam;>

#include "NeoFOAM/core/executor/executor.hpp"

#include "fvCFD.H" // include after NeoFOAM to avoid ambiguous sqrt error

#include "FoamAdapter/meshAdapter.hpp"

namespace Foam
{

std::tuple<bool, scalar, scalar> timeControls(const Time& runTime);

scalar calculateCoNum(const surfaceScalarField& phi);

void setDeltaT(Time& runTime, scalar maxCo, scalar CoNum, scalar maxDeltaT);

std::unique_ptr<fvccNeoMesh> createMesh(const NeoFOAM::Executor& exec, const Time& runTime);

std::unique_ptr<fvMesh> createMesh(const Time& runTime);

NeoFOAM::Executor createExecutor(const dictionary& dict);


} // namespace Foam
