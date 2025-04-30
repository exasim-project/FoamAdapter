// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 FoamAdapter authors
#pragma once

#include <tuple>
#include <memory>

#include "FoamAdapter/meshAdapter.hpp"

#include "NeoN/NeoN.hpp"


namespace Foam
{

std::tuple<bool, scalar, scalar> timeControls(const Time& runTime);

scalar calculateCoNum(const surfaceScalarField& phi);

void setDeltaT(Time& runTime, scalar maxCo, scalar CoNum, scalar maxDeltaT);

std::unique_ptr<MeshAdapter> createMesh(const NeoN::Executor& exec, const Time& runTime);

std::unique_ptr<fvMesh> createMesh(const Time& runTime);

NeoN::Executor createExecutor(const dictionary& dict);


} // namespace Foam
