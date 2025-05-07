// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 FoamAdapter authors
#pragma once

#include <tuple>

#include "NeoN/NeoN.hpp"

#include "fvc.H"

namespace FoamAdapter
{

std::tuple<bool, Foam::scalar,Foam::scalar> timeControls(const Foam::Time& runTime);

void setDeltaT(Foam::Time& runTime, Foam::scalar maxCo, Foam::scalar CoNum, Foam::scalar maxDeltaT);

NeoN::Executor createExecutor(const Foam::dictionary& dict);


} // namespace Foam
