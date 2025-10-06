// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 FoamAdapter authors
#pragma once

#include <tuple>

#include "NeoN/NeoN.hpp"
#include "FoamAdapter/datastructures/runTime.hpp"

#include "fvc.H"

namespace FoamAdapter
{

std::tuple<bool, Foam::scalar,Foam::scalar> timeControls(const Foam::Time& runTime);

void setDeltaT(Foam::Time& runTime, Foam::scalar maxCo, Foam::scalar CoNum, Foam::scalar maxDeltaT);

NeoN::Executor createExecutor(const Foam::dictionary& dict);

/* @brief create the commonly required objects for a simulation
 * @return a tuple of the executor, the controlDict, the schemesDict, the  solutionDict*/
RunTime createAdapterRunTime(const Foam::Time& runTime, const NeoN::Executor exec);

/* @brief create a struct holding the commonly required objects for a simulation
 * @return the RunTime instance*/
RunTime createAdapterRunTime(const Foam::Time& runTime);

} // namespace Foam
