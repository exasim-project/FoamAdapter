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

/*@brief based on the Courant number this function synchronizes the deltaT value in both runtimes*/
void setDeltaT(Foam::Time& ofRunTime, RunTime& nfRunTime, Foam::scalar coNum);

NeoN::Executor createExecutor(const Foam::dictionary& dict);

RunTime createAdapterRunTime(const Foam::Time& runTime, const NeoN::Executor exec);

/* @brief create the commonly required objects for a simulation
 * @return a tuple of the executor, the controlDict, the schemesDict, the  solutionDict*/
RunTime createAdapterRunTime(const Foam::Time& runTime, const NeoN::Executor exec);

/* @brief create a struct holding the commonly required objects for a simulation
 * @return the RunTime instance*/
RunTime createAdapterRunTime(const Foam::Time& runTime);

>>>>>>> a7c7017 (wip refactor solver)

} // namespace Foam
