// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 FoamAdapter authors
#pragma once

#include <tuple>

#include "NeoN/NeoN.hpp"
#include "FoamAdapter/datastructures/runTime.hpp"

#include "fvc.H"

namespace FoamAdapter
{

/*@brief based on the Courant number this function synchronizes the deltaT value in both runtimes*/
void setDeltaT(Foam::Time& ofRunTime, RunTime& nfRunTime, Foam::scalar coNum);

/* @brief create a NeoN executor from a dictionary
 * @return the Neon::Executor
 */
NeoN::Executor createExecutor(const Foam::dictionary& dict);

/* @brief create a NeoN executor from its name
 * @return the Neon::Executor
 */
NeoN::Executor createExecutor(const Foam::word& execName);

/* @brief create the commonly required objects for a simulation
 * @return a tuple of the executor, the controlDict, the schemesDict, the  solutionDict*/
RunTime createAdapterRunTime(const Foam::Time& runTime, const NeoN::Executor exec);

/* @brief create a struct holding the commonly required objects for a simulation
 * @return the RunTime instance*/
RunTime createAdapterRunTime(const Foam::Time& runTime);

>>>>>>> a7c7017 (wip refactor solver)

} // namespace Foam
