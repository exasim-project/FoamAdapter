// SPDX-License-Identifier: GPL-3.0-or-later
//
// SPDX-FileCopyrightText: 2023 FoamAdapter authors
/* This file implements comparison operator to compare OpenFOAM and corresponding FoamAdapter fields
 * TODO the comparison operator only make sense for testing purposes
 * so this should be part of the tests
 */


#include "FoamAdapter/compatibility/fvSolution.hpp"
#include <map>


namespace FoamAdapter
{

void updateSolver(NeoN::Dictionary& solverDict)
{
    // Map OpenFOAM solver names to NeoN/Ginkgo solver names and types
    static const std::map<std::string, std::pair<std::string, std::string>> solverMap = {
        {"PCG", {"Ginkgo", "solver::Cg"}},
        {"PBiCG", {"Ginkgo", "solver::Bicg"}},
        {"PBiCGStab", {"Ginkgo", "solver::Bicgstab"}},
        {"GAMG", {"Ginkgo", "solver::Multigrid"}},
    };

    std::string& solverName = solverDict.get<std::string>("solver");
    auto it = solverMap.find(solverName);
    if (it != solverMap.end())
    {
        solverName = it->second.first;
        solverDict.insert("type", it->second.second);
    }
}

void updatePreconditioner(NeoN::Dictionary& preconditionerDict)
{
    // Map OpenFOAM preconditioner types to NeoN/Ginkgo preconditioner types
    static const std::map<std::string, std::string> preconditionerMap = {
        {"DIC", "Diagonal"},
        {"ILU", "IncompleteLU"},
    };

    if (preconditionerDict.isDict("type"))
    {
        std::string& preconditionerType = preconditionerDict.get<std::string>("type");
        auto it = preconditionerMap.find(preconditionerType);
        if (it != preconditionerMap.end())
        {
            preconditionerType = it->second;
        }
    }
}


NeoN::Dictionary mapFvSolution(const NeoN::Dictionary& solverDict)
{
    NeoN::Dictionary modSolverDict = solverDict;

    // std::cout << "Mapping solver: " << solverName << std::endl;
    auto keys = modSolverDict.keys();

    for (const auto& key : keys)
    {
        std::cout << "Key: " << key << std::endl;
    }
    updateSolver(modSolverDict);

    if (!modSolverDict.contains("criteria"))
    {
        modSolverDict.insert("criteria", NeoN::Dictionary());
        NeoN::Dictionary& criteriaDict = modSolverDict.get<NeoN::Dictionary>("criteria");
        criteriaDict.insert("relative_residual_norm", 1e-6);
        criteriaDict.insert("iteration", 1000);
    }
    // {
    //     NeoN::Dictionary& preconditionerDict = modSolverDict.subDict("preconditioner");
    //     if (preconditionerDict.isDict("type"))
    //     {
    //         std::string preconditionerType = preconditionerDict.get<std::string>("type");
    //         if (preconditionerType == "DIC")
    //         {
    //             preconditionerDict.insert("type", "Diagonal");
    //         }
    //         else if (preconditionerType == "ILU")
    //         {
    //             preconditionerDict.insert("type", "IncompleteLU");
    //         }
    //     }
    // }
    // if (solverDict.has(solverName))
    // {
    //     return solverDict.get<NeoN::Dictionary>(solverName);
    // }
    // else
    // {
    //     throw std::runtime_error("Solver " + solverName + " not found in the dictionary.");
    // }

    return modSolverDict;
}

} // namespace FoamAdapter
