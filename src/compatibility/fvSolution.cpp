// SPDX-License-Identifier: GPL-3.0-or-later
//
// SPDX-FileCopyrightText: 2023 FoamAdapter authors
/* This file implements comparison operator to compare OpenFOAM and corresponding FoamAdapter fields
 * TODO the comparison operator only make sense for testing purposes
 * so this should be part of the tests
 */


#include "FoamAdapter/compatibility/fvSolution.hpp"
#include <map>
#include <NeoN/core/primitives/scalar.hpp>
#include <NeoN/core/primitives/label.hpp>


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

    std::cout << __FILE__ << __LINE__ << "1\n";
    std::string& solverName = solverDict.get<std::string>("solver");
    std::cout << __FILE__ << __LINE__ << "2\n";
    auto it = solverMap.find(solverName);
    if (it != solverMap.end())
    {
        solverName = it->second.first;
        if (solverName == "GAMG")
        {
            throw std::runtime_error(
                "GAMG is not supported in FoamAdapter, please use a different solver."
            );
        }
        solverDict.insert("type", it->second.second);
    }
}

void updatePreconditioner(NeoN::Dictionary& solverDict)
{
    // Map OpenFOAM preconditioner types to NeoN/Ginkgo preconditioner types
    static const std::map<std::string, NeoN::Dictionary> preconditionerMap = {
        {"DIC",
         NeoN::Dictionary(
             {{std::string("type"), std::string("preconditioner::Jacobi")},
              {std::string("max_block_size"), 8}}
         )},
        {"DILU",
         NeoN::Dictionary(
             {{std::string("type"), std::string("preconditioner::Ilu")},
              {std::string("reverse_apply"), false},
              {std::string("factorization"),
               NeoN::Dictionary({{std::string("type"), std::string("factorization::ParIlu")}})}}
         )},
    };


    if (solverDict.isDict("preconditioner"))
    {
        NeoN::Dictionary& preconditionerDict = solverDict.subDict("preconditioner");
        if (preconditionerDict.isDict("type"))
        {
            throw std::runtime_error(
                "GAMG is not supported in FoamAdapter, please use a different preconditioner."
            );
            // std::string& preconditionerType = preconditionerDict.get<std::string>("type");
            // auto it = preconditionerMap.find(preconditionerType);
            // if (it != preconditionerMap.end())
            // {
            //     preconditionerType = it->second;
            // }
        }
    }
    else
    {
        // If no preconditioner is specified, we can insert a default one
        std::string& preconditionerName = solverDict.get<std::string>("preconditioner");
        auto it = preconditionerMap.find(preconditionerName);
        if (it != preconditionerMap.end())
        {
            solverDict.insert("preconditioner", it->second);
        }
    }
}

void updateCriteria(NeoN::Dictionary& solverDict)
{
    // Ensure the criteria dictionary exists
    if (!solverDict.contains("criteria"))
    {
        solverDict.insert("criteria", NeoN::Dictionary());
    }

    NeoN::Dictionary& criteriaDict = solverDict.subDict("criteria");

    // set default max iteration count
    criteriaDict.insert("type", std::string("Iteration"));
    criteriaDict.insert("max_iters", 1000);

    // Set default values for relative residual norm and iteration count
    // if (solverDict.contains("relTol"))
    // {
    //     criteriaDict.insert("relative_residual_norm", solverDict.get<NeoN::scalar>("relTol"));
    // }
    // if (solverDict.contains("maxIter"))
    // {
    //     criteriaDict.insert("iteration", solverDict.get<NeoN::label>("maxIter"));
    // }
    // if (solverDict.contains("tolerance"))
    // {
    //     criteriaDict.insert("absolute_residual_norm", solverDict.get<NeoN::scalar>("tolerance"));
    // }
}


NeoN::Dictionary mapFvSolution(const NeoN::Dictionary& solverDict)
{
    NeoN::Dictionary modSolverDict = solverDict;

    updateSolver(modSolverDict);
    updatePreconditioner(modSolverDict);
    updateCriteria(modSolverDict);

    return modSolverDict;
}

} // namespace FoamAdapter
