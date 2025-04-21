// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include <cstddef>
#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main
#include <unordered_set>
#include <set>

#include "common.hpp"

#include "gaussConvectionScheme.H"


namespace fvcc = NeoN::finiteVolume::cellCentred;
namespace dsl = NeoN::dsl;

extern Foam::Time* timePtr;    // A single time object
extern Foam::argList* argsPtr; // Some forks want argList access at createMesh.H
extern Foam::fvMesh* meshPtr;  // A single mesh object


TEST_CASE("sparsityPattern")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    // NeoN::Executor exec = GENERATE(
    //     NeoN::Executor(NeoN::SerialExecutor {}),
    //     NeoN::Executor(NeoN::CPUExecutor {}),
    //     NeoN::Executor(NeoN::GPUExecutor {})
    // );
    NeoN::Executor exec = NeoN::SerialExecutor {};

    std::string execName = std::visit([](auto e) { return e.name(); }, exec);

    auto meshPtr = Foam::createMesh(exec, runTime);
    Foam::MeshAdapter& mesh = *meshPtr;
    auto& nfMesh = mesh.nfMesh();

    SECTION("sparsityPattern_" + execName)
    {
        fvcc::SparsityPattern pattern(nfMesh);
        const auto colIdxs = pattern.colIdxs().view();
        const auto rowPtrs = pattern.rowPtrs().view();

        forAll(mesh.cellCells(), celli)
        {
            std::set<Foam::label> stencilCells;
            stencilCells.insert(celli);
            for (auto neiCelli : mesh.cellCells()[celli])
            {
                stencilCells.insert(neiCelli);
            }

            for (std::size_t i = rowPtrs[celli]; i < rowPtrs[celli + 1]; i++)
            {
                auto colIdx = colIdxs[i];
                REQUIRE(stencilCells.contains(colIdx));
            }
        }
    }
}
