// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include <cstddef>
#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main
#include <unordered_set>
#include <set>

#include "NeoFOAM/NeoFOAM.hpp"

#include "gaussConvectionScheme.H"

#include "common.hpp"

namespace fvcc = NeoFOAM::finiteVolume::cellCentred;
namespace dsl = NeoFOAM::dsl;

extern Foam::Time* timePtr;    // A single time object
extern Foam::argList* argsPtr; // Some forks want argList access at createMesh.H
extern Foam::fvMesh* meshPtr;  // A single mesh object


TEST_CASE("sparsityPattern")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    // NeoFOAM::Executor exec = GENERATE(
    //     NeoFOAM::Executor(NeoFOAM::SerialExecutor {}),
    //     NeoFOAM::Executor(NeoFOAM::CPUExecutor {}),
    //     NeoFOAM::Executor(NeoFOAM::GPUExecutor {})
    // );
    NeoFOAM::Executor exec = NeoFOAM::SerialExecutor {};

    std::string execName = std::visit([](auto e) { return e.name(); }, exec);

    auto meshPtr = Foam::createMesh(exec, runTime);
    Foam::MeshAdapter& mesh = *meshPtr;
    auto& nfMesh = mesh.nfMesh();

    SECTION("Can create linear system and SparsityPattern from mesh" + execName)
    {
        auto sp = fvcc::SparsityPattern(nfMesh);
        auto ls = NeoFOAM::la::createEmptyLinearSystem<NeoFOAM::scalar, NeoFOAM::localIdx>(sp);

        const auto colIdxs = ls.matrix().colIdxs();
        const auto rowPtrs = ls.matrix().rowPtrs();

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
