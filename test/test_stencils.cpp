// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main
#include <unordered_set>
#include "NeoN/finiteVolume/cellCentred/stencil/cellToFaceStencil.hpp"


#include "gaussConvectionScheme.H"
#include "NeoN/core/input.hpp"
#include "NeoN/dsl/explicit.hpp"

#include "common.hpp"

namespace fvcc = NeoN::finiteVolume::cellCentred;
namespace dsl = NeoN::dsl;

extern Foam::Time* timePtr;    // A single time object
extern Foam::argList* argsPtr; // Some forks want argList access at createMesh.H
extern Foam::fvMesh* meshPtr;  // A single mesh object


TEST_CASE("cell To Face Stencil")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    NeoN::Executor exec = GENERATE(NeoN::Executor(NeoN::SerialExecutor {})
                                   // NeoN::Executor(NeoN::CPUExecutor {}),
                                   // NeoN::Executor(NeoN::GPUExecutor {})
    );

    std::string execName = std::visit([](auto e) { return e.name(); }, exec);

    auto meshPtr = Foam::createMesh(exec, runTime);
    Foam::MeshAdapter& mesh = *meshPtr;
    auto& nfMesh = mesh.nfMesh();

    SECTION("cellToFaceStencil_" + execName)
    {
        fvcc::CellToFaceStencil cellToFaceStencil(nfMesh);
        NeoN::SegmentedVector<NeoN::localIdx, NeoN::localIdx> stencil =
            cellToFaceStencil.computeStencil();

        auto stencilView = stencil.view();
        Foam::label nFaces = mesh.nFaces();

        forAll(mesh.cells(), celli)
        {
            std::unordered_set<Foam::label> faceSet;
            REQUIRE(stencilView.span(celli).size() == mesh.cells()[celli].size());
            for (auto facei : mesh.cells()[celli])
            {
                faceSet.insert(facei);
            }

            for (auto facei : stencilView.span(celli))
            {
                REQUIRE(faceSet.contains(facei));
            }
        }
    }
}
