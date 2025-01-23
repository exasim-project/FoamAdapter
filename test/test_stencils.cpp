// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main
#include <unordered_set>
#include "NeoFOAM/finiteVolume/cellCentred/stencil/cellToFaceStencil.hpp"

#define namespaceFoam // Suppress <using namespace Foam;>
#include "gaussConvectionScheme.H"
#include "NeoFOAM/core/input.hpp"
#include "NeoFOAM/dsl/explicit.hpp"

#include "common.hpp"

namespace fvcc = NeoFOAM::finiteVolume::cellCentred;
namespace dsl = NeoFOAM::dsl;

extern Foam::Time* timePtr;    // A single time object
extern Foam::argList* argsPtr; // Some forks want argList access at createMesh.H
extern Foam::fvMesh* meshPtr;  // A single mesh object


TEST_CASE("cell To Face Stencil")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    NeoFOAM::Executor exec = GENERATE(NeoFOAM::Executor(NeoFOAM::SerialExecutor {})
                                      // NeoFOAM::Executor(NeoFOAM::CPUExecutor {}),
                                      // NeoFOAM::Executor(NeoFOAM::GPUExecutor {})
    );

    std::string execName = std::visit([](auto e) { return e.name(); }, exec);

    auto meshPtr = Foam::createMesh(exec, runTime);
    Foam::MeshAdapter& mesh = *meshPtr;
    auto& nfMesh = mesh.nfMesh();

    SECTION("cellToFaceStencil_" + execName)
    {
        fvcc::CellToFaceStencil cellToFaceStencil(nfMesh);
        NeoFOAM::SegmentedField<NeoFOAM::localIdx, NeoFOAM::localIdx> stencil =
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
