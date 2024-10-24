// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include <vector>

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include "catch2/common.hpp"

#include "FoamAdapter/meshAdapter.hpp"
#include "FoamAdapter/writers/writers.hpp"
#include "FoamAdapter/setup.hpp"

#define namespaceFoam // Suppress <using namespace Foam;>

namespace fvcc = NeoFOAM::finiteVolume::cellCentred;

extern Foam::Time* timePtr;   // A single time object
extern Foam::fvMesh* meshPtr; // A single mesh object


TEST_CASE("unstructuredMesh")
{
    NeoFOAM::Executor exec = GENERATE(
        NeoFOAM::Executor(NeoFOAM::CPUExecutor {}),
        NeoFOAM::Executor(NeoFOAM::SerialExecutor {}),
        NeoFOAM::Executor(NeoFOAM::GPUExecutor {})
    );

    std::string execName = std::visit([](auto e) { return e.print(); }, exec);

    std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, *timePtr);
    Foam::MeshAdapter& mesh = *meshPtr;
    const NeoFOAM::UnstructuredMesh& nfMesh = mesh.nfMesh();

    SECTION("Fields" + execName)
    {
        const int32_t nCells = nfMesh.nCells();

        REQUIRE(nCells == mesh.nCells());

        const int32_t nInternalFaces = nfMesh.nInternalFaces();

        REQUIRE(nInternalFaces == mesh.nInternalFaces());

        SECTION("points") { checkField(nfMesh.points(), mesh.points()); }

        SECTION("cellVolumes") { checkField(nfMesh.cellVolumes(), mesh.cellVolumes()); }

        SECTION("cellCentres") { checkField(nfMesh.cellCentres(), mesh.cellCentres()); }

        SECTION("faceCentres") { checkField(nfMesh.faceCentres(), mesh.faceCentres()); }

        SECTION("faceAreas") { checkField(nfMesh.faceAreas(), mesh.faceAreas()); }

        SECTION("magFaceAreas")
        {
            Foam::scalarField magSf(mag(mesh.faceAreas()));
            checkField(nfMesh.magFaceAreas(), magSf);
        }

        SECTION("faceOwner") { checkField(nfMesh.faceOwner(), mesh.faceOwner()); }

        SECTION("faceNeighbour") { checkField(nfMesh.faceNeighbour(), mesh.faceNeighbour()); }
    }

    SECTION("boundaryMesh" + execName)
    {
        const Foam::fvBoundaryMesh& bMeshOF = mesh.boundary();
        const NeoFOAM::BoundaryMesh& bMesh = nfMesh.boundaryMesh();
        const auto& offset = bMesh.offset();

        SECTION("offset")
        {
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                const std::string patchName = patchOF.name();
                REQUIRE(patchOF.size() == bMesh.faceCells(patchi).size());
            }
        }

        // TODO: prettify the following tests
        SECTION("faceCells")
        {
            auto faceCellsHost = bMesh.faceCells().copyToHost();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pFaceCells = faceCellsHost.span({start, end});
                forAll(patchOF.faceCells(), i)
                {
                    REQUIRE(pFaceCells[i] == patchOF.faceCells()[i]);
                }
            }
        }

        SECTION("Cf")
        {
            const auto& cfHost = bMesh.cf().copyToHost();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pCf = cfHost.span({start, end});
                forAll(patchOF.Cf(), i)
                {
                    REQUIRE(pCf[i][0] == patchOF.Cf()[i][0]);
                    REQUIRE(pCf[i][1] == patchOF.Cf()[i][1]);
                    REQUIRE(pCf[i][2] == patchOF.Cf()[i][2]);
                }
            }
        }

        SECTION("Cn")
        {
            const auto& cnHost = bMesh.cn().copyToHost();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pCn = cnHost.span({start, end});
                forAll(patchOF.Cn()(), i)
                {
                    REQUIRE(pCn[i][0] == patchOF.Cn()()[i][0]);
                    REQUIRE(pCn[i][1] == patchOF.Cn()()[i][1]);
                    REQUIRE(pCn[i][2] == patchOF.Cn()()[i][2]);
                }
            }
        }

        SECTION("Sf")
        {
            const auto& sFHost = bMesh.sf().copyToHost();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pSf = sFHost.span({start, end});
                forAll(patchOF.Sf(), i)
                {
                    REQUIRE(pSf[i][0] == patchOF.Sf()[i][0]);
                    REQUIRE(pSf[i][1] == patchOF.Sf()[i][1]);
                    REQUIRE(pSf[i][2] == patchOF.Sf()[i][2]);
                }
            }
        }

        SECTION("magSf")
        {
            const auto& magSfHost = bMesh.magSf().copyToHost();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pMagSf = magSfHost.span({start, end});
                forAll(patchOF.magSf(), i)
                {
                    REQUIRE(pMagSf[i] == patchOF.magSf()[i]);
                }
            }
        }

        SECTION("nf")
        {
            const auto& nfHost = bMesh.nf().copyToHost();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pNf = nfHost.span({start, end});
                forAll(patchOF.nf()(), i)
                {
                    REQUIRE(pNf[i][0] == patchOF.nf()()[i][0]);
                    REQUIRE(pNf[i][1] == patchOF.nf()()[i][1]);
                    REQUIRE(pNf[i][2] == patchOF.nf()()[i][2]);
                }
            }
        }

        SECTION("delta")
        {
            const auto& deltaHost = bMesh.delta().copyToHost();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pDelta = deltaHost.span({start, end});
                forAll(patchOF.delta()(), i)
                {
                    REQUIRE(pDelta[i][0] == patchOF.delta()()[i][0]);
                    REQUIRE(pDelta[i][1] == patchOF.delta()()[i][1]);
                    REQUIRE(pDelta[i][2] == patchOF.delta()()[i][2]);
                }
            }
        }

        SECTION("weights")
        {
            const auto& weightsHost = bMesh.weights().copyToHost();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pWeights = weightsHost.span({start, end});
                forAll(patchOF.weights(), i)
                {
                    REQUIRE(pWeights[i] == patchOF.weights()[i]);
                }
            }
        }

        SECTION("deltaCoeffs")
        {
            const auto& deltaCoeffsHost = bMesh.deltaCoeffs().copyToHost();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pDeltaCoeffs = deltaCoeffsHost.span({start, end});
                forAll(patchOF.deltaCoeffs(), i)
                {
                    REQUIRE(pDeltaCoeffs[i] == patchOF.deltaCoeffs()[i]);
                }
            }
        }
    }
}


TEST_CASE("fvccGeometryScheme")
{
    NeoFOAM::Executor exec = GENERATE(
        NeoFOAM::Executor(NeoFOAM::CPUExecutor {}),
        NeoFOAM::Executor(NeoFOAM::SerialExecutor {}),
        NeoFOAM::Executor(NeoFOAM::GPUExecutor {})
    );

    std::string execName = std::visit([](auto e) { return e.print(); }, exec);

    std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, *timePtr);
    Foam::MeshAdapter& mesh = *meshPtr;
    const NeoFOAM::UnstructuredMesh& nfMesh = mesh.nfMesh();

    SECTION("BasicFvccGeometryScheme" + execName)
    {
        // update on construction
        auto scheme =
            fvcc::GeometryScheme(exec, nfMesh, std::make_unique<fvcc::BasicGeometryScheme>(nfMesh));
        scheme.update(); // make sure it uptodate
        auto foamWeights = mesh.weights();

        auto weightsHost = scheme.weights().internalField().copyToHost();
        std::span<Foam::scalar> sFoamWeights(
            foamWeights.primitiveFieldRef().data(), foamWeights.size()
        );
        REQUIRE_THAT(
            weightsHost.span({0, foamWeights.size()}),
            Catch::Matchers::RangeEquals(sFoamWeights, ApproxScalar(1e-16))
        );
    }

    SECTION("DefaultBasicFvccGeometryScheme" + execName)
    {
        // update on construction
        fvcc::GeometryScheme scheme(nfMesh);
        scheme.update(); // make sure it uptodate
        auto foamWeights = mesh.weights();

        auto weightsHost = scheme.weights().internalField().copyToHost();
        std::span<Foam::scalar> sFoamWeights(
            foamWeights.primitiveFieldRef().data(), foamWeights.size()
        );
        REQUIRE_THAT(
            weightsHost.span({0, foamWeights.size()}),
            Catch::Matchers::RangeEquals(sFoamWeights, ApproxScalar(1e-16))
        );
    }
}
