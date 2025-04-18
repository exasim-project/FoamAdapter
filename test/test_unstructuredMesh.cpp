// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include <vector>

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include "catch2/common.hpp"

#include "FoamAdapter/setup.hpp"
#include "FoamAdapter/comparison.hpp"
#include "FoamAdapter/meshAdapter.hpp"


namespace fvcc = NeoN::finiteVolume::cellCentred;

extern Foam::Time* timePtr;   // A single time object
extern Foam::fvMesh* meshPtr; // A single mesh object


TEST_CASE("unstructuredMesh")
{
    NeoN::Executor exec = GENERATE(
        NeoN::Executor(NeoN::CPUExecutor {}),
        NeoN::Executor(NeoN::SerialExecutor {}),
        NeoN::Executor(NeoN::GPUExecutor {})
    );

    std::string execName = std::visit([](auto e) { return e.name(); }, exec);

    std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, *timePtr);
    const Foam::fvMesh& ofMesh = *meshPtr;
    const NeoN::UnstructuredMesh& nfMesh = meshPtr->nfMesh();

    SECTION("Internal mesh" + execName)
    {
        REQUIRE(nfMesh.nCells() == ofMesh.nCells());

        REQUIRE(nfMesh.nInternalFaces() == ofMesh.nInternalFaces());

        SECTION("points") { REQUIRE(nfMesh.points() == ofMesh.points()); }

        SECTION("cellVolumes") { REQUIRE(nfMesh.cellVolumes() == ofMesh.cellVolumes()); }

        SECTION("cellCentres") { REQUIRE(nfMesh.cellCentres() == ofMesh.cellCentres()); }

        SECTION("faceCentres") { REQUIRE(nfMesh.faceCentres() == ofMesh.faceCentres()); }

        SECTION("faceAreas") { REQUIRE(nfMesh.faceAreas() == ofMesh.faceAreas()); }

        SECTION("magFaceAreas")
        {
            Foam::scalarField magSf(mag(ofMesh.faceAreas()));
            REQUIRE(nfMesh.magFaceAreas() == magSf);
        }

        // TODO This requires ofMesh.faceOwner to be field but faceOwner() is a list
        // SECTION("faceOwner") { REQUIRE(nfMesh.faceOwner() == ofMesh.faceOwner()); }
        // SECTION("faceNeighbour") { REQUIRE(nfMesh.faceNeighbour() == ofMesh.faceNeighbour()); }
    }

    SECTION("boundaryMesh " + execName)
    {
        const Foam::fvBoundaryMesh& ofBoundaryMesh = ofMesh.boundary();
        const NeoN::BoundaryMesh& bMesh = nfMesh.boundaryMesh();
        const auto& offset = bMesh.offset();

        SECTION("offset")
        {
            forAll(ofBoundaryMesh, patchi)
            {
                REQUIRE(ofBoundaryMesh[patchi].size() == bMesh.faceCells(patchi).size());
            }
        }

        // TODO: prettify the following tests
        SECTION("faceCells")
        {
            auto faceCellsHost = bMesh.faceCells().copyToHost();
            forAll(ofBoundaryMesh, patchi)
            {
                const Foam::fvPatch& patchOF = ofBoundaryMesh[patchi];
                NeoN::label start = bMesh.offset()[patchi];
                NeoN::label end = bMesh.offset()[patchi + 1];
                auto pFaceCells = faceCellsHost.view({start, end});
                forAll(patchOF.faceCells(), i)
                {
                    REQUIRE(pFaceCells[i] == patchOF.faceCells()[i]);
                }
            }
        }

        SECTION("Cf")
        {
            const auto& cfHost = bMesh.cf().copyToHost();
            forAll(ofBoundaryMesh, patchi)
            {
                const Foam::fvPatch& patchOF = ofBoundaryMesh[patchi];
                NeoN::label start = bMesh.offset()[patchi];
                NeoN::label end = bMesh.offset()[patchi + 1];
                auto pCf = cfHost.view({start, end});
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
            forAll(ofBoundaryMesh, patchi)
            {
                const Foam::fvPatch& patchOF = ofBoundaryMesh[patchi];
                NeoN::label start = bMesh.offset()[patchi];
                NeoN::label end = bMesh.offset()[patchi + 1];
                auto pCn = cnHost.view({start, end});
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
            forAll(ofBoundaryMesh, patchi)
            {
                const Foam::fvPatch& patchOF = ofBoundaryMesh[patchi];
                NeoN::label start = bMesh.offset()[patchi];
                NeoN::label end = bMesh.offset()[patchi + 1];
                auto pSf = sFHost.view({start, end});
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
            forAll(ofBoundaryMesh, patchi)
            {
                const Foam::fvPatch& patchOF = ofBoundaryMesh[patchi];
                NeoN::label start = bMesh.offset()[patchi];
                NeoN::label end = bMesh.offset()[patchi + 1];
                auto pMagSf = magSfHost.view({start, end});
                forAll(patchOF.magSf(), i)
                {
                    REQUIRE(pMagSf[i] == patchOF.magSf()[i]);
                }
            }
        }

        SECTION("nf")
        {
            const auto& nfHost = bMesh.nf().copyToHost();
            forAll(ofBoundaryMesh, patchi)
            {
                const Foam::fvPatch& patchOF = ofBoundaryMesh[patchi];
                NeoN::label start = bMesh.offset()[patchi];
                NeoN::label end = bMesh.offset()[patchi + 1];
                auto pNf = nfHost.view({start, end});
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
            forAll(ofBoundaryMesh, patchi)
            {
                const Foam::fvPatch& patchOF = ofBoundaryMesh[patchi];
                NeoN::label start = bMesh.offset()[patchi];
                NeoN::label end = bMesh.offset()[patchi + 1];
                auto pDelta = deltaHost.view({start, end});
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
            forAll(ofBoundaryMesh, patchi)
            {
                const Foam::fvPatch& patchOF = ofBoundaryMesh[patchi];
                NeoN::label start = bMesh.offset()[patchi];
                NeoN::label end = bMesh.offset()[patchi + 1];
                auto pWeights = weightsHost.view({start, end});
                forAll(patchOF.weights(), i)
                {
                    REQUIRE(pWeights[i] == patchOF.weights()[i]);
                }
            }
        }

        SECTION("deltaCoeffs")
        {
            const auto& deltaCoeffsHost = bMesh.deltaCoeffs().copyToHost();
            forAll(ofBoundaryMesh, patchi)
            {
                const Foam::fvPatch& patchOF = ofBoundaryMesh[patchi];
                NeoN::label start = bMesh.offset()[patchi];
                NeoN::label end = bMesh.offset()[patchi + 1];
                auto pDeltaCoeffs = deltaCoeffsHost.view({start, end});
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
    NeoN::Executor exec = GENERATE(
        NeoN::Executor(NeoN::CPUExecutor {}),
        NeoN::Executor(NeoN::SerialExecutor {}),
        NeoN::Executor(NeoN::GPUExecutor {})
    );

    std::string execName = std::visit([](auto e) { return e.name(); }, exec);

    std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, *timePtr);
    Foam::MeshAdapter& mesh = *meshPtr;
    const NeoN::UnstructuredMesh& nfMesh = mesh.nfMesh();

    SECTION("BasicFvccGeometryScheme" + execName)
    {
        // update on construction
        auto scheme =
            fvcc::GeometryScheme(exec, nfMesh, std::make_unique<fvcc::BasicGeometryScheme>(nfMesh));
        scheme.update(); // make sure it uptodate
        auto foamWeights = mesh.weights();

        auto weightsHost = scheme.weights().internalVector().copyToHost();
        std::span<Foam::scalar> sFoamWeights(
            foamWeights.primitiveFieldRef().data(),
            foamWeights.size()
        );
        REQUIRE_THAT(
            weightsHost.view({0, foamWeights.size()}),
            Catch::Matchers::RangeEquals(sFoamWeights, ApproxScalar(1e-16))
        );
    }

    SECTION("DefaultBasicFvccGeometryScheme" + execName)
    {
        // update on construction
        fvcc::GeometryScheme scheme(nfMesh);
        scheme.update(); // make sure it uptodate
        auto foamWeights = mesh.weights();

        auto weightsHost = scheme.weights().internalVector().copyToHost();
        std::span<Foam::scalar> sFoamWeights(
            foamWeights.primitiveFieldRef().data(),
            foamWeights.size()
        );
        REQUIRE_THAT(
            weightsHost.view({0, foamWeights.size()}),
            Catch::Matchers::RangeEquals(sFoamWeights, ApproxScalar(1e-16))
        );
    }
}
