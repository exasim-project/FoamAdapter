// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include <vector>

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include "catch2/common.hpp"

#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/writers/writers.hpp"
#include "FoamAdapter/setup/setup.hpp"
#include "FoamAdapter/fvcc/mesh/fvccNeoMesh.hpp"

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

    std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec, *timePtr);
    Foam::fvccNeoMesh& mesh = *meshPtr;
    const NeoFOAM::UnstructuredMesh& uMesh = mesh.uMesh();

    SECTION("Fields" + execName)
    {
        const int32_t nCells = uMesh.nCells();

        REQUIRE(nCells == mesh.nCells());

        const int32_t nInternalFaces = uMesh.nInternalFaces();

        REQUIRE(nInternalFaces == mesh.nInternalFaces());

        SECTION("points")
        {
            const auto points = uMesh.points().copyToHost().span();
            REQUIRE(uMesh.points().size() == mesh.points().size());
            for (int i = 0; i < points.size(); i++)
            {
                REQUIRE(points[i] == convert(mesh.points()[i]));
            }
        }


        SECTION("cellVolumes")
        {
            const auto cellVolumes = uMesh.cellVolumes().copyToHost().span();
            REQUIRE(cellVolumes.size() == mesh.cellVolumes().size());
            for (int i = 0; i < cellVolumes.size(); i++)
            {
                REQUIRE(cellVolumes[i] == mesh.cellVolumes()[i]);
            }
        }

        SECTION("cellCentres")
        {
            const auto cellCentres = uMesh.cellCentres().copyToHost().span();
            REQUIRE(cellCentres.size() == mesh.cellCentres().size());
            for (int i = 0; i < cellCentres.size(); i++)
            {
                REQUIRE(cellCentres[i] == convert(mesh.cellCentres()[i]));
            }
        }

        SECTION("faceCentres")
        {
            const auto faceCentres = uMesh.faceCentres().copyToHost().span();
            REQUIRE(faceCentres.size() == mesh.faceCentres().size());
            for (int i = 0; i < faceCentres.size(); i++)
            {
                REQUIRE(faceCentres[i] == convert(mesh.faceCentres()[i]));
            }
        }

        SECTION("faceAreas")
        {
            const auto faceAreas = uMesh.faceAreas().copyToHost().span();
            // REQUIRE(faceAreas.size() == mesh.Sf().size());
            for (int i = 0; i < mesh.Sf().size(); i++)
            {
                REQUIRE(faceAreas[i] == convert(mesh.Sf()[i]));
            }
        }


        SECTION("magFaceAreas")
        {
            Foam::scalarField magSf(mag(mesh.faceAreas()));
            const auto magFaceAreas = uMesh.magFaceAreas().copyToHost().span();
            REQUIRE(magFaceAreas.size() == magSf.size());
            for (int i = 0; i < magFaceAreas.size(); i++)
            {
                REQUIRE(magFaceAreas[i] == magSf[i]);
            }
        }


        SECTION("faceOwner")
        {
            const auto faceOwner = uMesh.faceOwner().copyToHost().span();
            REQUIRE(faceOwner.size() == mesh.faceOwner().size());
            for (int i = 0; i < faceOwner.size(); i++)
            {
                REQUIRE(faceOwner[i] == mesh.faceOwner()[i]);
            }
        }

        SECTION("faceNeighbour")
        {
            const auto faceNeighbour = uMesh.faceNeighbour().copyToHost().span();
            REQUIRE(faceNeighbour.size() == mesh.faceNeighbour().size());
            for (int i = 0; i < faceNeighbour.size(); i++)
            {
                REQUIRE(faceNeighbour[i] == mesh.faceNeighbour()[i]);
            }
        }
    }

    SECTION("boundaryMesh" + execName)
    {
        const Foam::fvBoundaryMesh& bMeshOF = mesh.boundary();
        const NeoFOAM::BoundaryMesh& bMesh = uMesh.boundaryMesh();
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
            auto faceCells = bMesh.faceCells().copyToHost().span();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pFaceCells = faceCells.subspan(start, end - start);
                forAll(patchOF.faceCells(), i)
                {
                    REQUIRE(pFaceCells[i] == patchOF.faceCells()[i]);
                }
            }
        }

        SECTION("Cf")
        {
            const auto& cf = bMesh.cf().copyToHost().span();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pCf = cf.subspan(start, end - start);
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
            const auto& cn = bMesh.cn().copyToHost().span();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pCn = cn.subspan(start, end - start);
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
            const auto& sF = bMesh.sf().copyToHost().span();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto psF = sF.subspan(start, end - start);
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
            const auto& magSf = bMesh.magSf().copyToHost().span();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pMagSf = magSf.subspan(start, end - start);
                forAll(patchOF.magSf(), i)
                {
                    REQUIRE(pMagSf[i] == patchOF.magSf()[i]);
                }
            }
        }

        SECTION("nf")
        {
            const auto& nf = bMesh.nf().copyToHost().span();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pNf = nf.subspan(start, end - start);
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
            const auto& delta = bMesh.delta().copyToHost().span();
            ;
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pDelta = delta.subspan(start, end - start);
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
            const auto& weights = bMesh.weights().copyToHost().span();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pWeights = weights.subspan(start, end - start);
                forAll(patchOF.weights(), i)
                {
                    REQUIRE(pWeights[i] == patchOF.weights()[i]);
                }
            }
        }

        SECTION("deltaCoeffs")
        {
            const auto& deltaCoeffs = bMesh.deltaCoeffs().copyToHost().span();
            ;
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto pDeltaCoeffs = deltaCoeffs.subspan(start, end - start);
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

    std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec, *timePtr);
    Foam::fvccNeoMesh& mesh = *meshPtr;
    const NeoFOAM::UnstructuredMesh& uMesh = mesh.uMesh();

    SECTION("BasicFvccGeometryScheme" + execName)
    {
        // update on construction
        auto scheme =
            fvcc::GeometryScheme(exec, uMesh, std::make_unique<fvcc::BasicGeometryScheme>(uMesh));
        scheme.update(); // make sure it uptodate
        auto foamWeights = mesh.weights();

        auto weights = scheme.weights().internalField().copyToHost().span();
        std::span<Foam::scalar> sFoamWeights(
            foamWeights.primitiveFieldRef().data(), foamWeights.size()
        );
        REQUIRE_THAT(
            weights.subspan(0, foamWeights.size()),
            Catch::Matchers::RangeEquals(sFoamWeights, ApproxScalar(1e-16))
        );
    }

    SECTION("DefaultBasicFvccGeometryScheme" + execName)
    {
        // update on construction
        fvcc::GeometryScheme scheme(uMesh);
        scheme.update(); // make sure it uptodate
        auto foamWeights = mesh.weights();

        auto weights = scheme.weights().internalField().copyToHost().span();
        std::span<Foam::scalar> sFoamWeights(
            foamWeights.primitiveFieldRef().data(), foamWeights.size()
        );
        REQUIRE_THAT(
            weights.subspan(0, foamWeights.size()),
            Catch::Matchers::RangeEquals(sFoamWeights, ApproxScalar(1e-16))
        );
    }
}
