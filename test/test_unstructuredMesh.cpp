// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

// Add necessary include paths
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
#include <vector>

namespace fvcc = NeoFOAM::finiteVolume::cellCentred;

extern Foam::Time* timePtr;    // A single time object
extern Foam::argList* argsPtr; // Some forks want argList access at createMesh.H
extern Foam::fvMesh* meshPtr;  // A single mesh object


TEST_CASE("unstructuredMesh")
{
    NeoFOAM::Executor exec = GENERATE(
        NeoFOAM::Executor(NeoFOAM::CPUExecutor {}),
        NeoFOAM::Executor(NeoFOAM::SerialExecutor {}),
        NeoFOAM::Executor(NeoFOAM::GPUExecutor {})
    );

    std::string exec_name = std::visit([](auto e) { return e.print(); }, exec);

    std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec, *timePtr);
    Foam::fvccNeoMesh& mesh = *meshPtr;
    const NeoFOAM::UnstructuredMesh& uMesh = mesh.uMesh();

    SECTION("Fields" + exec_name)
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

    SECTION("boundaryMesh" + exec_name)
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
                auto p_face_cells = faceCells.subspan(start, end - start);
                forAll(patchOF.faceCells(), i)
                {
                    REQUIRE(p_face_cells[i] == patchOF.faceCells()[i]);
                }
            }
        }

        SECTION("Cf")
        {
            const auto& Cf = bMesh.cf().copyToHost().span();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto p_Cf = Cf.subspan(start, end - start);
                forAll(patchOF.Cf(), i)
                {
                    REQUIRE(p_Cf[i][0] == patchOF.Cf()[i][0]);
                    REQUIRE(p_Cf[i][1] == patchOF.Cf()[i][1]);
                    REQUIRE(p_Cf[i][2] == patchOF.Cf()[i][2]);
                }
            }
        }

        SECTION("Cn")
        {
            const auto& Cn = bMesh.cn().copyToHost().span();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto p_Cn = Cn.subspan(start, end - start);
                forAll(patchOF.Cn()(), i)
                {
                    REQUIRE(p_Cn[i][0] == patchOF.Cn()()[i][0]);
                    REQUIRE(p_Cn[i][1] == patchOF.Cn()()[i][1]);
                    REQUIRE(p_Cn[i][2] == patchOF.Cn()()[i][2]);
                }
            }
        }

        SECTION("Sf")
        {
            const auto& Sf = bMesh.sf().copyToHost().span();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch& patchOF = bMeshOF[patchi];
                NeoFOAM::label start = bMesh.offset()[patchi];
                NeoFOAM::label end = bMesh.offset()[patchi + 1];
                auto p_Sf = Sf.subspan(start, end - start);
                forAll(patchOF.Sf(), i)
                {
                    REQUIRE(p_Sf[i][0] == patchOF.Sf()[i][0]);
                    REQUIRE(p_Sf[i][1] == patchOF.Sf()[i][1]);
                    REQUIRE(p_Sf[i][2] == patchOF.Sf()[i][2]);
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
                auto p_magSf = magSf.subspan(start, end - start);
                forAll(patchOF.magSf(), i)
                {
                    REQUIRE(p_magSf[i] == patchOF.magSf()[i]);
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
                auto p_nf = nf.subspan(start, end - start);
                forAll(patchOF.nf()(), i)
                {
                    REQUIRE(p_nf[i][0] == patchOF.nf()()[i][0]);
                    REQUIRE(p_nf[i][1] == patchOF.nf()()[i][1]);
                    REQUIRE(p_nf[i][2] == patchOF.nf()()[i][2]);
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
                auto p_delta = delta.subspan(start, end - start);
                forAll(patchOF.delta()(), i)
                {
                    REQUIRE(p_delta[i][0] == patchOF.delta()()[i][0]);
                    REQUIRE(p_delta[i][1] == patchOF.delta()()[i][1]);
                    REQUIRE(p_delta[i][2] == patchOF.delta()()[i][2]);
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
                auto p_weights = weights.subspan(start, end - start);
                forAll(patchOF.weights(), i)
                {
                    REQUIRE(p_weights[i] == patchOF.weights()[i]);
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
                auto p_deltaCoeffs = deltaCoeffs.subspan(start, end - start);
                forAll(patchOF.deltaCoeffs(), i)
                {
                    REQUIRE(p_deltaCoeffs[i] == patchOF.deltaCoeffs()[i]);
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

    std::string exec_name = std::visit([](auto e) { return e.print(); }, exec);

    std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec, *timePtr);
    Foam::fvccNeoMesh& mesh = *meshPtr;
    const NeoFOAM::UnstructuredMesh& uMesh = mesh.uMesh();

    SECTION("BasicFvccGeometryScheme" + exec_name)
    {
        // update on construction
        auto scheme =
            fvcc::GeometryScheme(exec, uMesh, std::make_unique<fvcc::BasicGeometryScheme>(uMesh));
        scheme.update(); // make sure it uptodate
        auto foam_weights = mesh.weights();

        auto weights = scheme.weights().internalField().copyToHost().span();
        std::span<Foam::scalar> s_foam_weights(
            foam_weights.primitiveFieldRef().data(), foam_weights.size()
        );
        REQUIRE_THAT(
            weights.subspan(0, foam_weights.size()),
            Catch::Matchers::RangeEquals(s_foam_weights, ApproxScalar(1e-16))
        );
    }

    SECTION("DefaultBasicFvccGeometryScheme" + exec_name)
    {
        // update on construction
        fvcc::GeometryScheme scheme(uMesh);
        scheme.update(); // make sure it uptodate
        auto foam_weights = mesh.weights();

        auto weights = scheme.weights().internalField().copyToHost().span();
        std::span<Foam::scalar> s_foam_weights(
            foam_weights.primitiveFieldRef().data(), foam_weights.size()
        );
        REQUIRE_THAT(
            weights.subspan(0, foam_weights.size()),
            Catch::Matchers::RangeEquals(s_foam_weights, ApproxScalar(1e-16))
        );
    }
}
