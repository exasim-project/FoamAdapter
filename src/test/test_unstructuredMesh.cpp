// SPDX-License-Identifier: GPL-2.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

// Add necessary include paths
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/catch_approx.hpp>

#include "NeoFOAM/fields/Field.hpp"
#include "NeoFOAM/fields/FieldOperations.hpp"
#include "NeoFOAM/fields/FieldTypeDefs.hpp"
#include "NeoFOAM/fields/comparisions/fieldComparision.hpp"

#include "NeoFOAM/fields/boundaryFields.hpp"
#include "NeoFOAM/fields/domainField.hpp"
#include "NeoFOAM/fields/operations/sum.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/fields/fvccVolField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/fvccBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/vol/scalar/fvccScalarFixedValueBoundaryField.hpp"

#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/writers/writers.hpp"
#define namespaceFoam // Suppress <using namespace Foam;>
#include "fvCFD.H"
#include <vector>

Foam::Time *timePtr;    // A single time object
Foam::argList *argsPtr; // Some forks want argList access at createMesh.H
Foam::fvMesh *meshPtr;  // A single mesh object

int main(int argc, char *argv[])
{

    // Initialize Catch2
    Kokkos::initialize(argc, argv);
    Catch::Session session;

    // Specify command line options
    int returnCode = session.applyCommandLine(argc, argv);
    if (returnCode != 0) // Indicates a command line error
        return returnCode;

#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
    argsPtr = &args;
    timePtr = &runTime;

    int result = session.run();

    // Run benchmarks if there are any
    Kokkos::finalize();

    return result;
}

TEST_CASE("unstructuredMesh")
{
    Foam::Time &runTime = *timePtr;
    Foam::argList &args = *argsPtr;
#include "createMesh.H"
    NeoFOAM::executor exec = NeoFOAM::CPUExecutor();
    Foam::Info << "reading mesh" << Foam::endl;
    NeoFOAM::unstructuredMesh uMesh = readOpenFOAMMesh(exec, mesh);

    const int32_t nCells = uMesh.nCells();

    REQUIRE(nCells == 9);

    const int32_t nInternalFaces = uMesh.nInternalFaces();

    REQUIRE(nInternalFaces == 12);

    const NeoFOAM::Field<NeoFOAM::Vector> points = uMesh.points();

    REQUIRE(points.size() == 32);

    const NeoFOAM::Field<NeoFOAM::scalar> cellVolumes = uMesh.cellVolumes();
    REQUIRE(cellVolumes.size() == 9);
    REQUIRE(NeoFOAM::sum(cellVolumes) == Catch::Approx(0.1));
    REQUIRE(cellVolumes.field()[0] == Catch::Approx(0.1 / 9.0));

    const NeoFOAM::Field<NeoFOAM::Vector> cellCentres = uMesh.cellCentres();

    REQUIRE(cellCentres.size() == 9);
    auto cellCenterSum = NeoFOAM::sum(cellCentres);
    REQUIRE(cellCenterSum[0] / cellVolumes.size() == Catch::Approx(0.5));
    REQUIRE(cellCenterSum[1] / cellVolumes.size() == Catch::Approx(0.5));
    REQUIRE(cellCenterSum[2] / cellVolumes.size() == Catch::Approx(0.05));

    const NeoFOAM::Field<NeoFOAM::Vector> faceCentres = uMesh.faceCentres();

    REQUIRE(faceCentres.size() == 42);
    auto faceCentresSum = NeoFOAM::sum(faceCentres);
    REQUIRE(faceCentresSum[0] / faceCentres.size() == Catch::Approx(0.5));
    REQUIRE(faceCentresSum[1] / faceCentres.size() == Catch::Approx(0.5));
    REQUIRE(faceCentresSum[2] / faceCentres.size() == Catch::Approx(0.05));

    const NeoFOAM::Field<NeoFOAM::Vector> faceAreas = uMesh.faceAreas();

    REQUIRE(faceAreas.size() == 42);

    const NeoFOAM::Field<NeoFOAM::scalar> magFaceAreas = uMesh.magFaceAreas();

    REQUIRE(magFaceAreas.size() == 42);
    REQUIRE(NeoFOAM::sum(magFaceAreas) == Catch::Approx(2.8));

    const NeoFOAM::Field<NeoFOAM::label> faceOwner = uMesh.faceOwner();

    REQUIRE(faceOwner.size() == 42);

    const NeoFOAM::Field<NeoFOAM::label> faceNeighbour = uMesh.faceNeighbour();

    REQUIRE(faceNeighbour.size() == 12);

    SECTION("boundaryMesh")
    {
        const Foam::fvBoundaryMesh &bMeshOF = mesh.boundary();
        const NeoFOAM::BoundaryMesh &bMesh = uMesh.boundaryMesh();
        const auto &offset = bMesh.offset();

        SECTION("offset")
        {
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
                const std::string patchName = patchOF.name();
                REQUIRE(patchOF.size() == bMesh.faceCells(patchi).size());
            }
        }

        SECTION("faceCells")
        {
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
                const auto &faceCells = bMesh.faceCells(patchi);
                forAll(faceCells, i)
                {
                    REQUIRE(faceCells[i] == patchOF.faceCells()[i]);
                }
            }
        }

        SECTION("Cf")
        {
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
                const auto &Cf = bMesh.Cf(patchi);
                forAll(Cf, i)
                {
                    REQUIRE(Cf[i][0] == patchOF.Cf()[i][0]);
                    REQUIRE(Cf[i][1] == patchOF.Cf()[i][1]);
                    REQUIRE(Cf[i][2] == patchOF.Cf()[i][2]);
                }
            }
        }

        SECTION("Cn")
        {
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
                const auto &Cn = bMesh.Cn(patchi);
                forAll(Cn, i)
                {
                    REQUIRE(Cn[i][0] == patchOF.Cn()()[i][0]);
                    REQUIRE(Cn[i][1] == patchOF.Cn()()[i][1]);
                    REQUIRE(Cn[i][2] == patchOF.Cn()()[i][2]);
                }
            }
        }

        SECTION("Sf")
        {
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
                const auto &Sf = bMesh.Sf(patchi);
                forAll(Sf, i)
                {
                    REQUIRE(Sf[i][0] == patchOF.Sf()[i][0]);
                    REQUIRE(Sf[i][1] == patchOF.Sf()[i][1]);
                    REQUIRE(Sf[i][2] == patchOF.Sf()[i][2]);
                }
            }
        }

        SECTION("magSf")
        {
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
                const auto &magSf = bMesh.magSf(patchi);
                forAll(magSf, i)
                {
                    REQUIRE(magSf[i] == patchOF.magSf()[i]);
                }
            }
        }

        SECTION("nf")
        {
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
                const auto &nf = bMesh.nf(patchi);
                forAll(nf, i)
                {
                    REQUIRE(nf[i][0] == patchOF.nf()()[i][0]);
                    REQUIRE(nf[i][1] == patchOF.nf()()[i][1]);
                    REQUIRE(nf[i][2] == patchOF.nf()()[i][2]);
                }
            }
        }

        SECTION("delta")
        {
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
                const auto &delta = bMesh.delta(patchi);
                forAll(delta, i)
                {
                    REQUIRE(delta[i][0] == patchOF.delta()()[i][0]);
                    REQUIRE(delta[i][1] == patchOF.delta()()[i][1]);
                    REQUIRE(delta[i][2] == patchOF.delta()()[i][2]);
                }
            }
        }

        SECTION("weights")
        {
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
                const auto &weights = bMesh.weights(patchi);
                forAll(weights, i)
                {
                    REQUIRE(weights[i] == patchOF.weights()[i]);
                }
            }
        }

        SECTION("deltaCoeffs")
        {
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
                const auto &deltaCoeffs = bMesh.deltaCoeffs(patchi);
                forAll(deltaCoeffs, i)
                {
                    REQUIRE(deltaCoeffs[i] == patchOF.deltaCoeffs()[i]);
                }
            }
        }
    }
}
