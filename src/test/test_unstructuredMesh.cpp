// SPDX-License-Identifier: GPL-2.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

// Add necessary include paths
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
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

#include "NeoFOAM/mesh/unstructuredMesh/unstructuredMesh.hpp"
#include "NeoFOAM/mesh/stencil/FvccGeometryScheme.hpp"
#include "NeoFOAM/mesh/stencil/BasicFvccGeometryScheme.hpp"

#include "FoamAdapter/readers/foamMesh.hpp"
#include "FoamAdapter/writers/writers.hpp"
#include "FoamAdapter/setup/setup.hpp"
#include "FoamAdapter/fvcc/mesh/fvccNeoMesh.hpp"

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
    // #include "createMesh.H"
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
    NeoFOAM::executor exec = GENERATE(
        NeoFOAM::executor(NeoFOAM::CPUExecutor{}),
        NeoFOAM::executor(NeoFOAM::OMPExecutor{}),
        NeoFOAM::executor(NeoFOAM::GPUExecutor{}));
    // NeoFOAM::executor exec = NeoFOAM::CPUExecutor {};
    // NeoFOAM::executor exec = NeoFOAM::CPUExecutor{};
    std::string exec_name = std::visit([](auto e)
                                       { return e.print(); },
                                       exec);

    std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec, runTime);
    Foam::fvccNeoMesh &mesh = *meshPtr;
    const NeoFOAM::unstructuredMesh &uMesh = mesh.uMesh();

    SECTION("Fields" + exec_name)
    {
        const int32_t nCells = uMesh.nCells();

        REQUIRE(nCells == 9);

        const int32_t nInternalFaces = uMesh.nInternalFaces();

        REQUIRE(nInternalFaces == 12);

        const NeoFOAM::Field<NeoFOAM::Vector> &points = uMesh.points();

        REQUIRE(points.size() == 32);

        const NeoFOAM::Field<NeoFOAM::scalar> &cellVolumes = uMesh.cellVolumes();
        REQUIRE(cellVolumes.size() == 9);
        REQUIRE(NeoFOAM::sum(cellVolumes) == Catch::Approx(0.1));
        REQUIRE(cellVolumes.copyToHost().field()[0] == Catch::Approx(0.1 / 9.0));

        const NeoFOAM::Field<NeoFOAM::Vector> &cellCentres = uMesh.cellCentres();

        REQUIRE(cellCentres.size() == 9);
        auto cellCenterSum = NeoFOAM::sum(cellCentres);
        REQUIRE(cellCenterSum[0] / cellVolumes.size() == Catch::Approx(0.5));
        REQUIRE(cellCenterSum[1] / cellVolumes.size() == Catch::Approx(0.5));
        REQUIRE(cellCenterSum[2] / cellVolumes.size() == Catch::Approx(0.05));

        const NeoFOAM::Field<NeoFOAM::Vector> &faceCentres = uMesh.faceCentres();

        REQUIRE(faceCentres.size() == 42);
        auto faceCentresSum = NeoFOAM::sum(faceCentres);
        REQUIRE(faceCentresSum[0] / faceCentres.size() == Catch::Approx(0.5));
        REQUIRE(faceCentresSum[1] / faceCentres.size() == Catch::Approx(0.5));
        REQUIRE(faceCentresSum[2] / faceCentres.size() == Catch::Approx(0.05));

        const NeoFOAM::Field<NeoFOAM::Vector> &faceAreas = uMesh.faceAreas();

        REQUIRE(faceAreas.size() == 42);

        const NeoFOAM::Field<NeoFOAM::scalar> &magFaceAreas = uMesh.magFaceAreas();

        REQUIRE(magFaceAreas.size() == 42);
        REQUIRE(NeoFOAM::sum(magFaceAreas) == Catch::Approx(2.8));

        const NeoFOAM::Field<NeoFOAM::label> faceOwner = uMesh.faceOwner();

        REQUIRE(faceOwner.size() == 42);

        const NeoFOAM::Field<NeoFOAM::label> faceNeighbour = uMesh.faceNeighbour();

        REQUIRE(faceNeighbour.size() == 12);
    }

    SECTION("boundaryMesh" + exec_name)
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

        // TODO: prettify the following tests
        SECTION("faceCells")
        {
            auto faceCells = bMesh.faceCells().copyToHost().field();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
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
            const auto &Cf = bMesh.Cf().copyToHost().field();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
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
            const auto &Cn = bMesh.Cn().copyToHost().field();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
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
            const auto &Sf = bMesh.Sf().copyToHost().field();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
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
            const auto &magSf = bMesh.magSf().copyToHost().field();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
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
            const auto &nf = bMesh.nf().copyToHost().field();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
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
            const auto &delta = bMesh.delta().copyToHost().field();;
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
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
            const auto &weights = bMesh.weights().copyToHost().field();
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
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
            const auto &deltaCoeffs = bMesh.deltaCoeffs().copyToHost().field();;
            forAll(bMeshOF, patchi)
            {
                const Foam::fvPatch &patchOF = bMeshOF[patchi];
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

struct ApproxScalar
{
    Foam::scalar margin;
    bool operator()(double rhs, double lhs) const
    {
        return Catch::Approx(rhs).margin(margin) == lhs;
    }
};

TEST_CASE("fvccGeometryScheme")
{
    Foam::Time &runTime = *timePtr;
    Foam::argList &args = *argsPtr;
    NeoFOAM::executor exec = GENERATE(
        NeoFOAM::executor(NeoFOAM::CPUExecutor {}),
        NeoFOAM::executor(NeoFOAM::OMPExecutor {}), 
        NeoFOAM::executor(NeoFOAM::GPUExecutor {})
    );
    // NeoFOAM::executor exec = NeoFOAM::CPUExecutor{};
    std::string exec_name = std::visit([](auto e) { return e.print(); },exec);

    std::unique_ptr<Foam::fvccNeoMesh> meshPtr = Foam::createMesh(exec, runTime);
    Foam::fvccNeoMesh &mesh = *meshPtr;
    const NeoFOAM::unstructuredMesh& uMesh = mesh.uMesh();

    SECTION("BasicFvccGeometryScheme" + exec_name)
    {
        // update on construction
        NeoFOAM::FvccGeometryScheme scheme(exec,uMesh, std::make_unique<NeoFOAM::BasicFvccGeometryScheme>(uMesh)); 
        scheme.update(); // make sure it uptodate
        auto foam_weights = mesh.weights();

        auto weights = scheme.weights().internalField().copyToHost().field();
        std::span<Foam::scalar> s_foam_weights(foam_weights.primitiveFieldRef().data(), foam_weights.size());
        REQUIRE_THAT(weights.subspan(0, foam_weights.size()), Catch::Matchers::RangeEquals(s_foam_weights, ApproxScalar(1e-16)));
    }

    SECTION("DefaultBasicFvccGeometryScheme" + exec_name)
    {
        // update on construction
        NeoFOAM::FvccGeometryScheme scheme(uMesh); 
        scheme.update(); // make sure it uptodate
        auto foam_weights = mesh.weights();

        auto weights = scheme.weights().internalField().copyToHost().field();
        std::span<Foam::scalar> s_foam_weights(foam_weights.primitiveFieldRef().data(), foam_weights.size());
        REQUIRE_THAT(weights.subspan(0, foam_weights.size()), Catch::Matchers::RangeEquals(s_foam_weights, ApproxScalar(1e-16)));
    }

}
