// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "FoamAdapter/meshAdapter.hpp"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{

template<typename FieldT>
FieldT flatBCField(const fvMesh& mesh, std::function<FieldT(const fvPatch&)> f)
{
    FieldT result(computeNBoundaryFaces(mesh));
    const fvBoundaryMesh& bMesh = mesh.boundary();
    label idx = 0;
    forAll(bMesh, patchI)
    {
        const fvPatch& patch = bMesh[patchI];
        auto pResult = f(patch);
        forAll(pResult, i)
        {
            result[idx] = pResult[i];
            idx++;
        }
    }
    return result;
}

defineTypeNameAndDebug(MeshAdapter, 0);

std::vector<NeoFOAM::localIdx> computeOffset(const fvMesh& mesh)
{
    std::vector<NeoFOAM::localIdx> result;
    const fvBoundaryMesh& bMesh = mesh.boundary();
    result.push_back(0);
    forAll(bMesh, patchI)
    {
        NeoFOAM::localIdx curOffset = result.back();
        const fvPatch& patch = bMesh[patchI];
        result.push_back(curOffset + patch.size());
    }
    return result;
}

int32_t computeNBoundaryFaces(const fvMesh& mesh)
{
    const fvBoundaryMesh& bMesh = mesh.boundary();
    int32_t nBoundaryFaces = 0;
    forAll(bMesh, patchI)
    {
        const fvPatch& patch = bMesh[patchI];
        nBoundaryFaces += patch.size();
    }
    return nBoundaryFaces;
}

NeoFOAM::UnstructuredMesh readOpenFOAMMesh(const NeoFOAM::Executor exec, fvMesh& mesh)
{
    const int32_t nCells = mesh.nCells();
    const int32_t nInternalFaces = mesh.nInternalFaces();
    const int32_t nBoundaryFaces = computeNBoundaryFaces(mesh);
    const int32_t nBoundaries = mesh.boundary().size();
    const int32_t nFaces = mesh.nFaces();

    scalarField magFaceAreas(mag(mesh.faceAreas()));

    labelList faceCells =
        flatBCField<labelList>(mesh, [](const fvPatch& patch) { return patch.faceCells(); });
    vectorField cf =
        flatBCField<vectorField>(mesh, [](const fvPatch& patch) { return patch.Cf(); });
    vectorField cn = flatBCField<vectorField>(
        mesh, [](const fvPatch& patch) { return vectorField(patch.Cn()); }
    );
    vectorField sf =
        flatBCField<vectorField>(mesh, [](const fvPatch& patch) { return patch.Sf(); });
    scalarField magSf =
        flatBCField<scalarField>(mesh, [](const fvPatch& patch) { return patch.magSf(); });
    vectorField nf = flatBCField<vectorField>(
        mesh, [](const fvPatch& patch) { return vectorField(patch.nf()); }
    );
    vectorField delta = flatBCField<vectorField>(
        mesh, [](const fvPatch& patch) { return vectorField(patch.delta()); }
    );
    scalarField weights =
        flatBCField<scalarField>(mesh, [](const fvPatch& patch) { return patch.weights(); });
    scalarField deltaCoeffs =
        flatBCField<scalarField>(mesh, [](const fvPatch& patch) { return patch.deltaCoeffs(); });
    std::vector<NeoFOAM::localIdx> offset = computeOffset(mesh);


    NeoFOAM::BoundaryMesh bMesh(
        exec,
        fromFoamField(exec, faceCells),
        fromFoamField(exec, cf),
        fromFoamField(exec, cn),
        fromFoamField(exec, sf),
        fromFoamField(exec, magSf),
        fromFoamField(exec, nf),
        fromFoamField(exec, delta),
        fromFoamField(exec, weights),
        fromFoamField(exec, deltaCoeffs),
        offset
    );

    NeoFOAM::UnstructuredMesh uMesh(
        fromFoamField(exec, mesh.points()),
        fromFoamField(exec, mesh.cellVolumes()),
        fromFoamField(exec, mesh.cellCentres()),
        fromFoamField(exec, mesh.faceAreas()),
        fromFoamField(exec, mesh.faceCentres()),
        fromFoamField(exec, magFaceAreas),
        fromFoamField(exec, mesh.faceOwner()),
        fromFoamField(exec, mesh.faceNeighbour()),
        nCells,
        nInternalFaces,
        nBoundaryFaces,
        nBoundaries,
        nFaces,
        bMesh
    );

    return uMesh;
}

MeshAdapter::MeshAdapter(const NeoFOAM::Executor exec, const IOobject& io, const bool doInit)
    : fvMesh(io, doInit), nfMesh_(readOpenFOAMMesh(exec, *this))
{
    if (doInit)
    {
        init(false); // do not initialise lower levels
    }
}


MeshAdapter::MeshAdapter(const NeoFOAM::Executor exec, const IOobject& io, const zero, bool syncPar)
    : fvMesh(io, zero {}, syncPar), nfMesh_(readOpenFOAMMesh(exec, *this))
{}


MeshAdapter::MeshAdapter(
    const NeoFOAM::Executor exec,
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
    : fvMesh(
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour),
        syncPar
    ),
      nfMesh_(readOpenFOAMMesh(exec, *this))
{}


MeshAdapter::MeshAdapter(
    const NeoFOAM::Executor exec,
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    cellList&& cells,
    const bool syncPar
)
    : fvMesh(io, std::move(points), std::move(faces), std::move(cells), syncPar),
      nfMesh_(readOpenFOAMMesh(exec, *this))
{}

}
