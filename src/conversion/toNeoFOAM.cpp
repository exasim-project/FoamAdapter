// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 NeoFOAM authors

#include "FoamAdapter/conversion/toNeoFOAM.hpp"

template<>
NeoFOAM::UnstructuredMesh toNeoFOAM(const NeoFOAM::Executor exec, Foam::fvMesh& mesh)
{
    const int32_t nCells = mesh.nCells();
    const int32_t nInternalFaces = mesh.nInternalFaces();
    const int32_t nBoundaryFaces = computeNBoundaryFaces(mesh);
    const int32_t nBoundaries = mesh.boundary().size();
    const int32_t nFaces = mesh.nFaces();

    Foam::scalarField magFaceAreas(mag(mesh.faceAreas()));

    Foam::labelList faceCells = flatBCField<Foam::labelList>(
        mesh, [](const Foam::fvPatch& patch) { return patch.faceCells(); }
    );
    Foam::vectorField cf =
        flatBCField<Foam::vectorField>(mesh, [](const Foam::fvPatch& patch) { return patch.Cf(); });
    Foam::vectorField cn = flatBCField<Foam::vectorField>(
        mesh, [](const Foam::fvPatch& patch) { return Foam::vectorField(patch.Cn()); }
    );
    Foam::vectorField sf =
        flatBCField<Foam::vectorField>(mesh, [](const Foam::fvPatch& patch) { return patch.Sf(); });
    Foam::scalarField magSf = flatBCField<Foam::scalarField>(
        mesh, [](const Foam::fvPatch& patch) { return patch.magSf(); }
    );
    Foam::vectorField nf = flatBCField<Foam::vectorField>(
        mesh, [](const Foam::fvPatch& patch) { return Foam::vectorField(patch.nf()); }
    );
    Foam::vectorField delta = flatBCField<Foam::vectorField>(
        mesh, [](const Foam::fvPatch& patch) { return Foam::vectorField(patch.delta()); }
    );
    Foam::scalarField weights = flatBCField<Foam::scalarField>(
        mesh, [](const Foam::fvPatch& patch) { return patch.weights(); }
    );
    Foam::scalarField deltaCoeffs = flatBCField<Foam::scalarField>(
        mesh, [](const Foam::fvPatch& patch) { return patch.deltaCoeffs(); }
    );
    std::vector<NeoFOAM::localIdx> offset = computeOffset(mesh);

    NeoFOAM::BoundaryMesh bMesh(
        exec,
        toNeoFOAM(exec, faceCells),
        toNeoFOAM(exec, cf),
        toNeoFOAM(exec, cn),
        toNeoFOAM(exec, sf),
        toNeoFOAM(exec, magSf),
        toNeoFOAM(exec, nf),
        toNeoFOAM(exec, delta),
        toNeoFOAM(exec, weights),
        toNeoFOAM(exec, deltaCoeffs),
        offset
    );

    NeoFOAM::UnstructuredMesh uMesh(
        toNeoFOAM(exec, mesh.points()),
        toNeoFOAM(exec, mesh.cellVolumes()),
        toNeoFOAM(exec, mesh.cellCentres()),
        toNeoFOAM(exec, mesh.faceAreas()),
        toNeoFOAM(exec, mesh.faceCentres()),
        toNeoFOAM(exec, magFaceAreas),
        toNeoFOAM(exec, mesh.faceOwner()),
        toNeoFOAM(exec, mesh.faceNeighbour()),
        nCells,
        nInternalFaces,
        nBoundaryFaces,
        nBoundaries,
        nFaces,
        bMesh
    );

    return uMesh;
}
