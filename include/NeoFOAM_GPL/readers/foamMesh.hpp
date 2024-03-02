// SPDX-License-Identifier: GPLv-3.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "fvMesh.H"
#include "NeoFOAM//core/executor/executor.hpp"
#include "foamFields.hpp"
#include "NeoFOAM/mesh/unstructuredMesh/unstructuredMesh.hpp"
#include "NeoFOAM/primitives/label.hpp"

NeoFOAM::unstructuredMesh readOpenFOAMMesh(const NeoFOAM::executor exec, Foam::fvMesh &mesh)
{
    const int32_t nCells = mesh.nCells();
    const int32_t nInternalFaces = mesh.nInternalFaces();

    Foam::scalarField magFaceAreas = mag(mesh.faceAreas());

    NeoFOAM::unstructuredMesh uMesh(
        fromFoamField<NeoFOAM::vector,Foam::vector>(exec, mesh.points()),
        fromFoamField<NeoFOAM::scalar,Foam::scalar>(exec, mesh.cellVolumes()),
        fromFoamField<NeoFOAM::vector,Foam::vector>(exec, mesh.cellCentres()),
        fromFoamField<NeoFOAM::vector,Foam::vector>(exec, mesh.faceAreas() ),
        fromFoamField<NeoFOAM::vector,Foam::vector>(exec, mesh.faceCentres()),
        fromFoamField<NeoFOAM::scalar,Foam::scalar>(exec, magFaceAreas),
        fromFoamField<NeoFOAM::label,Foam::label>(exec, mesh.faceOwner()),
        fromFoamField<NeoFOAM::label,Foam::label>(exec, mesh.faceNeighbour()),
        nCells,
        nInternalFaces
    );

    return uMesh;
}