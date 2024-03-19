// SPDX-License-Identifier: GPLv-3.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "fvMesh.H"
#include "NeoFOAM/core/executor/executor.hpp"
#include "foamFields.hpp"
#include "NeoFOAM/mesh/unstructuredMesh/unstructuredMesh.hpp"
#include "NeoFOAM/mesh/unstructuredMesh/BoundaryMesh.hpp"
#include "NeoFOAM/primitives/label.hpp"
#include <functional>

std::vector<NeoFOAM::localIdx> computeOffset(const Foam::fvMesh &mesh)
{
    std::vector<NeoFOAM::localIdx> result;
    const Foam::fvBoundaryMesh& bMesh = mesh.boundary();
    result.push_back(0);
    forAll(bMesh, patchI)
    {
        const Foam::fvPatch& patch = bMesh[patchI];
        result.push_back(patch.start());
    }
    return result;
}

template <typename FieldT>
FieldT flatBCField(const Foam::fvMesh &mesh, std::function<FieldT(const Foam::fvPatch&)> f)
{
    FieldT result(mesh.nFaces() - mesh.nInternalFaces());
    const Foam::fvBoundaryMesh& bMesh = mesh.boundary();
    Foam::label idx = 0;
    forAll(bMesh, patchI)
    {
        const Foam::fvPatch& patch = bMesh[patchI];
        auto pResult = f(patch);
        forAll(pResult, i)
        {
            result[idx] = pResult[i];
            idx++;
        }
    }
    return result;
}

NeoFOAM::unstructuredMesh readOpenFOAMMesh(const NeoFOAM::executor exec, Foam::fvMesh &mesh)
{
    const int32_t nCells = mesh.nCells();
    const int32_t nInternalFaces = mesh.nInternalFaces();

    Foam::scalarField magFaceAreas = mag(mesh.faceAreas());

    Foam::labelList faceCells = flatBCField<Foam::labelList>(mesh, [](const Foam::fvPatch& patch){ return patch.faceCells(); });
    Foam::vectorField Cf = flatBCField<Foam::vectorField>(mesh, [](const Foam::fvPatch& patch){ return patch.Cf(); });
    Foam::vectorField Cn = flatBCField<Foam::vectorField>(mesh, [](const Foam::fvPatch& patch){ return Foam::vectorField(patch.Cn()); });
    Foam::vectorField Sf = flatBCField<Foam::vectorField>(mesh, [](const Foam::fvPatch& patch){ return patch.Sf(); });
    Foam::scalarField magSf = flatBCField<Foam::scalarField>(mesh, [](const Foam::fvPatch& patch){ return patch.magSf(); });
    Foam::vectorField nf = flatBCField<Foam::vectorField>(mesh, [](const Foam::fvPatch& patch){ return Foam::vectorField (patch.nf()); });
    Foam::vectorField delta = flatBCField<Foam::vectorField>(mesh, [](const Foam::fvPatch& patch){ return Foam::vectorField (patch.delta()); });
    Foam::scalarField weights = flatBCField<Foam::scalarField>(mesh, [](const Foam::fvPatch& patch){ return patch.weights(); });
    Foam::scalarField deltaCoeffs = flatBCField<Foam::scalarField>(mesh, [](const Foam::fvPatch& patch){ return patch.deltaCoeffs(); });
    std::vector<NeoFOAM::localIdx> offset = computeOffset(mesh);

    NeoFOAM::BoundaryMesh bMesh(
        exec,
        fromFoamField<NeoFOAM::label , Foam::label> (exec, faceCells),
        fromFoamField<NeoFOAM::Vector, Foam::vector>(exec, Cf),
        fromFoamField<NeoFOAM::Vector, Foam::vector>(exec, Cn),
        fromFoamField<NeoFOAM::Vector, Foam::vector>(exec, Sf),
        fromFoamField<NeoFOAM::scalar, Foam::scalar>(exec, magSf),
        fromFoamField<NeoFOAM::Vector, Foam::vector>(exec, nf),
        fromFoamField<NeoFOAM::Vector, Foam::vector>(exec, delta),
        fromFoamField<NeoFOAM::scalar, Foam::scalar>(exec, weights),
        fromFoamField<NeoFOAM::scalar, Foam::scalar>(exec, deltaCoeffs),
        offset
    );
       
    NeoFOAM::unstructuredMesh uMesh(
        fromFoamField<NeoFOAM::Vector,Foam::vector>(exec, mesh.points()),
        fromFoamField<NeoFOAM::scalar,Foam::scalar>(exec, mesh.cellVolumes()),
        fromFoamField<NeoFOAM::Vector,Foam::vector>(exec, mesh.cellCentres()),
        fromFoamField<NeoFOAM::Vector,Foam::vector>(exec, mesh.faceAreas() ),
        fromFoamField<NeoFOAM::Vector,Foam::vector>(exec, mesh.faceCentres()),
        fromFoamField<NeoFOAM::scalar,Foam::scalar>(exec, magFaceAreas),
        fromFoamField<NeoFOAM::label,Foam::label>(exec, mesh.faceOwner()),
        fromFoamField<NeoFOAM::label,Foam::label>(exec, mesh.faceNeighbour()),
        nCells,
        nInternalFaces,
        bMesh
    );

    return uMesh;
}