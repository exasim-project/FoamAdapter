// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/mesh/unstructuredMesh/BoundaryMesh.hpp"
#include "NeoFOAM/mesh/unstructuredMesh/unstructuredMesh.hpp"
#include "NeoFOAM/primitives/label.hpp"
#include "foamFields.hpp"
#include "fvMesh.H"
#include <functional>

namespace Foam
{

std::vector<NeoFOAM::localIdx> computeOffset(const Foam::fvMesh& mesh);

int32_t computeNBoundaryFaces(const Foam::fvMesh& mesh);

template<typename FieldT>
FieldT flatBCField(const Foam::fvMesh& mesh, std::function<FieldT(const Foam::fvPatch&)> f)
{
    FieldT result(computeNBoundaryFaces(mesh));
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


NeoFOAM::unstructuredMesh readOpenFOAMMesh(const NeoFOAM::executor exec, Foam::fvMesh& mesh);

}; // namespace Foam
