// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#pragma once

#include <functional>

#include "fvMesh.H"

#include "NeoFOAM/core/primitives/label.hpp"
#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/mesh/unstructured.hpp"

#include "FoamAdapter/conversion/fields.hpp"

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


}; // namespace Foam
