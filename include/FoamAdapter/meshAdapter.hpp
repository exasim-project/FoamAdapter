// SPDX-License-Identifier: GPL-3.0-or-later
//
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#pragma once

#include <functional>

#include "fvMesh.H"

#include "NeoFOAM/mesh/unstructured.hpp"
#include "NeoFOAM/core/primitives/label.hpp"
#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/mesh/unstructured.hpp"

#include "readers/foamFields.hpp"

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

NeoFOAM::UnstructuredMesh readOpenFOAMMesh(const NeoFOAM::Executor exec, Foam::fvMesh& mesh);

/** @class fvccNeoMesh
 */
class fvccNeoMesh : public Foam::fvMesh
{
    // Private Data
    const NeoFOAM::Executor exec;

    NeoFOAM::UnstructuredMesh nfMesh_;

    // Private Member Functions

    //- No copy construct
    fvccNeoMesh(const fvccNeoMesh&) = delete;

    //- No copy assignment
    void operator=(const fvccNeoMesh&) = delete;

public:

    //- Runtime type information
    TypeName("fvccNeoMesh");

    // Constructors

    //- Construct from IOobject
    explicit fvccNeoMesh(
        const NeoFOAM::Executor exec, const IOobject& io, const bool doInit = true
    );

    //- Construct from IOobject or as zero-sized mesh
    //  Boundary is added using addFvPatches() member function
    fvccNeoMesh(
        const NeoFOAM::Executor exec, const IOobject& io, const Foam::zero, bool syncPar = true
    );

    //- Construct from components without boundary.
    //  Boundary is added using addFvPatches() member function
    fvccNeoMesh(
        const NeoFOAM::Executor exec,
        const IOobject& io,

        pointField&& points,
        faceList&& faces,
        labelList&& allOwner,
        labelList&& allNeighbour,
        const bool syncPar = true
    );

    //- Construct without boundary from cells rather than owner/neighbour.
    //  Boundary is added using addPatches() member function
    fvccNeoMesh(
        const NeoFOAM::Executor exec,
        const IOobject& io,
        pointField&& points,
        faceList&& faces,
        cellList&& cells,
        const bool syncPar = true
    );

    //- Destructor
    virtual ~fvccNeoMesh() = default;

    NeoFOAM::UnstructuredMesh& nfMesh() { return nfMesh_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}; // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
