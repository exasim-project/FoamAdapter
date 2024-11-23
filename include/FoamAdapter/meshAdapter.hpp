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

#include "readers.hpp"

namespace Foam
{

std::vector<NeoFOAM::localIdx> computeOffset(const fvMesh& mesh);

int32_t computeNBoundaryFaces(const fvMesh& mesh);

template<typename FieldT>
FieldT flatBCField(const fvMesh& mesh, std::function<FieldT(const fvPatch&)> f);

NeoFOAM::UnstructuredMesh readOpenFOAMMesh(const NeoFOAM::Executor exec, const fvMesh& mesh);

/** @class MeshAdapter
 */
class MeshAdapter : public fvMesh
{

    NeoFOAM::UnstructuredMesh nfMesh_;

    // Private Member Functions

    //- No copy construct
    MeshAdapter(const MeshAdapter&) = delete;

    //- No copy assignment
    void operator=(const MeshAdapter&) = delete;

public:

    //- Runtime type information
    TypeName("MeshAdapter");

    // Constructors

    //- Construct from IOobject
    explicit MeshAdapter(
        const NeoFOAM::Executor exec,
        const IOobject& io,
        const bool doInit = true
    );

    //- Construct from IOobject or as zero-sized mesh
    //  Boundary is added using addFvPatches() member function
    MeshAdapter(const NeoFOAM::Executor exec, const IOobject& io, const zero, bool syncPar = true);

    //- Construct from components without boundary.
    //  Boundary is added using addFvPatches() member function
    MeshAdapter(
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
    MeshAdapter(
        const NeoFOAM::Executor exec,
        const IOobject& io,
        pointField&& points,
        faceList&& faces,
        cellList&& cells,
        const bool syncPar = true
    );

    //- Destructor
    virtual ~MeshAdapter() = default;

    NeoFOAM::UnstructuredMesh& nfMesh() { return nfMesh_; }

    const NeoFOAM::UnstructuredMesh& nfMesh() const { return nfMesh_; }

    const NeoFOAM::Executor exec() const { return nfMesh().exec(); }
};

} // End namespace Foam
