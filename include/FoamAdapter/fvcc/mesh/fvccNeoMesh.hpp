// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "NeoFOAM/mesh/unstructured.hpp"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fvccNeoMesh Declaration
\*---------------------------------------------------------------------------*/

class fvccNeoMesh : public Foam::fvMesh
{
    // Private Data
    const NeoFOAM::Executor exec;

    NeoFOAM::UnstructuredMesh uMesh_;

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

    NeoFOAM::UnstructuredMesh& uMesh() { return uMesh_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}; // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
