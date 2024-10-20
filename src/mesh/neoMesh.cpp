// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "FoamAdapter/mesh/neoMesh.hpp"
#include "FoamAdapter/readers/foamMesh.hpp"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(fvccNeoMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvccNeoMesh::fvccNeoMesh(const NeoFOAM::Executor exec, const IOobject& io, const bool doInit)
    : fvMesh(io, doInit), nfMesh_(readOpenFOAMMesh(exec, *this))
{
    if (doInit)
    {
        init(false); // do not initialise lower levels
    }
}


Foam::fvccNeoMesh::fvccNeoMesh(
    const NeoFOAM::Executor exec, const IOobject& io, const Foam::zero, bool syncPar
)
    : fvMesh(io, Foam::zero {}, syncPar), nfMesh_(readOpenFOAMMesh(exec, *this))
{}


Foam::fvccNeoMesh::fvccNeoMesh(
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


Foam::fvccNeoMesh::fvccNeoMesh(
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


// ************************************************************************* //
