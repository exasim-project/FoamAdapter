// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "FoamAdapter/fvcc/mesh/MeshAdapter.hpp"
#include "FoamAdapter/readers/foamMesh.hpp"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(MeshAdapter, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MeshAdapter::fvccNeoMesh(const NeoFOAM::Executor exec, const IOobject& io, const bool doInit)
    : fvMesh(io, doInit), uMesh_(readOpenFOAMMesh(exec, *this))
{
    if (doInit)
    {
        init(false); // do not initialise lower levels
    }
}


Foam::MeshAdapter::fvccNeoMesh(
    const NeoFOAM::Executor exec, const IOobject& io, const Foam::zero, bool syncPar
)
    : fvMesh(io, Foam::zero {}, syncPar), uMesh_(readOpenFOAMMesh(exec, *this))
{}


Foam::MeshAdapter::fvccNeoMesh(
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
      uMesh_(readOpenFOAMMesh(exec, *this))
{}


Foam::MeshAdapter::fvccNeoMesh(
    const NeoFOAM::Executor exec,
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    cellList&& cells,
    const bool syncPar
)
    : fvMesh(io, std::move(points), std::move(faces), std::move(cells), syncPar),
      uMesh_(readOpenFOAMMesh(exec, *this))
{}


// ************************************************************************* //
