// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "FoamAdapter/mesh/meshAdapter.hpp"
#include "FoamAdapter/conversion/mesh.hpp"
#include "FoamAdapter/conversion/toNeoFOAM.hpp"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(MeshAdapter, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MeshAdapter::MeshAdapter(const NeoFOAM::Executor exec, const IOobject& io, const bool doInit)
    : fvMesh(io, doInit), mesh_(toNeoFOAM(exec, *this))
{
    if (doInit)
    {
        init(false); // do not initialise lower levels
    }
}


Foam::MeshAdapter::MeshAdapter(
    const NeoFOAM::Executor exec, const IOobject& io, const Foam::zero, bool syncPar
)
    : fvMesh(io, Foam::zero {}, syncPar), mesh_(toNeoFOAM(exec, *this))
{}


Foam::MeshAdapter::MeshAdapter(
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
      mesh_(toNeoFOAM(exec, *this))
{}


Foam::MeshAdapter::MeshAdapter(
    const NeoFOAM::Executor exec,
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    cellList&& cells,
    const bool syncPar
)
    : fvMesh(io, std::move(points), std::move(faces), std::move(cells), syncPar),
      mesh_(toNeoFOAM(exec, *this))
{}


// ************************************************************************* //
