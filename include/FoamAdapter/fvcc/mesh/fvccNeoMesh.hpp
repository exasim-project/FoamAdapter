// SPDX-License-Identifier: GPL-3.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "fvMesh.H"
#include "NeoFOAM/mesh/unstructuredMesh/unstructuredMesh.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fvccNeoMesh Declaration
\*---------------------------------------------------------------------------*/

class fvccNeoMesh
:
    public Foam::fvMesh
{
    // Private Data
        const NeoFOAM::executor exec;

        NeoFOAM::unstructuredMesh uMesh_;


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
        explicit fvccNeoMesh(const NeoFOAM::executor exec,const IOobject& io, const bool doInit=true);

        //- Construct from IOobject or as zero-sized mesh
        //  Boundary is added using addFvPatches() member function
        fvccNeoMesh(const NeoFOAM::executor exec,const IOobject& io, const Foam::zero, bool syncPar=true);

        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        fvccNeoMesh
        (
            const NeoFOAM::executor exec,
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            labelList&& allOwner,
            labelList&& allNeighbour,
            const bool syncPar = true
        );

        //- Construct without boundary from cells rather than owner/neighbour.
        //  Boundary is added using addPatches() member function
        fvccNeoMesh
        (
            const NeoFOAM::executor exec,
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            cellList&& cells,
            const bool syncPar = true
        );


    //- Destructor
    virtual ~fvccNeoMesh() = default;

    NeoFOAM::unstructuredMesh& uMesh()
    {
        return uMesh_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}; // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
