// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "fvMesh.H"

#include "NeoFOAM/mesh/unstructured.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/** @class MeshAdapter adapter class simplifying conversion between an unstructured
 * OpenFOAM mesh and the corresponding NeoFOAM datastructure.
 */
class MeshAdapter : public Foam::fvMesh
{
    // Private Data
    const NeoFOAM::Executor exec;

    NeoFOAM::UnstructuredMesh mesh_;

    //- No copy construct
    MeshAdapter(const MeshAdapter&) = delete;

    //- No copy assignment
    void operator=(const MeshAdapter&) = delete;

public:

    //- Runtime type information
    TypeName("MeshAdapter");

    // Constructors
    /* @brief Construct from IOobject
     * */
    explicit MeshAdapter(const NeoFOAM::Executor exec, const fvMesh& mesh);

    /* @brief Construct from IOobject
     * */
    explicit MeshAdapter(
        const NeoFOAM::Executor exec, const IOobject& io, const bool doInit = true
    );

    /* @brief Construct from IOobject or as zero-sized mesh
     * Boundary is added using addFvPatches() member function
     */
    MeshAdapter(
        const NeoFOAM::Executor exec, const IOobject& io, const Foam::zero, bool syncPar = true
    );

    /* @brief Construct from components without boundary.
     * Boundary is added using addFvPatches() member function
     * */
    MeshAdapter(
        const NeoFOAM::Executor exec,
        const IOobject& io,
        pointField&& points,
        faceList&& faces,
        labelList&& allOwner,
        labelList&& allNeighbour,
        const bool syncPar = true
    );

    /* @brief  Construct without boundary from cells rather than owner/neighbour.
     *
     * Boundary is added using addPatches() member function
     */
    MeshAdapter(
        const NeoFOAM::Executor exec,
        const IOobject& io,
        pointField&& points,
        faceList&& faces,
        cellList&& cells,
        const bool syncPar = true
    );

    /* @brief Destructor
     * */
    virtual ~MeshAdapter() = default;

    NeoFOAM::UnstructuredMesh& mesh() { return mesh_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}; // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
