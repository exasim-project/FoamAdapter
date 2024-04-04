// SPDX-License-Identifier: GPL-3.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "NeoFOAM/cellCentredFiniteVolume/surfaceInterpolation/surfaceInterpolation.hpp"

namespace Foam
{
using namespace NeoFOAM;
class surfaceInterpolationFactory
{
    public:
    //- Runtime type information
    // TypeName("surfaceInterpolationFactory");

    declareRunTimeSelectionTable
    (
        autoPtr,
        surfaceInterpolationKernel,
        dictionary,
        (
            const executor& exec,
            const unstructuredMesh& mesh
        ),
        (exec, mesh)
    );

    //- Return a reference to the selected phaseChange model
    static autoPtr<surfaceInterpolationKernel> New
    (
        const executor& exec,
        const unstructuredMesh& mesh
    );


};

}; // namespace NeoFOAM