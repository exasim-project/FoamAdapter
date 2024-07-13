// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "NeoFOAM/finiteVolume/interpolation/surfaceInterpolation.hpp"

namespace Foam
{
using namespace NeoFOAM;
class surfaceInterpolationFactory
{
public:

    //- Runtime type information
    // TypeName("surfaceInterpolationFactory");

    declareRunTimeSelectionTable(
        autoPtr,
        SurfaceInterpolationKernel,
        dictionary,
        (const Executor& exec, const UnstructuredMesh& mesh),
        (exec, mesh)
    );

    //- Return a reference to the selected phaseChange model
    static autoPtr<SurfaceInterpolationKernel>
    New(const Executor& exec, const UnstructuredMesh& mesh);
};

}; // namespace NeoFOAM
