// SPDX-License-Identifier: MPL-2.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors


#include "addToRunTimeSelectionTable.H"

#include "FoamAdapter/fvcc/surfaceInterpolation/surfaceInterpolationFactory.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/surfaceInterpolation/linear.hpp"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
using namespace NeoFOAM;

addNamedToRunTimeSelectionTable
(
    surfaceInterpolationFactory,
    linear,
    dictionary,
    linear
);

}