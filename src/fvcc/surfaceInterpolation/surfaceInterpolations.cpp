// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "addToRunTimeSelectionTable.H"

#include "FoamAdapter/fvcc/surfaceInterpolation/surfaceInterpolationFactory.hpp"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::addNamedToRunTimeSelectionTable(surfaceInterpolationFactory, linear, dictionary, linear);
