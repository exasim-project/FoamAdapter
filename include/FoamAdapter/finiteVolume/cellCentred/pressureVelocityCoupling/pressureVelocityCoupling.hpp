// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#pragma once


#include "NeoFOAM/finiteVolume/cellCentred/fields/volumeField.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/fields/surfaceField.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/dsl/expression.hpp"


namespace NeoFOAM::finiteVolume::cellCentred
{

void constrainHbyA(
    VolumeField<Vector>& HbyA,
    const VolumeField<Vector>& U,
    const VolumeField<scalar>& p
);

std::tuple<VolumeField<scalar>, VolumeField<Vector>>
discreteMomentumFields(const Expression<Vector>& expr);

void updateFaceVelocity(
    SurfaceField<scalar>& phi,
    const SurfaceField<scalar>& predictedPhi,
    const Expression<scalar>& expr
);

void updateVelocity(
    VolumeField<Vector>& U,
    const VolumeField<Vector>& HbyA,
    VolumeField<scalar>& rAU,
    VolumeField<scalar>& p
);


SurfaceField<scalar> flux(const VolumeField<Vector>& volField);

}
