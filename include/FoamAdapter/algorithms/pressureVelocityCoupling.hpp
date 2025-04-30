// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 FoamAdapter authors

#pragma once

#include "FoamAdapter/datastructures/expression.hpp"

#include "NeoN/NeoN.hpp"

namespace nnfvcc = NeoN::finiteVolume::cellCentred;
using scalar = NeoN::scalar;
using Vec3 = NeoN::Vec3;

namespace FoamAdapter
{

/* @brief
 *
 */
void constrainHbyA(
    nnfvcc::VolumeField<Vec3>& HbyA,
    const nnfvcc::VolumeField<Vec3>& U,
    const nnfvcc::VolumeField<scalar>& p
);

std::tuple<nnfvcc::VolumeField<scalar>, nnfvcc::VolumeField<Vec3>>
discreteMomentumFields(const Expression<Vec3>& expr);

void updateFaceVelocity(
    nnfvcc::SurfaceField<scalar>& phi,
    const nnfvcc::SurfaceField<scalar>& predictedPhi,
    const Expression<scalar>& expr
);

void updateVelocity(
    nnfvcc::VolumeField<Vec3>& U,
    const nnfvcc::VolumeField<Vec3>& HbyA,
    nnfvcc::VolumeField<scalar>& rAU,
    nnfvcc::VolumeField<scalar>& p
);


/* @brief Reimplementation of OpenFOAMs fvMatrix.flux()
 * @return flux surface field
 */
nnfvcc::SurfaceField<scalar> flux(const nnfvcc::VolumeField<Vec3>& volField);

}
