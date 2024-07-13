// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "NeoFOAM/finiteVolume/cellCentred.hpp"
#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/fields/field.hpp"

#include "FoamAdapter/conversion/convert.hpp"
#include "FoamAdapter/conversion/type_conversion.hpp"

namespace Foam
{

template<typename Type>
auto readVolBoundaryCondition(
    const NeoFOAM::UnstructuredMesh& uMesh, int patchi, const NeoFOAM::Dictionary& patchDict
)
{
    // return NeoFOAM::getBC<Type>(uMesh, patchi, patchDict);
}

template<typename Type>
auto readSurfaceBoundaryCondition(
    const NeoFOAM::UnstructuredMesh& uMesh, int patchi, const NeoFOAM::Dictionary& patchDict
)
{
    // return NeoFOAM::getSurfaceBC<Type>(uMesh, patchi, patchDict);
}

}; // namespace Foam
