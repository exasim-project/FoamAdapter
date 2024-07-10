// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "FoamAdapter/conversion/convert.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/fvccBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/vol/scalar/fvccScalarEmptyBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/vol/scalar/fvccScalarFixedValueBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/vol/scalar/fvccScalarZeroGradientBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/fields/fvccVolField.hpp"
#include "NeoFOAM/core/Dictionary.hpp"
#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/fields/FieldTypeDefs.hpp"

#include "FoamAdapter/conversion/type_conversion.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/fvccBoundaryFieldSelector.hpp"

namespace Foam
{

template<typename Type>
auto readVolBoundaryCondition(
    const NeoFOAM::unstructuredMesh& uMesh, int patchi, const NeoFOAM::Dictionary& patchDict
)
{
    return NeoFOAM::getBC<Type>(uMesh, patchi, patchDict);
}

template<typename Type>
auto readSurfaceBoundaryCondition(
    const NeoFOAM::unstructuredMesh& uMesh, int patchi, const NeoFOAM::Dictionary& patchDict
)
{
    return NeoFOAM::getSurfaceBC<Type>(uMesh, patchi, patchDict);
}

}; // namespace Foam
