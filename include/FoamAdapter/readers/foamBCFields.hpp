// SPDX-License-Identifier: GPLv-3.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/core/Dictionary.hpp"
#include "FoamAdapter/conversion/convert.hpp"
#include "NeoFOAM/fields/FieldTypeDefs.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/fields/fvccVolField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/fvccBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/vol/scalar/fvccScalarFixedValueBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/vol/scalar/fvccScalarZeroGradientBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/vol/scalar/fvccScalarEmptyBoundaryField.hpp"

#include "FoamAdapter/conversion/type_conversion.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/fvccBoundaryFieldSelector.hpp"

template <typename Type>
auto readBoundaryCondition(
    const NeoFOAM::unstructuredMesh& uMesh,
    int patchi,
    const NeoFOAM::Dictionary& patchDict)
{
    return NeoFOAM::getBC<Type>(uMesh, patchi, patchDict);
}

