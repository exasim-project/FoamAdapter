// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#pragma once

#include "FoamAdapter/conversion/convert.hpp"

#include "NeoFOAM/fields/field.hpp"
#include "NeoFOAM/fields/boundaryFields.hpp"
#include "NeoFOAM/fields/domainField.hpp"
#include "NeoFOAM/finiteVolume/cellCentred.hpp"

#include "NeoFOAM/finiteVolume/cellCentred/operators/gaussGreenGrad.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/interpolation/linear.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/interpolation/upwind.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/interpolation/surfaceInterpolation.hpp"

#include "NeoFOAM/finiteVolume/cellCentred/stencil/geometryScheme.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/stencil/basicGeometryScheme.hpp"
#include "NeoFOAM/mesh/unstructured.hpp"

#include "fvCFD.H"

#include <catch2/catch_approx.hpp>


struct ApproxScalar
{
    Foam::scalar margin;
    bool operator()(double rhs, double lhs) const
    {
        return Catch::Approx(rhs).margin(margin) == lhs;
    }
};

struct ApproxVector
{
    Foam::scalar margin;
    bool operator()(NeoFOAM::Vector rhs, Foam::vector lhs) const
    {
        NeoFOAM::Vector diff(rhs[0] - lhs[0], rhs[1] - lhs[1], rhs[2] - lhs[2]);

        return Catch::Approx(0).margin(margin) == mag(diff);
    }
};
