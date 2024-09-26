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

template<typename ValueType>
void checkField(const NeoFOAM::Field<ValueType>& field, ValueType value)
{
    auto fieldHost = field.copyToHost();
    auto fieldSpan = fieldHost.span();

    for (int i = 0; i < fieldSpan.size(); i++)
    {
        REQUIRE(fieldSpan[i] == value);
    }
}

template<typename NF_FIELD, typename OF_FIELD>
void checkField(const NF_FIELD& a, const OF_FIELD& b)
{
    auto aHost = a.copyToHost();
    auto aSpan = aHost.span();

    REQUIRE(a.size() == b.size());
    for (int i = 0; i < aSpan.size(); i++)
    {
        REQUIRE(aSpan[i] == convert(b[i]));
    }
}


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
