// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
// common test functions
#pragma once

#include <random>
#include <span>

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/catch_approx.hpp>
#include "catch2/common.hpp"

#include "FoamAdapter/meshAdapter.hpp"
#include "FoamAdapter/setup.hpp"
#include "FoamAdapter/comparison.hpp"

#include "fvm.H"
#include "fvc.H"
// #include "fvCFD.H"

namespace Foam
{

template<typename FieldType, typename RandomFunc>
FieldType createRandomField(const Time& runTime, const fvMesh& mesh, word name, RandomFunc rand)
{
    FieldType t(
        IOobject(name, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
        mesh
    );

    forAll(t, celli)
    {
        t[celli] = rand();
    }

    t.correctBoundaryConditions();
    return t;
}


/* function to create a volScalarField filled with random values for test purposes */
auto randomScalarField(const Time& runTime, const fvMesh& mesh, word name)
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);
    return createRandomField<volScalarField>(runTime, mesh, name, [&]() { return dis(gen); });
}

auto randomVectorField(const Time& runTime, const MeshAdapter& mesh, word name)
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);
    return createRandomField<volVectorField>(
        runTime,
        mesh,
        name,
        [&]() {
            return vector {dis(gen), dis(gen), dis(gen)};
        }
    );
}

/* comparison function for volumeFields */
template<typename NFFIELD, typename OFFIELD, typename Compare>
void compare(NFFIELD& a, OFFIELD& b, Compare comp, bool withBoundaries = true)
{
    auto aHost = a.internalVector().copyToHost();
    auto bSpan = std::span(b.primitiveFieldRef().data(), b.size());
    // nf a span might be shorter than bSpan for surface fields
    REQUIRE_THAT(aHost.view({0, bSpan.size()}), Catch::Matchers::RangeEquals(bSpan, comp));

    if (withBoundaries)
    {
        size_t start = 0;
        auto aBoundaryHost = a.boundaryVector().value().copyToHost();
        for (const auto& patch : b.boundaryVector())
        {
            auto bBoundarySpan = std::span(patch.cdata(), patch.size());
            REQUIRE_THAT(
                aBoundaryHost.view({start, start + patch.size()}),
                Catch::Matchers::RangeEquals(bBoundarySpan, comp)
            );
            start += patch.size();
        }
    }
}

}
