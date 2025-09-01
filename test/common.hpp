// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 FoamAdapter authors
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

#include "NeoN/NeoN.hpp"
#include "catch2/executorGenerator.hpp"
#include "FoamAdapter/FoamAdapter.hpp"

#include "fvm.H"
#include "fvc.H"
// #include "fvCFD.H"

namespace FoamAdapter
{

template<typename FieldType, typename RandomFunc>
FieldType createRandomField(
    const Foam::Time& runTime,
    const Foam::fvMesh& mesh,
    Foam::word name,
    RandomFunc rand
)
{
    FieldType t(
        Foam::IOobject(
            name,
            runTime.timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::AUTO_WRITE
        ),
        mesh
    );

    for (auto celli = 0; celli < t.size(); celli++)
    {
        t[celli] = rand();
    }

    t.correctBoundaryConditions();
    return t;
}


/* function to create a volScalarField filled with random values for test purposes */
auto randomScalarField(const Foam::Time& runTime, const Foam::fvMesh& mesh, Foam::word name)
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);
    return createRandomField<Foam::volScalarField>(runTime, mesh, name, [&]() { return dis(gen); });
}

auto randomVectorField(
    const Foam::Time& runTime,
    const FoamAdapter::MeshAdapter& mesh,
    Foam::word name
)
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);
    return createRandomField<Foam::volVectorField>(
        runTime,
        mesh,
        name,
        [&]() {
            return Foam::vector {dis(gen), dis(gen), dis(gen)};
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
        auto aBoundaryHost = a.boundaryData().value().copyToHost();
        for (const auto& patch : b.boundaryField())
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


template<typename NFFIELD, typename OFFIELD>
double max_abs_diff(NFFIELD& a, OFFIELD& b){

    auto aHost = a.internalVector().copyToHost();
    auto bSpan = std::span(b.primitiveFieldRef().data(), b.size());

    double maxAbsDiff = 0.0;
    for (std::size_t i = 0; i < aHost.size(); ++i) {
        double absDiff = std::abs(aHost.data()[i] - bSpan[i]);
        if (absDiff > maxAbsDiff) maxAbsDiff = absDiff;
    }
    return maxAbsDiff;
}

template<typename NFFIELD, typename OFFIELD>
double max_rel_diff(NFFIELD& a, OFFIELD& b){

    auto aHost = a.internalVector().copyToHost();
    auto bSpan = std::span(b.primitiveFieldRef().data(), b.size());

    double maxRelDiff = 0.0;
    for (std::size_t i = 0; i < aHost.size(); ++i) {
        double relDiff = std::abs(aHost.data()[i] - bSpan[i]) / std::abs(bSpan[i]);
        if (relDiff > maxRelDiff) maxRelDiff = relDiff;
    }
    return maxRelDiff;
}
} // end namespace FoamAdapter