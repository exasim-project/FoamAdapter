// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "catch2/common.hpp"

#include "NeoFOAM/fields/field.hpp"

#include "FoamAdapter/readers/foamDictionary.hpp"
#include "FoamAdapter/writers.hpp"

extern Foam::Time* timePtr; // A single time object

TEST_CASE("read dict")
{
    Foam::Time& runTime = *timePtr;
    Foam::Info << "\nReading testDict" << Foam::endl;

    Foam::dictionary testDict;
    testDict.add("label", 1);
    testDict.add("scalar", 2.1);
    testDict.add("scalar2", 2.0);
    testDict.add("vector", Foam::vector(1.0, 2.0, 3.0));
    testDict.add("word", "word");
    // testDict.add("divScheme", "Gauss linear"); // will throw error in NeoFoam as conversion is
    // not implemented

    Foam::dictionary subDict;
    subDict.add("subScalar", 4.1);
    subDict.add("subVector", Foam::vector(5.0, 6.0, 7.0));
    subDict.add("subWord", "subWord");
    testDict.add("subDict", subDict);

    NeoFOAM::Dictionary neoDict = Foam::readFoamDictionary(testDict);

    REQUIRE(neoDict.get<NeoFOAM::label>("label") == 1);
    REQUIRE(neoDict.get<NeoFOAM::scalar>("scalar") == 2.1);
    REQUIRE(neoDict.get<NeoFOAM::scalar>("scalar2") == 2.0);
    REQUIRE(neoDict.get<NeoFOAM::Vector>("vector") == NeoFOAM::Vector(1.0, 2.0, 3.0));
    REQUIRE(neoDict.get<std::string>("word") == "word");
    // REQUIRE(neoDict.get<std::string>("divScheme") == "2.0"); // not added to the dictionary


    NeoFOAM::Dictionary& subNeoDict = neoDict.subDict("subDict");
    REQUIRE(subNeoDict.get<NeoFOAM::scalar>("subScalar") == 4.1);
    REQUIRE(subNeoDict.get<NeoFOAM::Vector>("subVector") == NeoFOAM::Vector(5.0, 6.0, 7.0));
}
