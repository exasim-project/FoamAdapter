// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 FoamAdapter authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "common.hpp"

extern Foam::Time* timePtr; // A single time object

TEST_CASE("read dict")
{
    Foam::Info << "\nReading testDict" << Foam::endl;

    Foam::dictionary testDict;
    testDict.add("label", 1);
    testDict.add("scalar", 2.1);
    testDict.add("scalar2", 2.0);
    testDict.add("vector", Foam::vector(1.0, 2.0, 3.0));
    testDict.add("word", "word");

    Foam::dictionary subDict;
    subDict.add("subScalar", 4.1);
    subDict.add("subVector", Foam::vector(5.0, 6.0, 7.0));
    subDict.add("subWord", "subWord");
    testDict.add("subDict", subDict);

    NeoN::Dictionary neoDict = FoamAdapter::convert(testDict);

    REQUIRE(neoDict.get<NeoN::label>("label") == 1);
    REQUIRE(neoDict.get<NeoN::scalar>("scalar") == 2.1);
    REQUIRE(neoDict.get<NeoN::scalar>("scalar2") == 2.0);
    REQUIRE(neoDict.get<NeoN::Vec3>("vector") == NeoN::Vec3(1.0, 2.0, 3.0));
    REQUIRE(neoDict.get<std::string>("word") == "word");


    NeoN::Dictionary& subNeoDict = neoDict.subDict("subDict");
    REQUIRE(subNeoDict.get<NeoN::scalar>("subScalar") == 4.1);
    REQUIRE(subNeoDict.get<NeoN::Vec3>("subVector") == NeoN::Vec3(5.0, 6.0, 7.0));
}


TEST_CASE("read fvSchemes")
{
    Foam::Time& runTime = *timePtr;
    auto meshPtr = FoamAdapter::createMesh(NeoN::SerialExecutor {}, runTime);
    FoamAdapter::MeshAdapter& mesh = *meshPtr;

    Foam::dictionary fvSchemes = mesh.schemesDict();
    Foam::Info << "reading fvSchemes" << fvSchemes << Foam::endl;

    NeoN::Dictionary fvSchemesDict = FoamAdapter::convert(mesh.schemesDict());
    auto keys = fvSchemesDict.keys();

    REQUIRE(fvSchemesDict.subDict("ddtSchemes").get<std::string>("ddt(T)") == "Euler");
    auto gradSchemeKeys = fvSchemesDict.subDict("gradSchemes").keys();
    NeoN::TokenList limitedToken =
        fvSchemesDict.subDict("gradSchemes").get<NeoN::TokenList>("limited");
    REQUIRE(limitedToken.size() == 4);
    REQUIRE(limitedToken.get<std::string>(0) == "cellLimited");
    REQUIRE(limitedToken.get<std::string>(1) == "Gauss");
    REQUIRE(limitedToken.get<std::string>(2) == "linear");
    REQUIRE(limitedToken.get<NeoN::label>(3) == 1);

    NeoN::TokenList gradU = fvSchemesDict.subDict("gradSchemes").get<NeoN::TokenList>("grad(U)");
    REQUIRE(gradU.size() == 2);
    REQUIRE(gradU.get<std::string>(0) == "Gauss");
    REQUIRE(gradU.get<std::string>(1) == "linear");
}

TEST_CASE("read testDictionary")
{
    Foam::Time& runTime = *timePtr;
    Foam::IOdictionary ofTestDict(Foam::IOobject(
        "testDictionary",
        runTime.system(),
        runTime,
        Foam::IOobject::MUST_READ,
        Foam::IOobject::NO_WRITE
    ));

    NeoN::Dictionary nfTestDict = Foam::readFoamDictionary(ofTestDict);

    REQUIRE(nfTestDict.get<NeoN::label>("label") == 1);
    REQUIRE(nfTestDict.get<NeoN::scalar>("scalar") == 2.1);
    REQUIRE(nfTestDict.get<NeoN::scalar>("scalar2") == 2.0);
    REQUIRE(nfTestDict.get<NeoN::scalar>("scalarWriteAnsInt") == 2);
    REQUIRE(nfTestDict.get<NeoN::Vec3>("vector") == NeoN::Vec3(1.0, 2.0, 3.0));
    REQUIRE(nfTestDict.get<std::string>("word") == "word");

    NeoN::Dictionary& subNeoDict = nfTestDict.subDict("subDict");
    REQUIRE(subNeoDict.get<NeoN::scalar>("subScalar") == 4.1);
    REQUIRE(subNeoDict.get<NeoN::Vec3>("subVector") == NeoN::Vec3(5.0, 6.0, 7.0));
    REQUIRE(subNeoDict.get<std::string>("subWord") == "subWord");
}
