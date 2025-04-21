// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "common.hpp"

extern Foam::Time* timePtr; // A single time object

TEST_CASE("VolumeField")
{
    NeoN::Executor exec = GENERATE(
        NeoN::Executor(NeoN::CPUExecutor {}),
        NeoN::Executor(NeoN::SerialExecutor {}),
        NeoN::Executor(NeoN::GPUExecutor {})
    );
    std::string execName = std::visit([](auto e) { return e.name(); }, exec);

    Foam::Time& runTime = *timePtr;
    auto meshPtr = Foam::createMesh(exec, runTime);
    Foam::MeshAdapter& mesh = *meshPtr;
    auto nfMesh = mesh.nfMesh();

    auto ofT = randomScalarField(runTime, mesh, "T");
    auto ofU = randomVectorField(runTime, mesh, "U");

    SECTION("volumeScalarField " + execName)
    {
        auto nfT = constructFrom(exec, nfMesh, ofT);
        compare(nfT, ofT, ApproxScalar(1e-15));
    }

    SECTION("volumeVectorField " + execName)
    {
        auto nfU = constructFrom(exec, nfMesh, ofU);
        compare(nfU, ofU, ApproxVector(1e-15));
    }
}
