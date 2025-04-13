// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2023 NeoN authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "NeoN/NeoN.hpp"
#include "benchmarks/catch_main.hpp"
#include "NeoN/test/catch2/executorGenerator.hpp"

namespace fvcc = NeoN::finiteVolume::cellCentred;
namespace dsl = NeoN::dsl;

#include "fv.H"
#include "gaussConvectionScheme.H"
#include "NeoN/core/input.hpp"
#include "NeoN/dsl/explicit.hpp"

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

    forAll(t, celli)
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

TEST_CASE("DivOperator")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    SECTION("OpenFOAM")
    {
        std::unique_ptr<Foam::fvMesh> meshPtr = Foam::createMesh(runTime);
        Foam::fvMesh& mesh = *meshPtr;

        auto ofT = randomScalarField(runTime, mesh, "T");
        Foam::surfaceScalarField ofPhi(
            Foam::IOobject(
                "phi",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("phi", Foam::dimless, 0.0)
        );
        forAll(ofPhi, facei)
        {
            ofPhi[facei] = facei;
        }
        SECTION("with Allocation")
        {
            BENCHMARK(std::string("OpenFOAM"))
            {
                Foam::IStringStream is("linear");
                Foam::fv::gaussConvectionScheme<Foam::scalar> foamDivScalar(mesh, ofPhi, is);
                Foam::volScalarField ofDivT("ofDivT", foamDivScalar.fvcDiv(ofPhi, ofT));
                return;
            };
        }
    }


    SECTION("NeoN")
    {
        auto [execName, exec] = GENERATE(allAvailableExecutor());

        std::unique_ptr<Foam::MeshAdapter> meshPtr = Foam::createMesh(exec, runTime);
        Foam::MeshAdapter& mesh = *meshPtr;
        const auto& nfMesh = mesh.nfMesh();
        // linear interpolation hardcoded for now


        auto ofT = randomScalarField(runTime, mesh, "T");
        auto nfT = constructFrom(exec, nfMesh, ofT);
        nfT.correctBoundaryConditions();

        Foam::surfaceScalarField ofPhi(
            Foam::IOobject(
                "phi",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("phi", Foam::dimless, 0.0)
        );
        forAll(ofPhi, facei)
        {
            ofPhi[facei] = facei;
        }

        auto nfPhi = constructSurfaceField(exec, nfMesh, ofPhi);

        SECTION("with Allocation")
        {
            NeoN::TokenList scheme({std::string("linear")});

            BENCHMARK(std::string(execName))
            {
                fvcc::VolumeField<NeoN::scalar> divPhiT =
                    fvcc::GaussGreenDiv<NeoN::scalar>(exec, nfMesh, scheme)
                        .div(nfPhi, nfT, dsl::Coeff(1.0));
                return;
            };
        }

        SECTION("No allocation")
        {
            auto nfDivT = constructFrom(exec, nfMesh, ofT);
            NeoN::TokenList scheme({std::string("linear")});

            BENCHMARK(std::string(execName))
            {
                NeoN::fill(nfDivT.internalField(), 0.0);
                NeoN::fill(nfDivT.boundaryField().value(), 0.0);
                fvcc::GaussGreenDiv<NeoN::scalar>(exec, nfMesh, scheme)
                    .div(nfDivT, nfPhi, nfT, dsl::Coeff(1.0));
                Kokkos::fence();
                return;
            };
        }

        // SECTION("compute div from dsl::exp")
        // {
        //     NeoN::TokenList scheme = NeoN::TokenList({std::string("Gauss"),
        //     std::string("linear")});

        //     auto nfDivT = constructFrom(exec, nfMesh, ofDivT);
        //     NeoN::fill(nfDivT.internalField(), 0.0);
        //     NeoN::fill(nfDivT.boundaryField().value(), 0.0);
        //     dsl::SpatialOperator divOp = dsl::exp::div(nfPhi, nfT);
        //     divOp.build(scheme);
        //     divOp.explicitOperation(nfDivT.internalField());

        //     nfDivT.correctBoundaryConditions();

        //     compare(nfT, ofT, ApproxScalar(1e-15), false);

        //     compare(nfDivT, ofDivT, ApproxScalar(1e-15), false);
        // }
    }
}

// TEST_CASE("Field<scalar>::addition no copy", "[bench]")
// {
//     auto size = GENERATE(1 << 16, 1 << 17, 1 << 18, 1 << 19, 1 << 20);
//     auto [execName, exec] = GENERATE(allAvailableExecutor());

//     DYNAMIC_SECTION("" << size)
//     {
//         NeoN::Field<NeoN::scalar> cpuA(exec, size);
//         const auto viewA = cpuA.view();
//         NeoN::fill(cpuA, 1.0);
//         NeoN::Field<NeoN::scalar> cpuB(exec, size);
//         NeoN::fill(cpuB, 2.0);
//         const auto viewB = cpuB.view();
//         NeoN::Field<NeoN::scalar> cpuC(exec, size);
//         NeoN::fill(cpuC, 0.0);

//         BENCHMARK(std::string(execName)) {
//             NeoN::parallelFor(cpuC, KOKKOS_LAMBDA(const int i) {
//                 return viewB[i] + viewB[i];
//             });
//             Kokkos::fence();
//             return;
//             // return (cpuC = cpuA + cpuB);
//         };
//     }
// }

// TEST_CASE("Field<scalar>::multiplication", "[bench]")
// {
//     auto size = GENERATE(1 << 16, 1 << 17, 1 << 18, 1 << 19, 1 << 20);

//     auto [execName, exec] = GENERATE(allAvailableExecutor());

//     DYNAMIC_SECTION("" << size)
//     {
//         NeoN::Field<NeoN::scalar> cpuA(exec, size);
//         NeoN::fill(cpuA, 1.0);
//         NeoN::Field<NeoN::scalar> cpuB(exec, size);
//         NeoN::fill(cpuB, 2.0);
//         NeoN::Field<NeoN::scalar> cpuC(exec, size);
//         NeoN::fill(cpuC, 0.0);

//         BENCHMARK(std::string(execName)) { return (cpuC = cpuA * cpuB); };
//     }
// }
