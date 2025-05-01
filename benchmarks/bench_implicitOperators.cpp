// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2023 NeoN authors

#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "NeoN/NeoN.hpp"
#include "benchmarks/catch_main.hpp"
#include "test/catch2/executorGenerator.hpp"
#include "common.hpp"

namespace fvcc = NeoN::finiteVolume::cellCentred;
namespace dsl = NeoN::dsl;

#include "fv.H"
#include "fvc.H"
#include "gaussGrad.H"
#include "gaussConvectionScheme.H"
#include "gaussLaplacianScheme.H"


TEST_CASE("DivOperator")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    SECTION("OpenFOAM")
    {
        std::unique_ptr<Foam::fvMesh> meshPtr = FoamAdapter::createMesh(runTime);
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
        for (auto facei = 0; facei < ofPhi.size(); facei++)
        {
            ofPhi[facei] = facei;
        }
        SECTION("with Allocation")
        {
            BENCHMARK(std::string("OpenFOAM"))
            {
                Foam::IStringStream is("linear");
                Foam::fv::gaussConvectionScheme<Foam::scalar> foamDivScalar(mesh, ofPhi, is);
                auto fvmDivT = foamDivScalar.fvmDiv(ofPhi, ofT);
                return;
            };
        }
    }


    SECTION("NeoN")
    {
        auto [execName, exec] = GENERATE(allAvailableExecutor());

        std::unique_ptr<FoamAdapter::MeshAdapter> meshPtr = FoamAdapter::createMesh(exec, runTime);
        FoamAdapter::MeshAdapter& mesh = *meshPtr;
        const auto& nfMesh = mesh.nfMesh();

        // linear interpolation hardcoded for now
        auto ofT = randomScalarField(runTime, mesh, "T");
        auto nfT = FoamAdapter::constructFrom(exec, nfMesh, ofT);
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
        for (auto facei = 0; facei < ofPhi.size(); facei++)
        {
            ofPhi[facei] = facei;
        }

        auto nfPhi = FoamAdapter::constructSurfaceField(exec, nfMesh, ofPhi);

        SECTION("with Allocation")
        {
            NeoN::TokenList scheme({std::string("linear")});

            BENCHMARK(std::string(execName))
            {
                NeoN::la::LinearSystem<NeoN::scalar, NeoN::localIdx> ls(
                    la::createEmptyLinearSystem<
                        NeoN::scalar,
                        NeoN::localIdx,
                        fvcc::SparsityPattern>(*fvcc::SparsityPattern::readOrCreate(nfMesh).get())
                );
                fvcc::GaussGreenDiv<NeoN::scalar>(exec, nfMesh, scheme)
                    .div(ls, nfPhi, nfT, NeoN::dsl::Coeff(1.0));
                Kokkos::fence();
                return;
            };
        }

        SECTION("No allocation")
        {
            NeoN::la::LinearSystem<NeoN::scalar, NeoN::localIdx> ls(
                la::createEmptyLinearSystem<NeoN::scalar, NeoN::localIdx, fvcc::SparsityPattern>(
                    *fvcc::SparsityPattern::readOrCreate(nfMesh).get()
                )
            );
            NeoN::TokenList scheme({std::string("linear")});

            BENCHMARK(std::string(execName))
            {
                NeoN::fill(ls.matrix().values(), 0.0);
                NeoN::fill(ls.rhs(), 0.0);
                fvcc::GaussGreenDiv<NeoN::scalar>(exec, nfMesh, scheme)
                    .div(ls, nfPhi, nfT, dsl::Coeff(1.0));
                Kokkos::fence();
                return;
            };
        }

        // TODO dsl
    }
}


TEST_CASE("LaplacianOperator")
{
    Foam::Time& runTime = *timePtr;
    Foam::argList& args = *argsPtr;

    SECTION("OpenFOAM")
    {
        std::unique_ptr<Foam::fvMesh> meshPtr = FoamAdapter::createMesh(runTime);
        Foam::fvMesh& mesh = *meshPtr;

        auto ofT = randomScalarField(runTime, mesh, "T");
        Foam::surfaceScalarField ofGamma(
            Foam::IOobject(
                "Gamma",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("Gamma", Foam::dimless, 1.0)
        );

        SECTION("with Allocation")
        {
            BENCHMARK(std::string("OpenFOAM"))
            {
                Foam::IStringStream is("linear uncorrected");
                Foam::fv::gaussLaplacianScheme<Foam::scalar, Foam::scalar> foamLapScalar(mesh, is);
                auto fvmLapT = foamLapScalar.fvmLaplacian(ofGamma, ofT);
                return;
            };
        }
    }


    SECTION("NeoN")
    {
        auto [execName, exec] = GENERATE(allAvailableExecutor());

        std::unique_ptr<FoamAdapter::MeshAdapter> meshPtr = FoamAdapter::createMesh(exec, runTime);
        FoamAdapter::MeshAdapter& mesh = *meshPtr;
        const auto& nfMesh = mesh.nfMesh();
        // linear interpolation hardcoded for now

        auto ofT = randomScalarField(runTime, mesh, "T");
        auto nfT = FoamAdapter::constructFrom(exec, nfMesh, ofT);
        nfT.correctBoundaryConditions();

        Foam::surfaceScalarField ofGamma(
            Foam::IOobject(
                "Gamma",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("Gamma", Foam::dimless, 1.0)
        );

        auto nfGamma = FoamAdapter::constructSurfaceField(exec, nfMesh, ofGamma);

        SECTION("with Allocation")
        {
            NeoN::TokenList scheme({std::string("linear"), std::string("uncorrected")});

            BENCHMARK(std::string(execName))
            {
                NeoN::la::LinearSystem<NeoN::scalar, NeoN::localIdx> ls(
                    la::createEmptyLinearSystem<
                        NeoN::scalar,
                        NeoN::localIdx,
                        fvcc::SparsityPattern>(*fvcc::SparsityPattern::readOrCreate(nfMesh).get())
                );
                fvcc::GaussGreenLaplacian<NeoN::scalar>(exec, nfMesh, scheme)
                    .laplacian(ls, nfGamma, nfT, dsl::Coeff(1.0));
                Kokkos::fence();
                return;
            };
        }

        SECTION("No allocation")
        {
            NeoN::la::LinearSystem<NeoN::scalar, NeoN::localIdx> ls(
                la::createEmptyLinearSystem<NeoN::scalar, NeoN::localIdx, fvcc::SparsityPattern>(
                    *fvcc::SparsityPattern::readOrCreate(nfMesh).get()
                )
            );
            NeoN::TokenList scheme({std::string("linear"), std::string("uncorrected")});

            BENCHMARK(std::string(execName))
            {
                NeoN::fill(ls.matrix().values(), 0.0);
                NeoN::fill(ls.rhs(), 0.0);
                fvcc::GaussGreenLaplacian<NeoN::scalar>(exec, nfMesh, scheme)
                    .laplacian(ls, nfGamma, nfT, dsl::Coeff(1.0));
                Kokkos::fence();
                return;
            };
        }

        // TODO dsl
    }
}
