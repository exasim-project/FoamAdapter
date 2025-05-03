// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 FoamAdapter authors

#include "NeoN/NeoN.hpp"

#include "FoamAdapter/algorithms/pressureVelocityCoupling.hpp"
#include "Kokkos_Core.hpp"

namespace FoamAdapter
{

void constrainHbyA(
    nnfvcc::VolumeField<Vec3>& HbyA,
    const nnfvcc::VolumeField<Vec3>& U,
    const nnfvcc::VolumeField<scalar>& p
)
{
    // const UnstructuredMesh& mesh = HbyA.mesh();
    const auto pIn = p.internalVector().view();
    auto HbyAin = HbyA.internalVector().view();
    auto [HbyABcValue, UBcValue] = views(HbyA.boundaryData().value(), U.boundaryData().value());

    const auto& UBCs = U.boundaryConditions();

    for (auto patchi = 0; patchi < UBCs.size(); ++patchi)
    {
        bool assignable = UBCs[patchi].attributes().get<bool>("assignable");
        if (!assignable)
        {
            parallelFor(
                HbyA.exec(),
                HbyA.boundaryData().range(patchi),
                KOKKOS_LAMBDA(const size_t bfacei) { HbyABcValue[bfacei] = UBcValue[bfacei]; }
            );
        }
    }
}

std::tuple<nnfvcc::VolumeField<scalar>, nnfvcc::VolumeField<Vec3>>
discreteMomentumFields(const Expression<Vec3>& expr)
{
    const nnfvcc::VolumeField<Vec3>& U = expr.getVector();
    const auto& mesh = U.mesh();
    const auto& sparsityPattern = expr.sparsityPattern();
    const auto& ls = expr.linearSystem();
    const auto vol = mesh.cellVolumes().view();
    const auto values = ls.matrix().values().view();
    const auto rhs = ls.rhs().view();
    const auto diagOffset = sparsityPattern.diagOffset().view();
    const auto rowPtrs = ls.matrix().rowOffs().view();

    auto rABCs = nnfvcc::createExtrapolatedBCs<nnfvcc::VolumeBoundary<scalar>>(mesh);
    auto rAU = nnfvcc::VolumeField<scalar>(expr.exec(), "rAU", mesh, rABCs);

    rAU.internalVector().apply(KOKKOS_LAMBDA(const size_t celli) {
        auto diagOffsetCelli = diagOffset[celli];
        // all the diagonal coefficients are the same
        return vol[celli] / (values[rowPtrs[celli] + diagOffsetCelli][0]);
    });

    auto OffDiagonalSourceBCs = nnfvcc::createExtrapolatedBCs<nnfvcc::VolumeBoundary<Vec3>>(mesh);
    auto HbyA = nnfvcc::VolumeField<Vec3>(expr.exec(), "HbyA", mesh, OffDiagonalSourceBCs);
    NeoN::fill(HbyA.internalVector(), NeoN::zero<Vec3>());
    const auto nInternalFaces = mesh.nInternalFaces();

    const auto exec = U.exec();
    const auto [owner, neighbour, surfFaceCells, ownOffs, neiOffs, internalU, internalRAU] = views(
        mesh.faceOwner(),
        mesh.faceNeighbour(),
        mesh.boundaryMesh().faceCells(),
        sparsityPattern.ownerOffset(),
        sparsityPattern.neighbourOffset(),
        U.internalVector(),
        rAU.internalVector()
    );

    auto internalHbyA = HbyA.internalVector().view();

    NeoN::parallelFor(
        exec,
        {0, nInternalFaces},
        KOKKOS_LAMBDA(const size_t facei) {
            auto own = owner[facei];
            auto nei = neighbour[facei];

            auto rowNeiStart = rowPtrs[nei];
            auto rowOwnStart = rowPtrs[own];

            auto Lower = values[rowNeiStart + neiOffs[facei]];
            auto Upper = values[rowOwnStart + ownOffs[facei]];

            Kokkos::atomic_sub(&internalHbyA[nei], Lower[0] * internalU[own]);
            Kokkos::atomic_sub(&internalHbyA[own], Upper[0] * internalU[nei]);
        }
    );

    NeoN::parallelFor(
        exec,
        {0, internalHbyA.size()},
        KOKKOS_LAMBDA(const size_t celli) {
            internalHbyA[celli] += rhs[celli];
            internalHbyA[celli] *= internalRAU[celli] / vol[celli];
        }
    );

    HbyA.correctBoundaryConditions();
    rAU.correctBoundaryConditions();

    return {rAU, HbyA};
}


void updateFaceVelocity(
    nnfvcc::SurfaceField<scalar>& phi,
    const nnfvcc::SurfaceField<scalar>& predictedPhi,
    const Expression<scalar>& expr
)
{
    const auto& mesh = phi.mesh();
    const auto& p = expr.getVector();
    const auto sparsityPattern = expr.sparsityPattern();
    const auto nInternalFaces = mesh.nInternalFaces();
    const auto exec = phi.exec();
    const auto [owner, neighbour, surfFaceCells, ownOffs, neiOffs, internalP] = views(
        mesh.faceOwner(),
        mesh.faceNeighbour(),
        mesh.boundaryMesh().faceCells(),
        sparsityPattern.ownerOffset(),
        sparsityPattern.neighbourOffset(),
        p.internalVector()
    );

    const auto& ls = expr.linearSystem();

    const auto rowPtrs = ls.matrix().rowOffs().view();
    const auto colIdxs = ls.matrix().colIdxs().view();
    auto values = ls.matrix().values().view();
    auto rhs = ls.rhs().view();

    auto [iPhi, iPredPhi] = views(phi.internalVector(), predictedPhi.internalVector());

    NeoN::parallelFor(
        exec,
        {0, nInternalFaces},
        KOKKOS_LAMBDA(const size_t facei) {
            auto own = static_cast<std::size_t>(owner[facei]);
            auto nei = static_cast<std::size_t>(neighbour[facei]);

            auto rowNeiStart = rowPtrs[nei];
            auto rowOwnStart = rowPtrs[own];

            auto Upper = values[rowNeiStart + neiOffs[facei]];
            auto Lower = values[rowOwnStart + ownOffs[facei]];

            iPhi[facei] = iPredPhi[facei] - (Upper * internalP[nei] - Lower * internalP[own]);
        }
    );
}

void updateVelocity(
    nnfvcc::VolumeField<Vec3>& U,
    const nnfvcc::VolumeField<Vec3>& HbyA,
    nnfvcc::VolumeField<scalar>& rAU,
    nnfvcc::VolumeField<scalar>& p
)
{
    nnfvcc::VolumeField<Vec3> gradP = nnfvcc::GaussGreenGrad(p.exec(), p.mesh()).grad(p);
    auto [iHbyA, iRAU, iGradP] =
        views(HbyA.internalVector(), rAU.internalVector(), gradP.internalVector());

    U.internalVector().apply(KOKKOS_LAMBDA(const std::size_t celli) {
        return iHbyA[celli] - iRAU[celli] * iGradP[celli];
    });
}

nnfvcc::SurfaceField<scalar> flux(const nnfvcc::VolumeField<Vec3>& volField)
{
    const auto exec = volField.exec();

    const auto& mesh = volField.mesh();
    const auto nInternalFaces = mesh.nInternalFaces();
    NeoN::Input input = NeoN::TokenList({std::string("linear")});
    auto linear = nnfvcc::SurfaceInterpolation<Vec3>(exec, mesh, input);
    const auto weight = linear.weight(volField);

    auto surfaceBCs = nnfvcc::createCalculatedBCs<nnfvcc::SurfaceBoundary<scalar>>(mesh);
    auto faceFlux = nnfvcc::SurfaceField<scalar>(exec, "out", mesh, surfaceBCs);

    NeoN::fill(faceFlux.internalVector(), NeoN::zero<scalar>());
    NeoN::fill(faceFlux.boundaryData().value(), NeoN::zero<scalar>());
    const auto [owner, neighbour, weightIn, faceAreas, volFieldIn, volFieldBc, bSf] = views(
        mesh.faceOwner(),
        mesh.faceNeighbour(),
        weight.internalVector(),
        mesh.faceAreas(),
        volField.internalVector(),
        volField.boundaryData().value(),
        mesh.boundaryMesh().sf()
    );

    auto [faceFluxIn, bvalue] = views(faceFlux.internalVector(), faceFlux.boundaryData().value());

    NeoN::parallelFor(
        exec,
        {0, nInternalFaces},
        KOKKOS_LAMBDA(const size_t facei) {
            auto own = static_cast<std::size_t>(owner[facei]);
            auto nei = static_cast<std::size_t>(neighbour[facei]);

            faceFluxIn[facei] =
                faceAreas[facei]
                & (weightIn[facei] * (volFieldIn[own] - volFieldIn[nei]) + volFieldIn[nei]);
        }
    );

    NeoN::parallelFor(
        exec,
        {nInternalFaces, faceFluxIn.size()},
        KOKKOS_LAMBDA(const size_t facei) {
            auto faceBCI = facei - nInternalFaces;

            faceFluxIn[facei] = bSf[faceBCI] & volFieldBc[faceBCI];
            bvalue[faceBCI] = bSf[faceBCI] & volFieldBc[faceBCI];
        }
    );

    return faceFlux;
}

}
