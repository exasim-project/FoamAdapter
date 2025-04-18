// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors


#include "FoamAdapter/finiteVolume/cellCentred/pressureVelocityCoupling/pressureVelocityCoupling.hpp"
#include "NeoN/finiteVolume/cellCentred/boundary.hpp"
#include "NeoN/finiteVolume/cellCentred/operators/gaussGreenGrad.hpp"
#include "Kokkos_Core.hpp"

namespace NeoN::finiteVolume::cellCentred
{

void constrainHbyA(
    VolumeField<Vec3>& HbyA,
    const VolumeField<Vec3>& U,
    const VolumeField<scalar>& p
)
{
    // const UnstructuredMesh& mesh = HbyA.mesh();
    const auto pIn = p.internalField().view();
    auto HbyAin = HbyA.internalField().view();
    auto [HbyABcValue, UBcValue] = spans(HbyA.boundaryField().value(), U.boundaryField().value());

    const std::vector<VolumeBoundary<Vec3>>& HbyABCs = HbyA.boundaryConditions();

    for (std::size_t patchi = 0; patchi < HbyABCs.size(); ++patchi)
    {
        parallelFor(
            HbyA.exec(),
            HbyA.boundaryField().range(patchi),
            KOKKOS_LAMBDA(const size_t bfacei) { HbyABcValue[bfacei] = UBcValue[bfacei]; }
        );
    }
}

std::tuple<VolumeField<scalar>, VolumeField<Vec3>>
discreteMomentumFields(const Expression<Vec3>& expr)
{
    const VolumeField<Vec3>& U = expr.getField();
    const UnstructuredMesh& mesh = U.mesh();
    const SparsityPattern& sparsityPattern = expr.sparsityPattern();
    const auto& ls = expr.linearSystem();
    const auto vol = mesh.cellVolumes().view();
    const auto values = ls.matrix().values().view();
    const auto rhs = ls.rhs().view();
    const auto diagOffset = sparsityPattern.diagOffset().view();
    const auto rowPtrs = ls.matrix().rowPtrs().view();

    auto rABCs = createExtrapolatedBCs<VolumeBoundary<scalar>>(mesh);
    VolumeField<scalar> rAU = VolumeField<scalar>(expr.exec(), "rAU", mesh, rABCs);

    rAU.internalField().apply(KOKKOS_LAMBDA(const size_t celli) {
        auto diagOffsetCelli = diagOffset[celli];
        // all the diagonal coefficients are the same
        return vol[celli] / (values[rowPtrs[celli] + diagOffsetCelli][0]);
    });

    auto OffDiagonalSourceBCs = createExtrapolatedBCs<VolumeBoundary<Vec3>>(mesh);
    VolumeField<Vec3> HbyA = VolumeField<Vec3>(expr.exec(), "HbyA", mesh, OffDiagonalSourceBCs);
    fill(HbyA.internalField(), zero<Vec3>());
    const std::size_t nInternalFaces = mesh.nInternalFaces();

    const auto exec = U.exec();
    const auto [owner, neighbour, surfFaceCells, ownOffs, neiOffs, internalU, internalRAU] = spans(
        mesh.faceOwner(),
        mesh.faceNeighbour(),
        mesh.boundaryMesh().faceCells(),
        sparsityPattern.ownerOffset(),
        sparsityPattern.neighbourOffset(),
        U.internalField(),
        rAU.internalField()
    );

    auto internalHbyA = HbyA.internalField().view();

    parallelFor(
        exec,
        {0, nInternalFaces},
        KOKKOS_LAMBDA(const size_t facei) {
            std::size_t own = static_cast<std::size_t>(owner[facei]);
            std::size_t nei = static_cast<std::size_t>(neighbour[facei]);

            std::size_t rowNeiStart = rowPtrs[nei];
            std::size_t rowOwnStart = rowPtrs[own];

            auto Lower = values[rowNeiStart + neiOffs[facei]];
            auto Upper = values[rowOwnStart + ownOffs[facei]];

            Kokkos::atomic_sub(&internalHbyA[nei], Lower[0] * internalU[own]);
            Kokkos::atomic_sub(&internalHbyA[own], Upper[0] * internalU[nei]);
        }
    );

    parallelFor(
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
    SurfaceField<scalar>& phi,
    const SurfaceField<scalar>& predictedPhi,
    const Expression<scalar>& expr
)
{
    const UnstructuredMesh& mesh = phi.mesh();
    const VolumeField<scalar>& p = expr.getField();
    const SparsityPattern sparsityPattern = expr.sparsityPattern();
    const std::size_t nInternalFaces = mesh.nInternalFaces();
    const auto exec = phi.exec();
    const auto [owner, neighbour, surfFaceCells, ownOffs, neiOffs, internalP] = spans(
        mesh.faceOwner(),
        mesh.faceNeighbour(),
        mesh.boundaryMesh().faceCells(),
        sparsityPattern.ownerOffset(),
        sparsityPattern.neighbourOffset(),
        p.internalField()
    );

    const auto& ls = expr.linearSystem();

    const auto rowPtrs = ls.matrix().rowPtrs().view();
    const auto colIdxs = ls.matrix().colIdxs().view();
    auto values = ls.matrix().values().view();
    auto rhs = ls.rhs().view();

    auto [iPhi, iPredPhi] = spans(phi.internalField(), predictedPhi.internalField());

    parallelFor(
        exec,
        {0, nInternalFaces},
        KOKKOS_LAMBDA(const size_t facei) {
            std::size_t own = static_cast<std::size_t>(owner[facei]);
            std::size_t nei = static_cast<std::size_t>(neighbour[facei]);

            std::size_t rowNeiStart = rowPtrs[nei];
            std::size_t rowOwnStart = rowPtrs[own];

            auto Upper = values[rowNeiStart + neiOffs[facei]];
            auto Lower = values[rowOwnStart + ownOffs[facei]];

            iPhi[facei] = iPredPhi[facei] - (Upper * internalP[nei] - Lower * internalP[own]);
        }
    );
}

void updateVelocity(
    VolumeField<Vec3>& U,
    const VolumeField<Vec3>& HbyA,
    VolumeField<scalar>& rAU,
    VolumeField<scalar>& p
)
{
    VolumeField<Vec3> gradP = GaussGreenGrad(p.exec(), p.mesh()).grad(p);
    auto [iHbyA, iRAU, iGradP] =
        spans(HbyA.internalField(), rAU.internalField(), gradP.internalField());

    U.internalField().apply(KOKKOS_LAMBDA(const std::size_t celli) {
        return iHbyA[celli] - iRAU[celli] * iGradP[celli];
    });
}

SurfaceField<scalar> flux(const VolumeField<Vec3>& volField)
{
    const auto exec = volField.exec();

    const UnstructuredMesh& mesh = volField.mesh();
    const std::size_t nInternalFaces = mesh.nInternalFaces();
    Input input = TokenList({std::string("linear")});
    auto linear = SurfaceInterpolation<Vec3>(exec, mesh, input);
    const SurfaceField<scalar> weight = linear.weight(volField);

    auto surfaceBCs = fvcc::createCalculatedBCs<fvcc::SurfaceBoundary<scalar>>(mesh);
    auto faceFlux = SurfaceField<scalar>(exec, "out", mesh, surfaceBCs);

    fill(faceFlux.internalField(), zero<scalar>());
    fill(faceFlux.boundaryField().value(), zero<scalar>());
    const auto [owner, neighbour, weightIn, faceAreas, volFieldIn, volFieldBc, bSf] = spans(
        mesh.faceOwner(),
        mesh.faceNeighbour(),
        weight.internalField(),
        mesh.faceAreas(),
        volField.internalField(),
        volField.boundaryField().value(),
        mesh.boundaryMesh().sf()
    );

    auto [faceFluxIn, bvalue] = spans(faceFlux.internalField(), faceFlux.boundaryField().value());

    parallelFor(
        exec,
        {0, nInternalFaces},
        KOKKOS_LAMBDA(const size_t facei) {
            std::size_t own = static_cast<std::size_t>(owner[facei]);
            std::size_t nei = static_cast<std::size_t>(neighbour[facei]);

            faceFluxIn[facei] =
                faceAreas[facei]
                & (weightIn[facei] * (volFieldIn[own] - volFieldIn[nei]) + volFieldIn[nei]);
        }
    );

    parallelFor(
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
