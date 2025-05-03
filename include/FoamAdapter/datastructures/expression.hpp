// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 FoamAdapter authors
// TODO: move to cellCenred dsl?

#pragma once

#include "NeoN/NeoN.hpp"

namespace dsl = NeoN::dsl;

namespace FoamAdapter
{

/*@brief extends expression by giving access to assembled matrix
 * @note used in neoIcoFOAM directly instead of dsl::expression
 * TODO: implement flag if matrix is assembled or not -> if not assembled call assemble
 * for dependent operations like discrete momentum fields
 * needs storage for assembled matrix? and whether update is needed like for rAU and HbyA
 */
template<typename ValueType, typename IndexType = NeoN::localIdx>
class Expression
{
public:

    Expression(
        dsl::Expression<ValueType> expr,
        NeoN::finiteVolume::cellCentred::VolumeField<ValueType>& psi,
        const NeoN::Dictionary& fvSchemes,
        const NeoN::Dictionary& fvSolution
    )
        : psi_(psi)
        , expr_(expr)
        , fvSchemes_(fvSchemes)
        , fvSolution_(fvSolution)
        , sparsityPattern_(
              NeoN::finiteVolume::cellCentred::SparsityPattern::readOrCreate(psi.mesh())
          )
        , ls_(NeoN::la::createEmptyLinearSystem<
              ValueType,
              NeoN::localIdx,
              NeoN::finiteVolume::cellCentred::SparsityPattern>(*sparsityPattern_.get()))
    {
        expr_.read(fvSchemes_);
        // assemble();
    };

    Expression(const Expression& ls)
        : psi_(ls.psi_)
        , expr_(ls.expr_)
        , fvSchemes_(ls.fvSchemes_)
        , fvSolution_(ls.fvSolution_)
        , ls_(ls.ls_)
        , sparsityPattern_(ls.sparsityPattern_) {};

    ~Expression() = default;

    [[nodiscard]] NeoN::la::LinearSystem<ValueType, IndexType>& linearSystem() { return ls_; }
    [[nodiscard]] NeoN::finiteVolume::cellCentred::SparsityPattern& sparsityPattern()
    {
        if (!sparsityPattern_)
        {
            NF_THROW(std::string("fvcc:LinearSystem:sparsityPattern: sparsityPattern is null"));
        }
        return *sparsityPattern_;
    }

    NeoN::finiteVolume::cellCentred::VolumeField<ValueType>& getVector() { return this->psi_; }

    const NeoN::finiteVolume::cellCentred::VolumeField<ValueType>& getVector() const
    {
        return this->psi_;
    }

    [[nodiscard]] const NeoN::la::LinearSystem<ValueType, IndexType>& linearSystem() const
    {
        return ls_;
    }
    [[nodiscard]] const NeoN::finiteVolume::cellCentred::SparsityPattern& sparsityPattern() const
    {
        if (!sparsityPattern_)
        {
            NF_THROW("fvcc:LinearSystem:sparsityPattern: sparsityPattern is null");
        }
        return *sparsityPattern_;
    }

    const NeoN::Executor& exec() const { return ls_.exec(); }


    void assemble(NeoN::scalar t, NeoN::scalar dt)
    {
        auto vol = psi_.mesh().cellVolumes().view();
        auto expSource = expr_.explicitOperation(psi_.mesh().nCells());
        expr_.explicitOperation(expSource, t, dt);
        auto expSourceView = expSource.view();
        fill(ls_.rhs(), NeoN::zero<ValueType>());
        fill(ls_.matrix().values(), NeoN::zero<ValueType>());
        expr_.implicitOperation(ls_);
        // TODO rename implicitOperation -> assembleLinearSystem
        expr_.implicitOperation(ls_, t, dt);
        auto rhs = ls_.rhs().view();
        // we subtract the explicit source term from the rhs
        NeoN::parallelFor(
            exec(),
            {0, rhs.size()},
            KOKKOS_LAMBDA(const NeoN::localIdx i) { rhs[i] -= expSourceView[i] * vol[i]; }
        );
    }

    void assemble()
    {
        if (expr_.temporalOperators().size() == 0 && expr_.spatialOperators().size() == 0)
        {
            NF_ERROR_EXIT("No temporal or implicit terms to solve.");
        }

        if (expr_.temporalOperators().size() > 0)
        {
            // integrate equations in time
            // NeoN::timeIntegration::TimeIntegration<VolumeField<ValueType>> timeIntegrator(
            //     fvSchemes_.subDict("ddtSchemes"), fvSolution_
            // );
            // timeIntegrator.solve(expr_, psi_, t, dt);
        }
        else
        {
            // solve sparse matrix system
            auto vol = psi_.mesh().cellVolumes().view();
            auto expSource = expr_.explicitOperation(psi_.mesh().nCells());
            auto expSourceView = expSource.view();

            ls_ = expr_.implicitOperation();
            auto rhs = ls_.rhs().view();
            // we subtract the explicit source term from the rhs
            NeoN::parallelFor(
                exec(),
                {0, rhs.size()},
                KOKKOS_LAMBDA(const NeoN::localIdx i) { rhs[i] -= expSourceView[i] * vol[i]; }
            );
        }
    }

    // TODO unify with dsl/solver.hpp
    void solve(NeoN::scalar, NeoN::scalar)
    {
        // dsl::solve(expr_, psi_, t, dt, fvSchemes_, fvSolution_);
        if (expr_.temporalOperators().size() == 0 && expr_.spatialOperators().size() == 0)
        {
            NF_ERROR_EXIT("No temporal or implicit terms to solve.");
        }
        if (expr_.temporalOperators().size() > 0)
        {
            NF_ERROR_EXIT("Not implemented");
            //     // integrate equations in time
            //     NeoN::timeIntegration::TimeIntegration<VolumeField<ValueType>> timeIntegrator(
            //         fvSchemes_.subDict("ddtSchemes"), fvSolution_
            //     );
            //     timeIntegrator.solve(expr_, psi_, t, dt);
        }
        else
        {
            auto exec = psi_.exec();
            auto solver = NeoN::la::Solver(exec, fvSolution_);
            solver.solve(ls_, psi_.internalVector());
            // NF_ERROR_EXIT("No linear solver is available, build with -DNeoN_WITH_GINKGO=ON");
        }
    }

    void setReference(const IndexType refCell, ValueType refValue)
    {
        // TODO currently assumes that matrix is already assembled
        const auto diagOffset = sparsityPattern_->diagOffset().view();
        const auto rowOffs = ls_.matrix().rowOffs().view();
        auto rhs = ls_.rhs().view();
        auto values = ls_.matrix().values().view();
        NeoN::parallelFor(
            ls_.exec(),
            {refCell, refCell + 1},
            KOKKOS_LAMBDA(const std::size_t refCelli) {
                auto diagIdx = rowOffs[refCelli] + diagOffset[refCelli];
                auto diagValue = values[diagIdx];
                rhs[refCelli] += diagValue * refValue;
                values[diagIdx] += diagValue;
            }
        );
    }

private:

    NeoN::finiteVolume::cellCentred::VolumeField<ValueType>& psi_;
    NeoN::dsl::Expression<ValueType> expr_;
    const NeoN::Dictionary& fvSchemes_;
    const NeoN::Dictionary& fvSolution_;
    std::shared_ptr<NeoN::finiteVolume::cellCentred::SparsityPattern> sparsityPattern_;
    NeoN::la::LinearSystem<ValueType, IndexType> ls_;
};

template<typename ValueType, typename IndexType = NeoN::localIdx>
NeoN::Vector<ValueType> diag(
    const la::LinearSystem<ValueType, IndexType>& ls,
    const NeoN::finiteVolume::cellCentred::SparsityPattern& sparsityPattern
)
{
    NeoN::Vector<ValueType> diagonal(ls.exec(), sparsityPattern.diagOffset().size(), 0.0);
    auto diagView = diagonal.view();

    const auto diagOffset = sparsityPattern.diagOffset().view();
    const auto [matrix, b] = ls.view();
    NeoN::parallelFor(
        ls.exec(),
        {0, diagOffset.size()},
        KOKKOS_LAMBDA(const std::size_t celli) {
            auto diagOffsetCelli = diagOffset[celli];
            diagView[celli] = matrix.values[matrix.rowOffs[celli] + diagOffsetCelli];
        }
    );
    return diagonal;
}


template<typename ValueType, typename IndexType = NeoN::localIdx>
NeoN::finiteVolume::cellCentred::VolumeField<ValueType> applyOperator(
    const la::LinearSystem<ValueType, IndexType>& ls,
    const NeoN::finiteVolume::cellCentred::VolumeField<ValueType>& psi
)
{
    NeoN::finiteVolume::cellCentred::VolumeField<ValueType> resultVector(
        psi.exec(),
        "ls_" + psi.name,
        psi.mesh(),
        psi.internalVector(),
        psi.boundaryData(),
        psi.boundaryConditions()
    );

    auto [result, x] = views(resultVector.internalVector(), psi.internalVector());
    const auto [matrix, b] = ls.view();

    NeoN::parallelFor(
        resultVector.exec(),
        {0, result.size()},
        KOKKOS_LAMBDA(const std::size_t rowi) {
            IndexType rowStart = matrix.rowOffs[rowi];
            IndexType rowEnd = matrix.rowOffs[rowi + 1];
            ValueType sum = 0.0;
            for (IndexType coli = rowStart; coli < rowEnd; coli++)
            {
                sum += matrix.values[coli] * x[matrix.colIdxs[coli]];
            }
            result[rowi] = sum - b[rowi];
        }
    );

    return resultVector;
}


template<typename ValueType, typename IndexType = NeoN::localIdx>
NeoN::finiteVolume::cellCentred::VolumeField<ValueType> operator&(
    const Expression<ValueType, IndexType> expr,
    const NeoN::finiteVolume::cellCentred::VolumeField<ValueType>& psi
)
{
    return applyOperator(expr.linearSystem(), psi);
}

}
